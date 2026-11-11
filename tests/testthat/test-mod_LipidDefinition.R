#------------------------------------------------------------------------------
# Fixtures
#------------------------------------------------------------------------------

make_chain <- function(...) {
  chain <- list(
    chain_id = "AAAAA",
    C_min = 16, C_max = 18, even_steps = TRUE,
    DB_min = 0, DB_max = 3,
    priority = "basic",
    modifications = list()
  )
  overrides <- list(...)
  chain[names(overrides)] <- overrides  
  chain
}



make_validation <- function(...) {
  vals <- modifyList(
    list(
      name = "Validation name",
      adjustments = "-NH3",
      val_priority = "basic",
      selectedComponents = list("Backbone", "AAAAA")
    ),
    list(...)
  )
  vals
}

make_settings <- function(n_chains = 1) {
  list(
    head_group = "C6H12O6",
    shorthand = "GLC",
    chains = lapply(seq_len(n_chains), function(i) {
      make_chain(chain_id = paste0("CH", i, "XYZ"))
    }),
    validations = list(make_validation())
  )
}

#------------------------------------------------------------------------------
# Test functions
#------------------------------------------------------------------------------

test_that("create_unique_ids returns unique ids of the requested length", {
  ids <- create_unique_ids(n = 100, char_len = 5)
  expect_length(ids, 100)
  expect_equal(length(unique(ids)), 100)
  expect_true(all(nchar(ids) == 5))
})

test_that("build_chain_svg draws one circle per carbon", {
  svg <- build_chain_svg(c_min = 10, c_max = 20, db_min = 0, db_max = 2)
  n_circles <- lengths(regmatches(svg, gregexpr("<circle", svg)))
  expect_equal(n_circles, 20)
})

test_that("build_chain_svg survives edge cases", {
  expect_no_error(build_chain_svg(1, 1, 0, 0))
  expect_no_error(build_chain_svg(10, 10, 0, 0))   # c_min == c_max
  expect_no_error(build_chain_svg(2, 2, 0, 1))     # smallest chain with a double bond
})

#------------------------------------------------------------------------------
# Test chain block server
#------------------------------------------------------------------------------

chain_args <- function(init_chain = make_chain(), on_remove = function(id) NULL) {
  list(
    id = "chain_1",
    module_id = init_chain$chain_id,
    init_chain = init_chain,
    on_remove = on_remove
  )
}

test_that("chain starts with the values it was initialised with", {
  testServer(chain_block_server, args = chain_args(), {
    session$flushReact()
    out <- session$returned()
    expect_equal(out$chain_id, "AAAAA")
    expect_equal(out$C_min, 16)
    expect_equal(out$C_max, 18)
    expect_true(out$even_steps)
    expect_equal(out$DB_max, 3)
  })
})

test_that("c_max cannot go below c_min", {
  testServer(chain_block_server, args = chain_args(), {
    session$flushReact()
    session$setInputs(c_min = 12)
    session$setInputs(c_max = 8)
    expect_equal(session$returned()$C_max, 12)
  })
})

test_that("c_min and c_max are clamped to 0..100", {
  testServer(chain_block_server, args = chain_args(), {
    session$flushReact()
    session$setInputs(c_min = -5)
    expect_equal(chain_state$C_min, 0)
    session$setInputs(c_max = 500)
    expect_equal(chain_state$C_max, 100)
  })
})

test_that("double bonds are capped at floor(C_max / 2)", {
  testServer(chain_block_server, args = chain_args(), {
    session$flushReact()
    session$setInputs(c_max = 20)
    session$setInputs(db_max = 15)
    expect_equal(chain_state$DB_max, 10)
  })
})

test_that("db_max cannot go below db_min", {
  testServer(chain_block_server, args = chain_args(), {
    session$flushReact()
    session$setInputs(db_min = 2)
    session$setInputs(db_max = 0)
    expect_equal(session$returned()$DB_max, 2)
  })
})

test_that("even_steps toggle is reflected in the state", {
  testServer(chain_block_server, args = chain_args(), {
    session$flushReact()
    session$setInputs(even_steps = FALSE)
    expect_false(session$returned()$even_steps)
    session$setInputs(even_steps = TRUE)
    expect_true(session$returned()$even_steps)
  })
})

test_that("priority changes are reflected in the state", {
  testServer(chain_block_server, args = chain_args(), {
    session$flushReact()
    session$setInputs(chain_priority = "important")
    expect_equal(session$returned()$priority, "important")
  })
})

test_that("modifications can be added and edited", {
  testServer(chain_block_server, args = chain_args(), {
    session$flushReact()
    session$setInputs(add_mod = 1)
    mods <- session$returned()$modifications
    expect_length(mods, 1)
    
    session$setInputs(mods_1_formula = "OH")
    session$setInputs(mods_1_mod_optional = FALSE)
    
    mod <- session$returned()$modifications[[1]]
    expect_equal(mod$formula, "OH")
    expect_false(mod$optional)
  })
})

test_that("modifications from the init chain are loaded", {
  init <- make_chain(modifications = list(list(formula = "OH", optional = TRUE)))
  testServer(chain_block_server, args = chain_args(init_chain = init), {
    session$flushReact()
    mods <- session$returned()$modifications
    expect_length(mods, 1)
    expect_equal(mods[[1]]$formula, "OH")
  })
})

#------------------------------------------------------------------------------
# mod_LipidDefinition_server: coordination
#------------------------------------------------------------------------------

test_that("chains can be added and each gets a unique id", {
  testServer(mod_LipidDefinition_server,
             args = list(load_trigger = reactiveVal(NULL)), {
               session$flushReact()
               session$setInputs(add_chain = 1)
               session$setInputs(add_chain = 2)
               
               chains <- session$returned()$chains
               expect_length(chains, 2)
               ids <- vapply(chains, function(x) x$chain_id, character(1))
               expect_equal(length(unique(ids)), 2)
             })
})

test_that("loading settings fills chains and head group", {
  trigger <- reactiveVal(NULL)
  testServer(mod_LipidDefinition_server, args = list(load_trigger = trigger), {
    session$flushReact()
    trigger(make_settings(n_chains = 2))
    session$flushReact()
    
    out <- session$returned()
    expect_equal(out$head_group, "C6H12O6")
    expect_equal(out$shorthand, "GLC")
    expect_length(out$chains, 2)
    expect_equal(out$chains[[1]]$C_min, 16)
  })
})

test_that("loading a second dataset replaces the first", {
  trigger <- reactiveVal(NULL)
  testServer(mod_LipidDefinition_server, args = list(load_trigger = trigger), {
    session$flushReact()
    trigger(make_settings(n_chains = 3))
    session$flushReact()
    expect_length(session$returned()$chains, 3)
    
    trigger(make_settings(n_chains = 1))
    session$flushReact()
    expect_length(session$returned()$chains, 1)
  })
})

test_that("loaded chain priority is kept", {
  trigger <- reactiveVal(NULL)
  settings <- make_settings(1)
  settings$chains[[1]]$priority <- "important"
  
  testServer(mod_LipidDefinition_server, args = list(load_trigger = trigger), {
    session$flushReact()
    trigger(settings)
    session$flushReact()
    expect_equal(session$returned()$chains[[1]]$priority, "important")
  })
})

test_that("mod_LipidDefinition_ui ui works", {
  ui <- mod_LipidDefinition_ui(id = "test")
  golem::expect_shinytaglist(ui)
  # Check that formals have not been removed
  fmls <- formals(mod_LipidDefinition_ui)
  for (i in c("id")){
    expect_true(i %in% names(fmls))
  }
})
 