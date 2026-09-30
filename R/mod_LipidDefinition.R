#' LipidDefinition UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList textInput actionButton 
mod_LipidDefinition_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(
      class = "lipid-def-panel",
      h3("Lipid Definition"),
      div(
        class = "headgroup-row",
        div(
          class = "headgroup-box",
          textInput(ns("headgroup"), label = "Backbone", placeholder = "z.B. C3H8O3", width = "160px")
        ),
        actionButton(ns("add_chain"), "+ Chain", class = "btn-add-chain")
      ),
      div(id = ns("chains_container"), class = "chains-container"),
      tags$div(
        class = "precursor-validation",
        tags$details(
          tags$summary(
            class = "validation-summary",
            icon("flask"), "Validation: Precursor m/z"
          ),
          div(id = ns("validation_container")),
          actionButton(ns("add_validation"), "Add validation")
        )
      )
     
    )
  )
  
}



#' LipidDefinition Server Functions
#'
#' @noRd 
mod_LipidDefinition_server <- function(id, load_trigger, head_group){
  moduleServer(id, function(input, output, session){
    ns <- session$ns
    

    chain_modules <- reactiveValues()
    chain_ids <- reactiveVal(character(0))
    next_chain_id <- reactiveVal(1)
    
    
    observeEvent(input$headgroup, {
      head_group(input$headgroup)
      })
    
    observeEvent(head_group(), {
      updateTextInput(
        session = session,
        "headgroup",
        value = head_group()
      )
    })
    
    
    remove_chain <- function(cid) {
      cid <- as.character(cid)
      cat(paste("remove chain", cid))
      removeUI(selector = paste0("#", ns(paste0("chain_wrap_", cid))))
      chain_modules[[cid]] <- NULL

      chain_ids(setdiff(chain_ids(), cid))
    }
    
    observeEvent(input$add_chain, {
      add_chain_to_ui(init_chain = empty_chain)
    }, ignoreInit = TRUE)
    
    add_chain_to_ui <- function(init_chain){
      cid <- as.character(next_chain_id())
      chain_ids(c(chain_ids(), cid))
      
      child_id <- paste0("chain_", cid)
      
      insertUI(
        selector = paste0("#", ns("chains_container")),
        where = "beforeEnd",
        ui = div(
          id = ns(paste0("chain_wrap_", cid)),
          chain_block_ui(ns(child_id), init_chain)
        ),
        immediate = TRUE
      )
      
      chain_modules[[cid]] <- chain_block_server(child_id, cid, init_chain, remove_chain) 
      #validation_modules[[1]]$components <- c(validation_modules[[1]]$components, child_id)
      next_chain_id(next_chain_id() + 1)
      
    }
    
    
    validation_modules <- list()        
    validation_ids  <- reactiveVal(character(0))
    next_validation_id <- reactiveVal(1)
    
    remove_validation <- function(vid) {
      vid <- as.character(vid)
      removeUI(selector = paste0("#", ns(paste0("validation_wrap_", vid))),
               immediate = TRUE)
      validation_modules[[vid]] <<- NULL
      validation_ids(setdiff(validation_ids(), vid))
    }
    
    add_validation_to_ui <- function(init_validation, fixed = FALSE) {
      vid <- as.character(isolate(next_validation_id()))
      validation_ids(c(isolate(validation_ids()), vid))
      child_id <- paste0("validation_", vid)
      
      insertUI(
        selector = paste0("#", ns("validation_container")),
        where = "beforeEnd",
        ui = div(id = ns(paste0("validation_wrap_", vid)),      # <- unique prefix
                 validationUI(ns(child_id), init_validation, fixed)),
        immediate = TRUE
      )
      validation_modules[[vid]] <<- validationServer(
        child_id, init_validation, chain_ids, fixed, on_remove = if(!fixed) remove_validation(vid) else NULL)
      next_validation_id(as.numeric(vid) + 1)
    }

    
    # always-present precursor validation
    session$onFlushed(function() {
      isolate(add_validation_to_ui(list(name = "precursor validation",  adjustments = "", val_priority = "basic",  components = list()), fixed = TRUE))
    }, once = TRUE)
    
    observeEvent(input$add_validation, {
      add_validation_to_ui(empty_validation)
    })
    
    observeEvent(load_trigger(), {
      req(load_trigger())
      chains_to_load <- load_trigger()$chains
      cat("load_trigger fired, n chains:", length(chains_to_load), "\n")
      isolate({
        for (chain in chains_to_load) {
          tryCatch({
            add_chain_to_ui(chain)
          }, error = function(e) {
            cat("ERROR in add_chain_to_ui:", conditionMessage(e), "\n")
          })
        }
      })
    }, ignoreNULL = TRUE, ignoreInit = TRUE)
    
    reactive({list(
      chains = ({
        ids <- chain_ids()
        if (length(ids) == 0) return(list())
        lapply(ids, function(cid) chain_modules[[cid]]())
      }),
      validations = ({
        if (length(validation_ids()) == 0) return(list())
        lapply(validation_ids(), function(v){
          validation_modules[[v]]()
        })
      })
      
    )
    })
    
  })
}


#' UI-Block for one chain
#' @noRd
chain_block_ui <- function(id, init_chain) {
  ns <- NS(id)
  div(
    id = ns("chain_block_"),
    class = "chain-block",
    div(
      class = "chain-controls",
      actionButton(ns("remove_chain"), "\u00d7", class = "btn-remove-chain"),
      div(
        class = "control-group",
        tags$label("C-Atoms"),
        div(
          class = "control-inputs",
          numericInput(ns("c_min"), "min", value = init_chain$C_min, min = 1, max = 100, width = "62px"),
          numericInput(ns("c_max"), "max", value = init_chain$C_max, min = 1, max = 100, width = "62px")
        )
      ),
      div(
        class = "control-group",
        tags$label("Double Bounds"),
        div(
          class = "control-inputs",
          numericInput(ns("db_min"), "min", value = init_chain$DB_min, min = 0, max = 50, width = "62px"),
          numericInput(ns("db_max"), "max", value = init_chain$DB_max, min = 0, max = 50, width = "62px")
        )
      ),
      div(
        class = "control-group",
        tags$label("Modifikation"),
        actionButton(ns("add_mod"), "+", class = "btn-modification")
      )
    ),
    div(id = ns("mods_above"), class = "mods-row mods-above"),
    div(
      class = "chain-row",
      uiOutput(ns("chain_priority")),
      div(
        class = "chain-visual",
        div(class = "chain-connector"),
        uiOutput(ns("chain_svg"))
      ),
      uiOutput(ns("mz_box_ui"))
    ),
    div(id = ns("mods_below"), class = "mods-row mods-below")
  )
}

chain_block_server <- function(server_id, id, init_chain, on_remove ) {
  moduleServer(server_id, function(input, output, session) {
    ns <- session$ns
    cat("chain_block_server id =", server_id, "\n")
    
    chain_state <- reactiveValues(
      C_min = init_chain$C_min,
      C_max = init_chain$C_max,
      DB_min = init_chain$DB_min,
      DB_max = init_chain$DB_max,
      priority = "basic",
      modifications = init_chain$modifications
    )
    mod_ids <- reactiveVal(1)
    
    chain_mz_values <- reactive({
      calc_lipid_chain_masses(reactiveValuesToList(chain_state))
    })
    
    output$chain_priority <- renderUI({
      createPriorityRadioButtons(ns("chain_priority"), chain_state$priority)
    })
    
    observeEvent(input$c_min, {
      if(input$c_min < 0) updateNumericInput(session, "c_min", value = 0)
      if(input$c_min > 100) updateNumericInput(session, "c_min", value = 100)
      chain_state$C_min <- input$c_min
      updateNumericInput(
        session, "c_max",
        min = chain_state$C_min,
        value = max(input$c_max, chain_state$C_min)
      )
    })
    
    observeEvent(input$c_max, {
      if(input$c_max < 0) updateNumericInput(session, "c_max", value = 0)
      if(input$c_max > 100) updateNumericInput(session, "c_max", value = 100)
      chain_state$C_max <- input$c_max
      test_maxDB()
    })
    
    observeEvent(input$db_min, {
      if(input$db_min < 0) updateNumericInput(session, "db_min", value = 0)
      if(input$db_min > 50) updateNumericInput(session, "db_min", value = 50)
      chain_state$DB_min <- input$db_min
      updateNumericInput(
        session, "db_max",
        min = chain_state$DB_min,
        value = max(input$db_max, chain_state$DB_min)
      )
    })
    
    observeEvent(input$db_max, {
      if(input$db_max < 0) updateNumericInput(session, "db_max", value = 0)
      if(input$db_max > 50) updateNumericInput(session, "db_max", value = 50)
      chain_state$DB_max <- input$db_max
    })
    
    test_maxDB <- function(){
      max_db <- floor(chain_state$C_max / 2)
      updateNumericInput(
        session, "db_max",
        max = max_db,
        value = min(input$db_max, max_db)
      )
      updateNumericInput(
        session, "db_min",
        max = max_db,
        value = min(input$db_min, max_db)
      )
    }
    
    observeEvent(input$remove_chain, on_remove(id), ignoreInit = TRUE)
    
    observeEvent(input$add_mod,{
      init_mod <- list(formula = "", optional = TRUE)
      add_mod_to_ui(init_mod)
      
    })
    
    register_mod_observers <- function(mid, init_mod) {
      chain_state$modifications[[mid]] <- list(
        formula = init_mod$formula,
        optional = init_mod$optional
      )
      observeEvent(input[[paste0(mid, "_formula")]], {
        chain_state$modifications[[mid]]$formula <- input[[paste0(mid, "_formula")]]
      })
      observeEvent(input[[paste0(mid, "_mod_optional")]], {
        chain_state$modifications[[mid]]$optional <- input[[paste0(mid, "_mod_optional")]]
      })
      observeEvent(input[[paste0(mid, "_remove")]], {
        chain_state$modifications[[mid]] <- NULL
        removeUI(selector = paste0("#", ns(paste0("mod_chip_", mid))))
      })
    }
    
    add_mod_to_ui <- function(init_mod) {
      n_mods <- mod_ids()
      mid <- paste0("mods_", n_mods)
      position <- if (n_mods %% 2 == 0) "above" else "below"
      target <- paste0("mods_", position)
      
      insertUI(
        selector = paste0("#", ns(target)),
        where = "beforeEnd",
        ui = mod_chip_ui(ns, mid, init_mod)
      )
      register_mod_observers(mid, init_mod)
      mod_ids(mod_ids() + 1)
    }
    
    isolate({
      init_mods <- init_chain$modifications
      for (init_mod in init_mods) {
        add_mod_to_ui(init_mod)
        }
    })
    
    
    output$chain_svg <- renderUI({
      req(chain_state$C_min, chain_state$C_max)
      HTML(
        build_chain_svg(
          chain_state$C_min,
          chain_state$C_max,
          chain_state$DB_min,
          chain_state$DB_max
        )
      )
    })
    
    infoModalButtonServer("chain_mz_info", df = chain_mz_values)
    #' box with min and max mz values of chain
    #' @noRd
    output$mz_box_ui <- renderUI({
      values <- chain_mz_values()
      req(values)
      div(class = "mz_box_distance",
        div(
          class = "mz-box",
          tags$span(class = "mz-label", "m/z"),
          tags$span(class = "mz-value", sprintf("%.4f", min(values$mass))),
          tags$span(class = "mz-sep", "\u2013"),
          tags$span(class = "mz-value", sprintf("%.4f", max(values$mass)))
        ),
        infoModalButtonUI(ns("chain_mz_info"))
      )
        }
    ) 
    
    reactive({
      list(
        C_min = chain_state$C_min,
        C_max = chain_state$C_max,
        DB_min = chain_state$DB_min,
        DB_max = chain_state$DB_max,
        priority = chain_state$priority,
        modifications = chain_state$modifications
      )
    })
    
  })
}

#' UI-Chip for one modification
#' @noRd
mod_chip_ui <- function(ns, mid, init_mod) {
  div(
    id = ns(paste0("mod_chip_", mid)),
    class = "mod-chip",
    textInput(ns(paste0(mid, "_formula")), NULL, placeholder = "e.g. OH", value = init_mod$formula,  width = "72px"),
    checkboxInput(ns(paste0(mid,"_mod_optional")), "Optional", value = init_mod$optional),
    actionButton(ns(paste0(mid, "_remove")), "\u00d7", class = "btn-remove-mod")
  )
}

mod_chip_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    reactive({
      list(
        formula = input$formula,
        optional = input$optional
      )
    })
  })
}

# UI of table pop-up
infoModalButtonUI <- function(id) {
  ns <- NS(id)
  actionLink(
    ns("show_info"),
    label = NULL,
    icon = icon("info-circle"),
    class = "infoModalButtonUI-style"
  )
}

# Server of table pop-up
infoModalButtonServer <- function(id, df, title = NULL) {
  moduleServer(id, function(input, output, session) {
    observeEvent(input$show_info, {
      showModal(modalDialog(
        title = title,
        DT::dataTableOutput(session$ns("info_table")),
        easyClose = TRUE,
        size = "l",
        footer = modalButton("Close")
      ))
    })
    
    output$info_table <- DT::renderDataTable({
      DT::datatable(df(),selection = 'none', options = list(dom = "t", paging = FALSE, scrollY = "400px"))
    })
  })
}

validationUI <- function(id, validation, fixed = FALSE) {
  ns <- NS(id)
  div(
    class = "validation-block",
    div(class = "validation-header",
        if(fixed){
          tags$span(class = "formula-label", paste(validation$name, "="))
        }else{
          tagList(textInput(ns("val_name"), label = NULL,
                            value = validation$name,
                            placeholder = "Validation name",
                            width = "auto") |>
            tagAppendAttributes(class = "validation-name"),
            tags$span(class = "formula-label", "="))
        }),
    
    if (!fixed) tagList(
      div(id = ns("selector_container")),
      actionButton(ns("add_component"), "+")
    ) else tags$span(class = "formula-term", "Backbone + all chains"),
    div(
      class = "adjustment-row",
      icon("pen", class = "adjustment-icon"),
      textInput(ns("adjustment"), label = NULL,
                value = validation$adjustments %||% "",
                placeholder = "individual adjustment")
    ),
    if(!fixed){
      actionButton(ns("remove_validation"), label = NULL, icon = icon("trash"), class = "btn-remove-validation")
    },
    createPriorityRadioButtons(ns("val_priority"),
                               validation$val_priority,
                               horizontal = TRUE)
  )
}

validationServer <- function(id, init_validation, chain_ids, fixed = FALSE, on_remove = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    sel_ids  <- reactiveVal(character(0))
    next_sel <- reactiveVal(1)
    
    add_component_to_ui <- function(selection = NULL) {
      i <- isolate(next_sel())
      next_sel(i + 1)
      sel_ids(c(isolate(sel_ids()), as.character(i)))
      insertUI(
        selector = paste0("#", ns("selector_container")),
        where = "beforeEnd",
        ui = selectInput(ns(paste0("sel_", i)), label = NULL,
                         choices = c("Backbone", isolate(chain_ids())),
                         selected = selection),
        immediate = TRUE
      )
    }
    
    if (!fixed) {
      observeEvent(input$add_component, add_component_to_ui())
      for (com in init_validation$components) add_component_to_ui(com)
      
      # keep the choices in sync when chains are added/removed
      observe({
        choices <- c("Backbone", chain_ids())
        for (i in sel_ids()) {
          cur <- isolate(input[[paste0("sel_", i)]])
          updateSelectInput(session, paste0("sel_", i), choices = choices,
                            selected = if (!is.null(cur) && cur %in% choices) cur else NULL)
        }
      })
    }
    
    components <- reactive({
      if (fixed) {
        c("Backbone", chain_ids())
      } else {
        unique(unlist(lapply(sel_ids(), function(i) input[[paste0("sel_", i)]])))
      }
    })
    
    observeEvent(input$remove_validation,{
      if(is.function(on_remove)) on_remove()
    }, ignoreInit = TRUE, once = TRUE)
    reactive(list(
      name         = input$val_name %||% init_validation$name,
      adjustments  = input$adjustment %||% init_validation$adjustments,
      val_priority = input$val_priority %||% init_validation$val_priority,
      components   = components()
    ))
  })
}





#' Baut das SVG-Markup fuer eine Acylkette
#'
#' @param c_min,c_max max bestimmt die gezeichnete Kettenlaenge, ab min wird
#'   die Kette grau (Bereich der variablen Laenge).
#' @param db_min,db_max db_max Doppelbindungen werden gleichmaessig ueber die
#'   Kette verteilt, davon sind die ersten db_min schwarz, der Rest grau.
#' @noRd
build_chain_svg <- function(c_min, c_max, db_min = 0, db_max = 0) {
  c_max  <- max(c_max, c_min, 1)
  if(c_max > 100)return()
  if(db_max > 50)return() 
  db_max <- max(db_max, db_min, 0)
  
  dx <- 24; dy <- 16
  x0 <- 16; y0 <- 46
  n  <- c_max
  
  xs <- x0 + (0:(n - 1)) * dx
  ys <- y0 + rep(c(0, -dy), length.out = n)
  width <- x0 + n * dx + 20
  
  n_segs <- n - 1
  db_positions <- if (db_max > 0 && n_segs > 0) {
    if (db_max == 1) 1L else round(seq(1, n_segs, length.out = db_max))
  } else integer(0)
  
  segs <- character(0)
  for (i in seq_len(n_segs)) {
    black <- i <= (c_min - 1)
    color <- if (black) "#1a1a1a" else "#b0b0b0"
    x1 <- xs[i]; y1 <- ys[i]; x2 <- xs[i + 1]; y2 <- ys[i + 1]
    
    segs <- c(segs, sprintf(
      '<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="2.5" stroke-linecap="round"/>',
      x1, y1, x2, y2, color
    ))
    
    db_idx <- which(db_positions == i)
    if (length(db_idx) > 0) {
      dxs <- x2 - x1; dys <- y2 - y1
      len <- sqrt(dxs^2 + dys^2)
      ox <- -dys / len * 4
      oy <- dxs / len * 4
      db_color <- if (db_idx <= db_min) color else "#c9c9c9"
      segs <- c(segs, sprintf(
        '<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="2" stroke-linecap="round"/>',
        x1 + ox, y1 + oy, x2 + ox, y2 + oy, db_color
      ))
    }
  }
  
  dots <- character(0)
  for (i in seq_len(n)) {
    color <- if (i <= c_min) "#1a1a1a" else "#b0b0b0"
    dots <- c(dots, sprintf('<circle cx="%.1f" cy="%.1f" r="2.2" fill="%s"/>', xs[i], ys[i], color))
  }
  sprintf(
    '<svg viewBox="0 0 %d 70" width="100%%" height="70" class="lipid-chain-svg">%s%s</svg>',
    width, paste(segs, collapse = ""), paste(dots, collapse = "")
  )
}

#' @import Rdisop
#' @noRd
calc_chain_mz_range <- function(chain_state){
  
  mz_mod_min <- 0
  mz_mod_max <- 0
  for(mod in chain_state$modifications){
    if(mod$formula == "")next
    molecule <- Rdisop::getMolecule(mod$formula)
    mz_mod_max <- mz_mod_max + Rdisop::getMass(molecule) 
    if(!mod$optional)mz_mod_min <- mz_mod_min + Rdisop::getMass(molecule) 
  }
  chain_mass_min <- ((chain_state$C_min * atom_masses$C) 
                     + (chain_state$C_min * 2 * atom_masses$H) 
                     - (sum(vapply(chain_state$modifications,function(x) isFALSE(x$optional) & x$formula != "",logical(1)))* atom_masses$H)
                     - (chain_state$DB_max * 2 * atom_masses$H)
                     + mz_mod_min)
  
  chain_mass_max <- ((chain_state$C_max * atom_masses$C) 
                     + (chain_state$C_max * 2 * atom_masses$H) 
                     - (sum(vapply(chain_state$modifications,function(x) x$formula != "",logical(1)))* atom_masses$H)
                     - (chain_state$DB_min * 2 * atom_masses$H)
                     + mz_mod_max)
  
  
  list(min= chain_mass_min, max = chain_mass_max)
}

atom_masses <- list(
  H = 1.00783,
  C = 12.00000,
  N = 14.00307,
  O = 15.99491,
  P = 30.97376,
  S = 31.97207
)

empty_chain <- list(
  C_min = 10,
  C_max = 20,
  DB_min = 0,
  DB_max = 1,
  priority = "basic",
  modifications = list()
)

empty_validation <- list(
  name = "enter name",
  adjustments = NULL,
  val_priority = "basic",
  selectedComponents = list()
)
