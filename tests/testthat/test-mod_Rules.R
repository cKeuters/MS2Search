test_that("mod_Rules_server works", {
  id <- "test"
  
  initial <- list(
    id = 1,
    label = "Test rule",
    is_mzPeak = FALSE,
    input_type = "mass", 
    min_repeat = 1,
    max_repeat = 1,
    just_after_precursor = FALSE,
    filter_noise = FALSE,
    sequential = FALSE,
    next_card_id = 2,
    rule_priority = "basic",
    rule_mz_min = NA,
    rule_mz_max = NA,
    entries = list(),   
    next_entry_id = 1
  )
  
  testServer(
    mod_Rules_server,
    args = list(
      initial = initial,
      polarity = "positive",
      on_delete = function(...) NULL,
      on_duplicate = function(...) NULL,
      metabolite_choices = c("Cer", "HexCer")
      ),
    {
      expect_true(inherits(session$ns, "function"))
      expect_true(grepl("test", session$ns("test")))
    }
  )
})
