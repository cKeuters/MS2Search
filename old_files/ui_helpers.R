#-------------------------------------------------------------------------------
# functions for analysis screen
#-------------------------------------------------------------------------------
createPriorityRadioButtons <- function(id, selected, horizontal = TRUE){
  radioButtons(
    inputId = id,
    label = NULL,
    choiceNames = list(
      tags$span(
        tags$i(style = "color:#28a745;cursor:help;", `data-bs-toggle` = "tooltip", class = "fas fa-square", title = "The filter affects the score.")
      ),
      tags$span(
        tags$i(style = "color:#fd7e14;cursor:help;", `data-bs-toggle` = "tooltip", class = "fas fa-square", title = "The filter affects the score greatly")
      ),
      tags$span(
        tags$i(style = "color:#dc3545;cursor:help;", `data-bs-toggle` = "tooltip", class = "fas fa-square", title = "Spectra not included in the filter are removed")
      )
    ),
    choiceValues = c("basic", "important", "required"),
    selected = selected,
    inline = horizontal
  )
}

createInfoI <- function(id, info) {
  info <- as.character(info)
  
  info_js <- jsonlite::toJSON(info, auto_unbox = TRUE)
  
  tagList(
    tags$span(
      shiny::icon("info-circle"),
      id = id,
      style = "color:#007bff; cursor:pointer;"
    ),
    tags$script(HTML(paste0(
      "(function() {
         var initTooltip = function() {
           $('#", id, "').tooltip({
             html: true,
             sanitize: false,
             title: ", info_js, ",
             trigger: 'hover focus click'
           });
         };

         if (typeof Shiny !== 'undefined' && Shiny.shinyapp && Shiny.shinyapp.isConnected()) {
           initTooltip();
         } else {
           $(document).on('shiny:connected', initTooltip);
         }
       })();"
    )))
  )
}

#------------------------- Peak filter cards -----------------------------------

empty_entry <- function(id) list(id = id, value = "")

empty_card <- function(id) {
  list(
    id = id,
    label = paste("Peak filter", id),
    input_type = "mass",
    entries = list(empty_entry(1)),
    next_entry_sc_mod_id = 2,
    min_repeat = 0,
    max_repeat = 1,
    just_after_precursor = FALSE,
    sumOfSmall = FALSE,
    sequential = FALSE,
    priority = 1,
    next_card_id = NA
  )
}

resolve_mass <- function(card) {
  if (card$input_type == "mass") {
    m <- suppressWarnings(as.numeric(card$value))
    if (!is.na(m)) return(round(m, 4))
    return(NA)
  }
  idx <- which(METABOLITES$name == card$value)
  if (length(idx) == 1) return(METABOLITES$mass[idx])
  NA
}


#-------- Lipid Definition -----------------

#Side chain card
sidechainCardUI <- function(id, index) {
  ns <- NS(id)
  
  div(
    class = "sidechain-card",
    div(class = "sidechain-title", paste("Sidechain", index)),
    numericInput(
      inputId = ns(paste0("carbon_", index)),
      label = "Carbon atoms",
      value = 16,
      min = 0,
      step = 1
    ),
    numericInput(
      inputId = ns(paste0("doublebonds_", index)),
      label = "Double bonds",
      value = 0,
      min = 0,
      step = 1
    ),
    actionButton(ns("add_entry_sidechain_mod"), "Add Modification", icon = icon("plus"), class = "btn-sm"),
    
  )
}

#ui for lipid definition
lipidDefinitionUI <- function(id) {
  ns <- NS(id)
  
  div(
    class = "lipid-definition-card",
    div(class = "card-title", "Lipid Definition"),
    
    div(
      class = "lipid-form-row",
      div(
        class = "headgroup-input",
        textInput(
          inputId = ns("headgroup_formula"),
          label = "Headgroup formula",
          value = "",
          placeholder = "e.g. C3H9O6P"
        )
      ),
      div(
        class = "number-input",
        numericInput(
          inputId = ns("num_sidechains"),
          label = "Number of sidechains",
          value = 2,
          min = 0,
          max = 10,
          step = 1
        )
      )
    ),
    
    div(
      class = "lipid-modifications",
      textInput(
        inputId = ns("modifications"),
        label = "Add possible modifications to formula",
        value = "",
        placeholder = "e.g. +H2O, -NH3"
      )
    ),
    
    uiOutput(ns("sidechains_ui")),
    
    tags$details(
      class = "lipid-advanced",
      tags$summary("Advanced settings"),
      div(
        class = "advanced-content",
        div(
          class = "advanced-regression",
          numericInput(ns("reg_param1"), "Regression parameter 1", value = 0),
          numericInput(ns("reg_param2"), "Regression parameter 2", value = 0),
          numericInput(ns("reg_param3"), "Regression parameter 3", value = 0),
          numericInput(ns("reg_param4"), "Regression parameter 4", value = 0)
        ),
        div(
          class = "lipid-ionnotations",
          selectInput(
            inputId = ns("ionnotations"),
            label = "Possible ionnotations",
            choices = c("[M+H]+", "[M+Na]+", "[M-H]-", "[M+NH4]+", "[M+K]+"),
            selected = "[M+H]+",
            multiple = TRUE
          )
        )
      )
    )
  )
}

lipidDefinitionServer <- function(id, initial) {
  moduleServer(id, function(input, output, session) {
    
    entry_modules_sc_mod <- reactiveValues()
    entry_sc_mod_ids <- reactiveVal(character(0))
    next_entry_sc_mod_id <- reactiveVal(initial)
    
    add_entry_sidechain_mod <- function(entry_init) {
      message("add_entry_sidechain_mod called with: ", entry_init$id, " / ", entry_init$value)
      eid <- as.character(entry_init$id)
      insertUI(
        selector = paste0("#", ns("entries_container")), where = "beforeEnd",
        ui = div(id = ns(paste0("entry_wrap_", eid)),
                 entryModuleUI(ns(paste0("entry_", eid)), entry_init, rule_state$input_type, metabolite_choices))
      )
      entry_modules_sc_mod[[eid]] <- entryModuleServer(paste0("entry_", eid), entry_init, remove_entry)
      entry_sc_mod_ids(c(entry_sc_mod_ids(), eid))
    }
    
    remove_entry_sidechain_mod <- function(eid) {
      eid <- as.character(eid)
      if (length(entry_sc_mod_ids()) <= 1) return(invisible())  # mind. 1 Zeile behalten
      removeUI(selector = paste0("#", ns(paste0("entry_wrap_", eid))))
      entry_modules_sc_mod[[eid]] <- NULL
      entry_sc_mod_ids(setdiff(entry_sc_mod_ids(), eid))
    }
    
    current_entries_sc_mod <- reactive({
      lapply(entry_sc_mod_ids(), function(eid) entry_modules_sc_mod[[eid]]())
    })
    
    output$sidechains_ui <- renderUI({
      req(input$num_sidechains)
      n <- as.integer(input$num_sidechains)
      
      if (is.na(n) || n <= 0) return(NULL)
      
      div(
        class = "sidechains-container",
        lapply(seq_len(n), function(i) {
          sidechainCardUI(id, i)
        })
      )
    })
    
    # Optional: Reaktive Liste aller Sidechain-Werte zurückgeben
    sidechain_values <- reactive({
      req(input$num_sidechains)
      n <- as.integer(input$num_sidechains)
      
      lapply(seq_len(n), function(i) {
        list(
          carbons = input[[paste0("carbon_", i)]],
          double_bonds = input[[paste0("doublebonds_", i)]]
        )
      })
    })
    
    return(list(
      headgroup = reactive(input$headgroup_formula),
      modifications = reactive(input$modifications),
      num_sidechains = reactive(input$num_sidechains),
      sidechains = sidechain_values,
      regression_params = reactive(c(input$reg_param1, input$reg_param2, input$reg_param3, input$reg_param4)),
      ionnotations = reactive(input$ionnotations)
    ))
  })
}

#-------------------------------------------------------------------------------
# Save and Load parameters
#-------------------------------------------------------------------------------


create_settings_list <- function(state, basicFilterValues, toolSettings, cards){
  
  settings_list <- list(
    project_info = list(
      project_name = state$project_name,
      date = as.character(Sys.Date()),
      mode = state$mode,
      out_dir = state$out_dir,
      screen = state$screen
    ),
    LipidIN_settings = list(
      active = toolSettings$LipidIN,
      priority = toolSettings$LipidIN_priority,
      ppm1 = toolSettings$LipidIN_ppm1,
      ppm2 = toolSettings$LipidIN_ppm2
    ),
    DreaMS_settings = list(
      active = toolSettings$DreaMS,
      priority = toolSettings$DreaMS_priority
    ),
    basicFilters = list(
      mz_min = basicFilterValues$mz_min ,
      mz_max = basicFilterValues$mz_max ,
      mz_priority = basicFilterValues$mz_priority,
      rt_min = basicFilterValues$rt_min ,
      rt_max = basicFilterValues$rt_max ,
      rt_priority = basicFilterValues$rt_priority,
      n_peaks_min = basicFilterValues$n_peaks_min ,
      n_peaks_max = basicFilterValues$n_peaks_max,
      n_peaks_priority = basicFilterValues$n_peaks_priority,
      tolerance_ppm = basicFilterValues$tolerance_ppm,
      minimum_dalton = basicFilterValues$minimum_dalton
    )
  )
  settings_list <- c(
    settings_list,
    list(specific_filters = cards)
  )
  
  return(settings_list)
}


setup_and_save_project <- function(settings_list, output_zip_path_file) {

  tmp_dir <- tempfile("CeraSearch_save_dir")
  dir.create(tmp_dir)

  #arrow::write_parquet(ms2_data, file.path(tmp_dir, "ms2_data.parquet"))
  
  jsonlite::write_json(settings_list, file.path(tmp_dir, "settings.json"), auto_unbox = TRUE)
  
  #zip files
  files_to_zip <- list.files(tmp_dir, full.names = TRUE)
  zip::zip(zipfile = output_zip_path_file, files = files_to_zip, mode = "cherry-pick")
  
  print(paste("File saved at", output_zip_path_file))
  
  unlink(tmp_dir, recursive = TRUE)
}

load_project_settings <- function(zip_path) {
  tmp_dir <- tempfile("CeraSearch_load_dir")
  dir.create(tmp_dir)
  zip::unzip(zip_path, exdir = tmp_dir)
  
  settings_list <- jsonlite::fromJSON(
    file.path(tmp_dir, "settings.json"),
    simplifyVector = FALSE  
  )
  
  unlink(tmp_dir, recursive = TRUE)
  settings_list
}

load_user_project <- function(uploaded_zip_path) {
  tmp_dir <- tempfile("CeraSearch_save_dir")
  dir.create(tmp_dir)
  
  zip::unzip(uploaded_zip_path, exdir = tmp_dir)
  
  # Validation
  # required_files <- c("ms2_data.parquet", "settings.json")
  # if (!all(required_files %in% list.files(tmp_dir))) {
  #   stop("Files missing")
  # }
  
  settings <- jsonlite::read_json(file.path(tmp_dir, "settings.json"), simplifyVector = TRUE)
  
  unlink(tmp_dir, recursive = TRUE)
  
  return(list(ms1 = ms1, ms2 = ms2, settings = settings))
}
#ms2 <- arrow::read_parquet(file.path(tmp_dir, "ms2_data.parquet"))
#ms1_dataset <- arrow::open_dataset für zu größe Datensätze



#-------------------------------------------------------------------------------
# Data used by the peak filter
#-------------------------------------------------------------------------------
METABOLITES <- data.frame(
  name = c(
    "Glucose", "Fructose", "Sucrose", "Lactate", "Pyruvate",
    "Citrate", "Malate", "ATP", "ADP", "AMP"
  ),
  mass = c(
    180.0634, 180.0634, 342.1162, 89.0244, 87.0082,
    192.0270, 133.0137, 506.9957, 427.0294, 347.0631
  ),
  stringsAsFactors = FALSE
)

