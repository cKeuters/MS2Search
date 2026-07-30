library(shiny)
library(shinyFiles)
library(shinyWidgets)
library(bslib)

#loading and saving data
library(arrow)
library(jsonlite)
library(zip)

#run_analysis
library(Spectra)
library(MsBackendMgf)
library(SpectraQL)

options(shiny.maxRequestSize = 500 * 1024^2)
#shiny::runApp( appDir = "CeraSearch/CeraSearch_Test/", launch.browser = TRUE)

source("ui_helpers.R")
source("run_analysis.R")

#-------------------------------------------------------------------------------
# Main Server
#-------------------------------------------------------------------------------
function(input, output, session) {
  state <- reactiveValues(
    screen = "welcome",
    out_dir = path.expand("~"),
    mode = "pos",
    project_name = ""
  )
  
  basicFilterValues <- reactiveValues(
    mz_min = 0,
    mz_max = 1000,
    mz_priority = "required",
    rt_min = 0,
    rt_max = 1000,
    rt_priority = "required",
    n_peaks_min = 0,
    n_peaks_max = 1000,
    n_peaks_priority = "required",
    tolerance_ppm = 10,
    minimum_dalton = 0.01
  )
  
  ToolSettings <- reactiveValues(
    LipidIN = TRUE,
    LipidIN_ppm1 = 10,
    LipidIN_ppm2 = 50,
    LipidIN_priority = "basic",
    DreaMS = TRUE,
    DreaMS_priority = "basic"
  )
  
  lipid_data <- lipidDefinitionServer("lipid", 1)
  
  #------------------------------- Welcome Screen ---------------------------------------
  
  #roots for folder selection
  roots <- c(
    home = normalizePath(path.expand("~"), winslash = "/"),
    mount = "/mnt/ifs/data/Workspace"
  )
  shinyDirChoose(input, "out_dir", roots = roots, session = session)
  shinyDirChoose(input, "load_project", roots = roots, session = session)
  
  # show output path
  observeEvent(input$out_dir, {
    req(input$out_dir)
    path <- as.character(parseDirPath(roots, input$out_dir))
    if (length(path) > 0 && nzchar(path)) {
      state$out_dir <- path
    }
  })
  output$out_dir_text <- renderText({
    state$out_dir
  })
  
  # Load Project
  observeEvent(input$load_project_file, {
    req(input$load_project_file)
    settings_list <- load_project_settings(input$load_project_file$datapath)
    
    pi <- settings_list$project_info
    state$project_name <- pi$project_name
    state$mode <- pi$mode
    state$out_dir <- pi$out_dir
    state$screen <- pi$screen
    
    li <- settings_list$LipidIN_settings
    ToolSettings$LipidIN <- li$active
    ToolSettings$LipidIN_priority <- li$priority
    ToolSettings$LipidIN_ppm1 <- li$ppm1
    ToolSettings$LipidIN_ppm2 <- li$ppm2
    
    dr <- settings_list$DreaMS_settings
    ToolSettings$DreaMS <- dr$active
    ToolSettings$DreaMS_priority <- dr$priority
    
    bf <- settings_list$basicFilters
    basicFilterValues$mz_min <- bf$mz_min
    basicFilterValues$mz_max <- bf$mz_max
    basicFilterValues$mz_priority <- bf$mz_priority
    basicFilterValues$rt_min <- bf$rt_min
    basicFilterValues$rt_max <- bf$rt_max
    basicFilterValues$rt_priority <- bf$rt_priority
    basicFilterValues$n_peaks_min <- bf$n_peaks_min
    basicFilterValues$n_peaks_max <- bf$n_peaks_max
    basicFilterValues$n_peaks_priority <- bf$n_peaks_priority
    basicFilterValues$tolerance_ppm <- bf$tolerance_ppm
    basicFilterValues$minimum_dalton <- bf$minimum_dalton
    
    #remove existing cards
    isolate({
      for (cid in card_ids()) {
        removeUI(selector = paste0("#wrap_card_", cid))
        card_modules[[cid]] <- NULL
      }
      card_ids(character(0))
      
      #cards from json
      loaded_cards <- settings_list$specific_filters
      max_card_id <- 0
      
      if(!is.null(loaded_cards)){
        for (c in loaded_cards) {
          entries <- lapply(c$entries, function(e) {
            list(id = as.numeric(e$id), value = e$value)
          })
          max_entry_id <- if (length(entries) == 0) 0 else max(sapply(entries, `[[`, "id"))
          
          card_init <- list(
            id = as.numeric(c$id),
            label = c$label,
            input_type = c$input_type,
            entries = entries,
            next_entry_id = max_entry_id + 1,
            min_repeat = c$min_repeat,
            max_repeat = c$max_repeat,
            just_after_precursor = isTRUE(c$just_after_precursor),
            sumOfSmall = isTRUE(c$sumOfSmall),
            sequential = isTRUE(c$sequential),
            next_card_id = if (is.null(c$next_card_id)) NA else as.numeric(c$next_card_id)
          )
          
          add_card_to_ui(card_init)
          max_card_id <- max(max_card_id, card_init$id)
        }
      }
      next_card_id(max_card_id + 1)
    })
    state$screen <- "analysis"
  })
  
  # States
  observeEvent(input$mode,         { state$mode <- input$mode })
  observeEvent(input$project_name, { state$project_name <- input$project_name })
  
  # Start Rule and Filter selection
  observeEvent(input$start_analysis, {
    if(input$project_name == "") state$project_name <- paste0("project_",format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
    state$screen <- "analysis"
  })
  
  observeEvent(input$back_to_welcome, {
    state$screen <- "welcome"
  })

  #rule cards
  all_cards_data <- reactive({
    ids <- card_ids()
    if (length(ids) == 0) return(list())
    lapply(ids, function(cid) card_modules[[cid]]())
  })
  
  
  #Starting Analysis
  filtered_result <- eventReactive(input$run_analysis, {
    filterMS2data(input$ms2_file$datapath, all_cards_data())
  })
  
  output$result_text <- renderText({
    filtered_result()
  })
  
  observeEvent(input$run_analysis, {
    # setup_and_save_project(create_settings_list(
    #   state = reactiveValuesToList(state),
    #   basicFilterValues = reactiveValuesToList(basicFilterValues),
    #   toolSettings = reactiveValuesToList(ToolSettings),
    #   cards   = all_cards_data()), 
    #   paste0(state$out_dir, "/", state$project_name,".zip")
    # )
    state$screen <- "results"
  })

  observeEvent(input$save_project, {
    setup_and_save_project(create_settings_list(
      state = reactiveValuesToList(state),
      basicFilterValues = reactiveValuesToList(basicFilterValues),
      toolSettings = reactiveValuesToList(ToolSettings),
      cards   = all_cards_data()),
      paste0(state$out_dir, "/", state$project_name,".zip")
    )
  })
  
  MS2file_info <- reactive({
    if (is.null(input$ms2_file)) {
      NULL
    } else {
      input$ms2_file$name
    }
  })
  
  output$MS2file_info <- renderText({
    MS2file_info()
  })
  
  
  #--------------------------- observeEvent for variables -----------------------------
  
  #basic filter events
  observeEvent(input$mz_min, {
    basicFilterValues$mz_min <- input$mz_min
  })
  observeEvent(input$mz_max, {
    basicFilterValues$mz_max <- input$mz_max
  })
  observeEvent(input$mz_priority, {
    basicFilterValues$mz_priority <- input$mz_priority
  })
  observeEvent(input$rt_min, {
    basicFilterValues$rt_min <- input$rt_min
  })
  observeEvent(input$rt_max, {
    basicFilterValues$rt_max <- input$rt_max
  })
  observeEvent(input$rt_priority, {
    basicFilterValues$rt_priority <- input$rt_priority
  })
  observeEvent(input$n_peaks_min, {
    basicFilterValues$n_peaks_min <- input$n_peaks_min
  })
  observeEvent(input$n_peaks_max, {
    basicFilterValues$n_peaks_max <- input$n_peaks_max
  })
  observeEvent(input$n_peaks_max, {
    basicFilterValues$n_peaks_max <- input$n_peaks_max
  })
  observeEvent(input$n_peaks_priority, {
    basicFilterValues$n_peaks_priority <- input$n_peaks_priority
  })
  
  #Tools setting events
  observeEvent(input$LipidIN_active, {
    ToolSettings$LipidIN <- input$LipidIN_active
  })
  observeEvent(input$LipidIN_ppm1, {
    ToolSettings$LipidIN_ppm1 <- input$LipidIN_ppm1
  })
  observeEvent(input$LipidIN_ppm2, {
    ToolSettings$LipidIN_ppm2 <- input$LipidIN_ppm2
  })
  observeEvent(input$LipidIN_priority, {
    ToolSettings$LipidIN_priority <- input$LipidIN_priority
  })
  observeEvent(input$DreaMS_active, {
    ToolSettings$DreaMS <- input$DreaMS_active
  })
  observeEvent(input$DreaMS_priority, {
    ToolSettings$DreaMS_priority <- input$DreaMS_priority
  })
  
  #Tolerance
  
  observeEvent(input$tolerance_ppm, {
    basicFilterValues$tolerance_ppm <- input$tolerance_ppm
  })
  
  observeEvent(input$minimum_dalton, {
    basicFilterValues$minimum_dalton <- input$minimum_dalton
  })
  
  #------------------------------- Peak Filter Logic ---------------------------------------
  
  # Update card fields
  entryModuleUI <- function(id, entry, input_type, metabolite_choices) {
    ns <- NS(id)
    div(class = "entry-row", style = "display:flex; align-items:center; gap:8px;",
        if (input_type == "mass") {
          textInput(ns("value"), NULL, value = entry$value, placeholder = "e.g. 162.0528")
        } else {
          selectInput(ns("value"), NULL, choices = metabolite_choices, selected = entry$value)
        },
        actionButton(ns("remove"), NULL, icon = icon("times"), class = "btn-sm btn-outline-danger", style = "margin-bottom: 15px;")
    )
  }
  
  entryModuleServer <- function(id, initial, on_remove) {
    moduleServer(id, function(input, output, session) {
      state <- reactiveValues(id = initial$id, value = initial$value)
      observeEvent(input$value, state$value <- input$value, ignoreInit = TRUE)
      observeEvent(input$remove, on_remove(initial$id), ignoreInit = TRUE)
      reactive(list(id = state$id, value = state$value))
    })
  }
  
  
  # Render peak filter cards
  cardModuleUI <- function(id, card, metabolite_choices) {
    ns <- NS(id)
    card_class <- "peak-card"
    #if (card$just_after_precursor) card_class <- paste(card_class, "precursor-linked")
    
    div(class = "peak-card-wrap", id = paste0("wrap_", id),
        div(class = card_class,
            div(class = "card-header",
                div(class = "card-title",
                    textInput(ns("label"), NULL, value = card$label)
                ),
                div(class = "card-actions",
                    actionButton(ns("duplicate"), "Duplicate", icon = icon("copy")),
                    actionButton(ns("delete"), "Delete", icon = icon("trash"), class = "btn-danger")
                )
            ),
            div(class = "card-body",
                checkboxInput(ns("sequential"), "Sequential neutral losses", value = card$sequential),
                noUiSliderInput(
                  inputId = ns("slider"),
                  label = "Number of losses",
                  min = 0,
                  max = 10,
                  value = c(card$min_repeat, card$max_repeat),
                  step = 1,
                  format = wNumbFormat(decimals = 0),
                  tooltips = TRUE,
                  pips = list(mode = "values", values = list(0, 10), density = 10),
                  color = "#337ab7"
                ),
                selectInput(ns("input_type"), "Input type",
                            choices = c("Exact mass" = "mass", "Metabolite name" = "metabolite"),
                            selected = card$input_type),
                div(id = ns("entries_container"), class = "entries-container"),
                actionButton(ns("add_entry"), "Add value", icon = icon("plus"), class = "btn-sm"),
                uiOutput(ns("mass_preview")),
            
                checkboxInput(ns("precursor"), "Just after precursor", value = card$just_after_precursor),
                checkboxInput(ns("sumOfSmall"), "Validation by sum of smaller peaks", value = card$sumOfSmall),
                div(style = "text-align: center;", createPriorityRadioButtons(ns("priority"), "basic"))
                
            )
        )
    )
  }
  
  cardModuleServer <- function(id, initial, on_delete, on_duplicate, metabolite_choices, parent_session ) {
    moduleServer(id, function(input, output, session) {
      ns <- session$ns
      
      rule_state <- reactiveValues(
        id = initial$id, 
        label = initial$label, 
        input_type = initial$input_type,
        min_repeat = initial$min_repeat, 
        max_repeat = initial$max_repeat,
        just_after_precursor = initial$just_after_precursor, 
        sumOfSmall = initial$sumOfSmall, 
        sequential = initial$sequential,
        next_card_id = initial$next_card_id
      )
      
      entry_modules <- reactiveValues()
      entry_ids <- reactiveVal(character(0))
      next_entry_id <- reactiveVal(initial$next_entry_id)
      
        observeEvent(input$sequential, {
        shinyjs::toggle(id = "slider", condition = input$sequential)
      }, ignoreNULL = FALSE)
          
      add_entry <- function(entry_init) {
        message("add_entry called with: ", entry_init$id, " / ", entry_init$value)
        eid <- as.character(entry_init$id)
        insertUI(
          selector = paste0("#", ns("entries_container")), where = "beforeEnd",
          ui = div(id = ns(paste0("entry_wrap_", eid)),
                   entryModuleUI(ns(paste0("entry_", eid)), entry_init, rule_state$input_type, metabolite_choices))
        )
        entry_modules[[eid]] <- entryModuleServer(paste0("entry_", eid), entry_init, remove_entry)
        entry_ids(c(entry_ids(), eid))
      }
      
      remove_entry <- function(eid) {
        eid <- as.character(eid)
        if (length(entry_ids()) <= 1) return(invisible())  # mind. 1 Zeile behalten
        removeUI(selector = paste0("#", ns(paste0("entry_wrap_", eid))))
        entry_modules[[eid]] <- NULL
        entry_ids(setdiff(entry_ids(), eid))
      }
      
      observeEvent(input$input_type, {
        if (identical(input$input_type, isolate(rule_state$input_type))) return()
        
        rule_state$input_type <- input$input_type
        old_ids <- entry_ids()
        for (eid in old_ids) {
          removeUI(selector = paste0("#", ns(paste0("entry_wrap_", eid))))
          entry_modules[[eid]] <- NULL
        }
        entry_ids(character(0))
        if (length(old_ids) == 0) old_ids <- "1"
        for (eid in old_ids) add_entry(list(id = as.numeric(eid), value = ""))
      }, ignoreInit = TRUE, ignoreNULL = TRUE)
      
      isolate({ for (e in initial$entries) add_entry(e) })
      
      observeEvent(input$label, rule_state$label <- input$label, ignoreInit = TRUE)
      observeEvent(input$priority, rule_state$priority <- input$priority, ignoreInit = TRUE)
      observeEvent(input$precursor, rule_state$just_after_precursor <- isTRUE(input$precursor), ignoreInit = TRUE)
      observeEvent(input$sumOfSmall, rule_state$sumOfSmall <- isTRUE(input$sumOfSmall), ignoreInit = TRUE)
      observeEvent(input$sequential, rule_state$sequential <- isTRUE(input$sequential), ignoreInit = TRUE)
      
      observeEvent(input$slider, {
        req(input$slider)
        rule_state$min_repeat <- input$slider[1]
        rule_state$max_repeat <- input$slider[2]
      }, ignoreInit = TRUE, priority = 1)
      observeEvent(input$add_entry, {
        eid <- next_entry_id()
        add_entry(list(id = eid, value = ""))
        next_entry_id(eid + 1)
      })
      observeEvent(input$duplicate, on_duplicate(rule_state$id))
      observeEvent(input$delete, on_delete(rule_state$id))
      
      current_entries <- reactive({
        lapply(entry_ids(), function(eid) entry_modules[[eid]]())
      })
      
      output$mass_preview <- renderUI({
        masses <- vapply(current_entries(), function(e) {
          resolve_mass(list(input_type = rule_state$input_type, value = e$value))  # nutzt deine bestehende resolve_mass()
        }, numeric(1))
        masses <- masses[!is.na(masses)]
        if (length(masses) == 0) return(NULL)

        tagList(
          div(
            class = "mass-preview",
            "Resolved mass(es): ",
            createInfoI(
              "Mass_preview_info",
              paste0(
                "<span class='peak-badge positive'>Peak</span> ",
                "values are positive.<br>",
                "<span class='peak-badge negative'>Neutral Loss</span> ",
                "values are negative."
              )
            ),
            lapply(masses, function(m) {
              span(
                class = if (m >= 0) "peak-badge positive" else "peak-badge negative",
                round(m, 4)
              )
            })
          ),
          br()
        )
      })

      
      reactive({
        list(id = rule_state$id, label = rule_state$label, input_type = rule_state$input_type,
             entries = current_entries(), next_entry_id = next_entry_id(),
             min_repeat = rule_state$min_repeat, max_repeat = rule_state$max_repeat,
             just_after_precursor = rule_state$just_after_precursor, sumOfSmall = rule_state$sumOfSmall, 
             sequential = rule_state$sequential, next_card_id = rule_state$next_card_id)
      })
    }, session = parent_session)
  }
  
  card_modules <- reactiveValues()
  card_ids <- reactiveVal(character(0))
  next_card_id <- reactiveVal(1)
  
  add_card_to_ui <- function(card_init) {
    cid <- as.character(card_init$id)
    insertUI(selector = "#cards_container", where = "beforeEnd",
             ui = cardModuleUI(paste0("card_", cid), card_init, METABOLITES$name))
    card_modules[[cid]] <- cardModuleServer(
      paste0("card_", cid), card_init,
      on_delete = delete_card_by_id,
      on_duplicate = duplicate_card_by_id,
      metabolite_choices = METABOLITES$name,
      parent_session = session    
    )
    card_ids(c(card_ids(), cid))
  }
  
  delete_card_by_id <- function(id) {
    cid <- as.character(id)
    removeUI(selector = paste0("#wrap_card_", cid))
    card_modules[[cid]] <- NULL
    card_ids(setdiff(card_ids(), cid))
  }
  
  duplicate_card_by_id <- function(id) {
    cid <- as.character(id)
    src <- card_modules[[cid]]()
    new_id <- next_card_id()
    new_card <- src
    new_card$id <- new_id
    new_card$label <- paste(src$label, "duplicate")
    add_card_to_ui(new_card)
    next_card_id(new_id + 1)
  }

  
  observeEvent(input$add_card, {
    id <- next_card_id()
    add_card_to_ui(empty_card(id))
    next_card_id(id + 1)
  })
  
  

  
  #------------------------------- Main UI ---------------------------------------
  output$main_screen <- renderUI({
    
    #------------------------------- Welcome Screen UI ---------------------------------------
    if (state$screen == "welcome") {
      tagList(
        div(class = "top-bar",
            h2(class = "app-title", "Cera Search Project - Setup"),
            div(id = "load_project_wrapper",
                fileInput("load_project_file", label = NULL, accept = ".zip",
                          buttonLabel = tagList(icon("folder-open"), "Load Project"),
                          placeholder = "")
            )
        ),
        
        fluidRow(
          column(6,
                 div(class = "card required-card",
                     div(class = "required-badge", "Required"),
                     div(class = "card-title", "Required Inputs"),
                     div(
                       style = "margin-left: auto; min-width: 150px;",
                       fileInput("ms2_file", "MS2 Data", placeholder = ifelse(is.null(MS2file_info()),"No file selected", MS2file_info()),
                                 accept = c(".mgf", ".msp", ".mzmine", ".mzHub")),
                       div(
                         style = "text-align: right;",
                         textOutput("MS2file_info", inline = TRUE)
                       )
                     ),
                     
                     selectInput("mode", "Ionization Mode",
                                 choices = c("Positive" = "pos", "Negative" = "neg"),
                                 selected = state$mode),
                     br(),
                     tags$label(class = "control-label", "Output / Project Directory"),
                     div(
                       style = "display:flex; align-items:center;",
                       shinyDirButton("out_dir", "Select output directory","Select output directory"),
                       div(
                         style = "margin-left:auto;",
                         textOutput("out_dir_text", inline = TRUE)
                       )
                     )
                     
                 ),
                actionButton("start_analysis", "Start Analysis",
                             class = "btn-primary btn-lg start-btn")
         
          ),
          
          column(6,
                 div(class = "card optional-card",
                     div(class = "optional-badge", "Optional"),
                     div(class = "card-title", "Optional Inputs"),
                     textInput("project_name", "Project Name",
                               value = state$project_name,
                               placeholder = "My Project"),
                     fileInput("ms1_file", "MS1 Data (CSV)", accept = c(".csv")),
                     
                     fileInput("sirius_file", "Sirius Annotation",
                               accept = c(".csv", ".tsv")),
                     
                     fileInput("rules_file", "Mass Rules (json)", accept = c(".json"))
                    
                 )
          )
        )
      )
    }
    #------------------------------- Analysis & Parameter Screen UI ---------------------------------------
    else if(state$screen == "analysis"){
      tagList(
        div(class = "top-bar",
            h2(class = "app-title", "Cera Search Project - Analysis"),
            div(
              actionButton("save_project", "Save Project",
                           icon = icon("folder-open"), class = "btn-default"),
              actionButton("back_to_welcome", "Return to Setup",
                           icon = icon("house"), class = "btn-default")
            )
        ),
        
        # ----- Legend Bar -----
        
        # Enable Bootstrap tooltips
        tags$script(HTML("
            $(function () {
              $('[data-bs-toggle=\"tooltip\"]').tooltip();
            });
          ")),
        
        div(
          style = "
            display:flex;
            justify-content:center;
            margin-bottom:20px;
          ",
          
          div(class = "analysis_legend",
              tags$span(
                tags$i(style = "color:#28a745;cursor:help;", `data-bs-toggle` = "tooltip", class = "fas fa-square", title = "The filter affects the score."),
                "Basic"
              ),
              tags$span(
                tags$i(style = "color:#fd7e14;cursor:help;", `data-bs-toggle` = "tooltip", class = "fas fa-square", title = "The filter affects the score greatly"),
                "Important"
              ),
              tags$span(
                tags$i(style = "color:#dc3545;cursor:help;", `data-bs-toggle` = "tooltip", class = "fas fa-square", title = "Spectra not included in the filter are removed"),
                "Required"
              )
          )
        ),
        
        # ----- Filters -----
        fluidRow(
          column(
            width = 12,
            lipidDefinitionUI("lipid")
          ),
          column(
            width = 4,
            div(class = "card analysis-card",
                tags$label(
                  class = "switch",
                  tags$input(
                    type = "checkbox",
                    id = "LipidIN_active",
                    checked = if (ToolSettings$LipidIN) "checked" else NULL
                  ),
                  tags$span(class = "slider")
                ),
                div(class = if (!isTRUE(input$LipidIN_active)) "disabled-div" , 
                div(class = "card-title", "LipidIN",
                    createInfoI("LipidIN_info","LipidIN is a tool")),
                createPriorityRadioButtons("LipidIN_priority", ToolSettings$LipidIN_priority),
                tags$details(
                  tags$summary(
                    class = "advanced-settings",
                    "Advanced Settings"
                  ),
                  
                  br(),
                  numericInput("LipidIN_ppm1", "ppm 1", value = ToolSettings$LipidIN_ppm1),
                  numericInput("LipidIN_ppm2", "ppm 2", value = ToolSettings$LipidIN_ppm2)
                )
                
            )
            ),
          ),
          column(
            width = 4,
            div(class = "card analysis-card",
                div(class = "card-title", "Basic filters"
          ),
                div(class = "flex_distance",
                    div(class = "input-block",numericInput("mz_min", "m/z minimum",  basicFilterValues$mz_min)),
                    div(class = "input-block",numericInput("mz_max", "m/z maximum", basicFilterValues$mz_max)),
                    createPriorityRadioButtons("mz_priority", basicFilterValues$mz_priority, horizontal = FALSE)),
                div(class = "flex_distance",
                    div(class = "input-block",numericInput("rt_min", "RT minimum (in seconds)", basicFilterValues$rt_min)),
                    div(class = "input-block",numericInput("rt_max", "RT maximum (in seconds)", basicFilterValues$rt_max)),
                    createPriorityRadioButtons("rt_priority",basicFilterValues$rt_priority, horizontal = FALSE)),
                div(class = "flex_distance",
                    div(class = "input-block", numericInput("n_peaks_min", "minimum number of peaks", basicFilterValues$n_peaks_min)),
                    div(class = "input-block", numericInput("n_peaks_max", "maximum number of peaks", basicFilterValues$n_peaks_max)),
                    createPriorityRadioButtons("n_peaks_priority",basicFilterValues$n_peaks_priority, horizontal = FALSE))
            ) 
          ),
          column(
            width = 4,
            div(class = "card analysis-card",
                tags$label(
                  class = "switch",
                  tags$input(
                    type = "checkbox",
                    id = "DreaMS_active",
                    checked = if (ToolSettings$DreaMS) "checked" else NULL
                  ),
                  tags$span(class = "slider")
                ),
                div(class = if (!isTRUE(input$DreaMS_active)) "disabled-div" , 
                div(class = "card-title", "DreaMS",
                    createInfoI("DreaMS_info","DreaMS is a tool")),
                
                createPriorityRadioButtons("DreaMS_priority", ToolSettings$DreaMS_priority),
                
                tags$details(
                  tags$summary(
                    class = "advanced-settings",
                    "Advanced Settings"
                  ),
                  
                  br(),
                  numericInput("ppm1", "ppm 1", value = 10),
                  numericInput("ppm2", "ppm 2", value = 50)
                )
                )
            )
          )
        ),

        
        fluidRow(
          column(
            width = 12,
            titlePanel("MS2 Peak Filter Builder"),

            div(class = "flex_distance",
              actionButton("add_card", "Add peak filter", icon = icon("plus"), class = "btn-primary"),
              selectInput("card_sort", "Order by",
                                        choices = c("Creation time" = "creation", "Mass" = "mass"),
                                        selected = "creation",
                          , width = "140px"),
              numericInput("tolerance_ppm", "Tolerance (ppm)", value = basicFilterValues$tolerance_ppm, width = "140px"),
              numericInput("minimum_dalton", "Minimum Dalton", value = basicFilterValues$minimum_dalton, width = "140px")),
            
            hr(),
            div(id = "cards_container", class = "cards-container")
            
          )
        ),
        fluidRow(
          column(12, align = "center",
                 actionButton("run_analysis", "Run Analysis",
                              class = "btn-primary btn-lg start-btn")
          )
        )
      )
    }
    #------------------------------- Results Screen UI ---------------------------------------
    else if(state$screen == "results"){
      tagList(
        div(class = "top-bar",
            h2(class = "app-title", "Cera Search Project - Results"),
            div(
              shinyDirButton("save_project", "Save Project",
                             title = "Save current project",
                             icon = icon("folder-open"), class = "btn-default"),
              actionButton("back_to_welcome", "Return to Setup",
                           icon = icon("house"), class = "btn-default")
            )
          ),
        textOutput("result_text", inline = TRUE)
      )
    }
  })
}
