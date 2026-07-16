library(shiny)
library(shinyFiles)
library(shinyWidgets)
library(bslib)
library(qs) #loading and saving data

options(shiny.maxRequestSize = 500 * 1024^2)
#shiny::runApp( appDir = "CeraSearch/CeraSearch_Test/", launch.browser = TRUE)

source("ui_helpers.R")

#-------------------------------------------------------------------------------
# Main Server
#-------------------------------------------------------------------------------
function(input, output, session) {
  state <- reactiveValues(
    screen = "welcome",
    out_dir = path.expand("~"),
    mode = "pos",
    project_name = "",
    MS2LDA = TRUE,
    DreaMS = TRUE
  )
  
  rv <- reactiveValues(
    cards = list())
  
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
  observeEvent(input$load_project, {
    req(input$load_project)
    path <- as.character(parseDirPath(roots, input$load_project))
    if (length(path) > 0 && nzchar(path)) {
      state$out_dir <- path
      state$project_name <- basename(path)
      showNotification(paste("Loaded project:", path), type = "message")
    }
  })
  
  # States
  observeEvent(input$mode,         { state$mode <- input$mode })
  observeEvent(input$project_name, { state$project_name <- input$project_name })
  
  # Start Analysis
  observeEvent(input$start_analysis, {
    state$screen <- "analysis"
  })
  
  observeEvent(input$back_to_welcome, {
    state$screen <- "welcome"
  })
  
  observeEvent(input$run_analysis, {
    save_parameters(
      state = reactiveValuesToList(state),
      rv    = reactiveValuesToList(rv),
      ppm1  = input$ppm1,
      ppm2  = input$ppm2,
      project_path = state$out_dir
    )
    
    state$screen <- "results"
  })
  
  observeEvent(input$save_project, {
    save_parameters(
      state = reactiveValuesToList(state),
      rv    = reactiveValuesToList(rv),
      ppm1  = input$ppm1,
      ppm2  = input$ppm2,
      project_path = state$out_dir
    )
  })
  
  observeEvent(input$MS2LDA_active, {
    state$MS2LDA <- input$MS2LDA_active
  })
  observeEvent(input$DreaMS_active, {
    state$DreaMS <- input$DreaMS_active
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
        actionButton(ns("remove"), NULL, icon = icon("times"), class = "btn-sm btn-outline-danger")
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
  
  # Sort cards reactive helper
  sorted_cards <- reactive({
    sort_mode <- if (!is.null(input$card_sort)) input$card_sort else "creation"
    cards <- rv$cards
    if (sort_mode == "mass") {
      masses <- sapply(cards, resolve_mass)
      ord <- order(masses, na.last = TRUE)
      cards <- cards[ord]
    }
    cards
  })
  
  
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
                createPriorityRadioButtons(ns("priority")),
                noUiSliderInput(
                  inputId = ns("slider"), label = "Number of losses",
                  min = 0, max = 10, value = c(card$min_repeat, card$max_repeat), step = 1,
                  format = wNumbFormat(decimals = 0), tooltips = TRUE,
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
                checkboxInput(ns("sumOfSmall"), "Validation by sum of smaller peaks", value = card$sumOfSmall)
                
            )
        )
    )
  }
  
  cardModuleServer <- function(id, initial, on_delete, on_duplicate, metabolite_choices, parent_session ) {
    moduleServer(id, function(input, output, session) {
      ns <- session$ns
      
      state <- reactiveValues(
        id = initial$id, label = initial$label, input_type = initial$input_type,
        min_repeat = initial$min_repeat, max_repeat = initial$max_repeat,
        just_after_precursor = initial$just_after_precursor, sumOfSmall = initial$sumOfSmall, next_card_id = initial$next_card_id
      )
      
      entry_modules <- reactiveValues()
      entry_ids <- reactiveVal(character(0))
      next_entry_id <- reactiveVal(initial$next_entry_id)
      
      add_entry <- function(entry_init) {
        message("add_entry called with: ", entry_init$id, " / ", entry_init$value)
        eid <- as.character(entry_init$id)
        insertUI(
          selector = paste0("#", ns("entries_container")), where = "beforeEnd",
          ui = div(id = ns(paste0("entry_wrap_", eid)),
                   entryModuleUI(ns(paste0("entry_", eid)), entry_init, state$input_type, metabolite_choices))
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
        if (identical(input$input_type, isolate(state$input_type))) return()
        
        state$input_type <- input$input_type
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
      
      observeEvent(input$label, state$label <- input$label, ignoreInit = TRUE)
      observeEvent(input$priority, state$priority <- input$priority, ignoreInit = TRUE)
      observeEvent(input$precursor, state$just_after_precursor <- isTRUE(input$precursor), ignoreInit = TRUE)
      observeEvent(input$sumOfSmall, state$sumOfSmall <- isTRUE(input$sumOfSmall), ignoreInit = TRUE)
      
      observeEvent(input$slider, {
        req(input$slider)
        state$min_repeat <- input$slider[1]
        state$max_repeat <- input$slider[2]
      }, ignoreInit = TRUE, priority = 1)
      observeEvent(input$add_entry, {
        eid <- next_entry_id()
        add_entry(list(id = eid, value = ""))
        next_entry_id(eid + 1)
      })
      observeEvent(input$duplicate, on_duplicate(state$id))
      observeEvent(input$delete, on_delete(state$id))
      
      current_entries <- reactive({
        lapply(entry_ids(), function(eid) entry_modules[[eid]]())
      })
      
      output$mass_preview <- renderUI({
        masses <- vapply(current_entries(), function(e) {
          resolve_mass(list(input_type = state$input_type, value = e$value))  # nutzt deine bestehende resolve_mass()
        }, numeric(1))
        masses <- masses[!is.na(masses)]
        if (length(masses) == 0) return(NULL)

        tagList(
          div(class = "mass-preview", HTML(paste("Resolved mass(es):", paste(round(masses, 4), collapse = ", ")))),
          br()
        )
      })

      
      reactive({
        list(id = state$id, label = state$label, input_type = state$input_type,
             entries = current_entries(), next_entry_id = next_entry_id(),
             min_repeat = state$min_repeat, max_repeat = state$max_repeat,
             just_after_precursor = state$just_after_precursor, sumOfSmall = state$sumOfSmall, next_card_id = state$next_card_id)
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
            shinyDirButton("load_project", "Load Project",
                           title = "Select existing project folder",
                           icon = icon("folder-open"), class = "btn-default")
        ),
        
        fluidRow(
          column(6,
                 div(class = "card required-card",
                     div(class = "required-badge", "Required"),
                     div(class = "card-title", "Required Inputs"),
                     
                     fileInput("ms2_file", "MS2 Data",
                               accept = c(".mgf", ".msp", ".mzmine")),
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
                     
                 )
          ),
          
          column(6,
                 div(class = "card optional-card",
                     div(class = "optional-badge", "Optional"),
                     div(class = "card-title", "Optional Inputs"),
                     
                     fileInput("ms1_file", "MS1 Data (CSV)", accept = c(".csv")),
                     
                     fileInput("sirius_file", "Sirius Annotation",
                               accept = c(".csv", ".tsv")),
                     
                     fileInput("rules_file", "Mass Rules (CSV)", accept = c(".csv")),
                     
                     textInput("project_name", "Project Name",
                               value = state$project_name,
                               placeholder = "My Project")
                 )
          )
        ),
        
        fluidRow(
          column(12, align = "center",
                 actionButton("start_analysis", "Start Analysis",
                              class = "btn-primary btn-lg start-btn")
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
              shinyDirButton("save_project", "Save Project",
                             title = "Save current project",
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
            width = 4,
            div(class = "card analysis-card",
                tags$label(
                  class = "switch",
                  tags$input(
                    type = "checkbox",
                    id = "MS2LDA_active",
                    checked = if (state$MS2LDA) "checked" else NULL
                  ),
                  tags$span(class = "slider")
                ),
                div(class = if (!isTRUE(input$MS2LDA_active)) "disabled-div" , 
                div(class = "card-title", "MS2LIPID",
                    createInfoI("MS2LIPID_info","MS2LIPID is a tool")),
                createPriorityRadioButtons("MS2Lipid_priority"),
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
            ),
          ),
          column(
            width = 4,
            div(class = "card analysis-card",
                div(class = "card-title", "Basic filters"
          ),
                div(class = "flex_distance",
                    div(class = "input-block",numericInput("mz_min", "m/z minimum", 0)),
                    div(class = "input-block",numericInput("mz_max", "m/z maximum", 10000)),
                    createPriorityRadioButtons("mz_priority", horizontal = FALSE)),
                div(class = "flex_distance",
                    div(class = "input-block",numericInput("rt_min", "RT minimum", 0)),
                    div(class = "input-block",numericInput("rt_max", "RT maximum", 10000)),
                    createPriorityRadioButtons("rt_priority", horizontal = FALSE)),
                div(class = "flex_distance",
                    div(class = "input-block", numericInput("peakCount_min", "minimum number of peaks", 0)),
                    div(class = "input-block", numericInput("peakCount_max", "maximum number of peaks", 10000)),
                    createPriorityRadioButtons("peakCount_priority", horizontal = FALSE))
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
                    checked = if (state$DreaMS) "checked" else NULL
                  ),
                  tags$span(class = "slider")
                ),
                div(class = if (!isTRUE(input$DreaMS_active)) "disabled-div" , 
                div(class = "card-title", "DreaMS",
                    createInfoI("DreaMS_info","DreaMS is a tool")),
                
                createPriorityRadioButtons("DreaMS_priority"),
                
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
                                        selected = "creation")),
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
            shinyDirButton("load_project", "Load Project",
                           title = "Select existing project folder",
                           icon = icon("folder-open"), class = "btn-default")
        )
      )
    }
  })
}
