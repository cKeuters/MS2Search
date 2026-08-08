#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#' @import shiny
#' @import shinyFiles
#' @import shinyWidgets
#' @import jsonlite
#' @import zip
#' @noRd
app_server <- function(input, output, session) {
  
  #-----------------------------------------------------------------------------
  # Reactive Values
  #-----------------------------------------------------------------------------
  
  # Application state values
  #
  # Stores global application state information used for navigation
  # and project configuration within the Shiny application.
  #
  # @details
  # The following values are maintained:
  # \describe{
  #   \item{screen}{Current application screen.}
  #   \item{out_dir}{Output directory for generated files.}
  #   \item{mode}{Ionization mode (e.g. `"pos"` or `"neg"`).}
  #   \item{project_name}{Name of the current project.}
  # }
  # @keywords internal
  state <- reactiveValues(
    screen = "welcome",
    out_dir = path.expand("~"),
    mode = "pos",
    project_name = ""
  )
  
  # Values of Basisfilters
  #
  # Contains the initial and reactive values for the basisfilters of
  # mz, rt, peak counts and tolerance.
  #
  # @details
  # Variables and values:
  # \describe{
  #   \item{mz_min}{Minimum m/z (initial 0)}
  #   \item{mz_max}{Maximum m/z (initial 1000)}
  #   \item{mz_priority}{priority of m/z-filter}
  #   \item{rt_min}{Minimal retention time in seconds (initial 0)}
  #   \item{rt_max}{Maximal retention time in seconds (initial 1000)}
  #   \item{rt_priority}{priority of the retention time filter}
  #   \item{n_peaks_min}{Minimum number of peaks in MS2 spectrum (initial 0)}
  #   \item{n_peaks_max}{Maximum number of peaks in MS2 spectrum (initial 1000)}
  #   \item{n_peaks_priority}{priority of peak filter}
  #   \item{tolerance_ppm}{Mass tolerance in ppm (initial 10)}
  #   \item{minimum_dalton}{Minimal absolute mass difference in dalton (initial 0.01)}
  # }
  # @keywords internal
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
  
  # Tool configuration settings
  #
  # Stores user-configurable settings and activation states for
  # external analysis tools used by the application.
  #
  # @details
  # The following settings are available:
  # \describe{
  #   \item{LipidIN}{Logical flag indicating whether LipidIN analysis is enabled.}
  #   \item{LipidIN_ppm1}{Primary mass tolerance (ppm) used by LipidIN.}
  #   \item{LipidIN_ppm2}{Secondary mass tolerance (ppm) used by LipidIN.}
  #   \item{LipidIN_priority}{Execution priority assigned to LipidIN.}
  #   \item{DreaMS}{Logical flag indicating whether DreaMS analysis is enabled.}
  #   \item{DreaMS_priority}{Execution priority assigned to DreaMS.}
  # }
  # @keywords internal
  ToolSettings <- reactiveValues(
    LipidIN = TRUE,
    LipidIN_ppm1 = 10,
    LipidIN_ppm2 = 50,
    LipidIN_priority = "basic",
    DreaMS = TRUE,
    DreaMS_priority = "basic"
  )
  
  analysis_results <- reactiveValues(
    spectra = NULL,
    results = NULL,
    running = FALSE
  )
  
  #-----------------------------------------------------------------------------
  # ObserveEvent of Reactive Variables
  #-----------------------------------------------------------------------------
  
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
  
  # states
  observeEvent(input$mode,         { state$mode <- input$mode })
  observeEvent(input$project_name, { state$project_name <- input$project_name })
  
  #-----------------------------------------------------------------------------
  # ObserveEvent of Outputs and Buttons
  #-----------------------------------------------------------------------------
  
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
  
  # Output MS2 file name
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
  
  # Button Start Analysis
  observeEvent(input$start_analysis, {
    if(input$project_name == "") state$project_name <- paste0("project_",format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
    state$screen <- "analysis"
  })
  
  # Button Return to Setup
  observeEvent(input$back_to_welcome, {
    add_log("----------- Return to Setup -------------")
    log_text <- ""
    state$screen <- "welcome"
  })
  
  #Button run analysis
  observeEvent(input$run_analysis, {
     
     state$screen <- "results"
   })
  
  
  #-----------------------------------------------------------------------------
  # Load and Save Project
  #-----------------------------------------------------------------------------
  
  #Select Directionary
  roots <- c(
    home = normalizePath(path.expand("~"), winslash = "/"),
    mount = "/mnt/ifs/data/Workspace"
  )
  shinyDirChoose(input, "out_dir", roots = roots, session = session)
  shinyDirChoose(input, "load_project", roots = roots, session = session)
  
  # Save Project
  observeEvent(input$save_project, {
    setup_and_save_project(create_settings_list(
      state = reactiveValuesToList(state),
      basicFilterValues = reactiveValuesToList(basicFilterValues),
      toolSettings = reactiveValuesToList(ToolSettings),
      cards   = all_cards_data()),
      paste0(state$out_dir, "/", state$project_name,".zip")
    )
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
  })
  
  #-----------------------------------------------------------------------------
  # Rule based interactions
  #-----------------------------------------------------------------------------
  
  card_modules <- reactiveValues()
  card_ids <- reactiveVal(character(0))
  next_card_id <- reactiveVal(1)
  
  all_cards_data <- reactive({
    ids <- card_ids()
    if (length(ids) == 0) return(list())
    lapply(ids, function(cid) card_modules[[cid]]())
  })
  
  observeEvent(input$add_card, {
    id <- next_card_id()
    add_card_to_ui(empty_card(id))
    next_card_id(id + 1)
  })
  
  
  add_card_to_ui <- function(card_init) {
    cid <- as.character(card_init$id)
    card_modules[[cid]] <- mod_Rules_server(
      paste0("card_", cid), card_init, state$mode,
      on_delete = delete_card_by_id,
      on_duplicate = duplicate_card_by_id,
      metabolite_choices = METABOLITES$name,
      parent_session = session    
    )
    card_ids(c(card_ids(), cid))
    
    insertUI(
      selector = "#cards_container",
      where = "beforeEnd",
      ui = div(id = paste0("wrap_card_", cid),
               mod_Rules_ui(paste0("card_", cid), card_init, METABOLITES$name))
    )
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
  
  # # UI rendern
  # output$cards_ui <- renderUI({
  #   ids <- card_ids()
  #   if (length(ids) == 0) return(NULL)
  #   
  #   lapply(ids, function(cid) {
  #     initial_card <- isolate(card_modules[[cid]]())
  #     mod_Rules_ui(paste0("card_", initial_card$id), initial_card, METABOLITES$name)
  #   })
  # })
  
  
  #-----------------------------------------------------------------------------
  # Results
  #-----------------------------------------------------------------------------
  
  output$log <- renderText(log_text())

  spectra <- eventReactive(input$run_analysis, {
    load_spectra(input$ms2_file$datapath)
  })
  
  results <- eventReactive(input$run_analysis, {
    settings <- create_settings_list(
      state = reactiveValuesToList(state),
      basicFilterValues = reactiveValuesToList(basicFilterValues),
      toolSettings = reactiveValuesToList(ToolSettings),
      cards   = all_cards_data())
    run_analysis(settings, spectra())
  })
  
  selected_feature_id <- reactive({
    req(nrow(results()) > 0)
    req(input$featureTable_rows_selected)
    res <- results()
    res$FEATURE_ID[input$featureTable_rows_selected]
  })
  
  output$spectrum_plot <- renderPlot({
    req(nrow(results()) > 0)
    req(selected_feature_id())
    plot_spectrum_by_id(spectra(), selected_feature_id(), id_column = "FEATURE_ID")
  })


  #-----------------------------------------------------------------------------
  # Screens hide/show
  #-----------------------------------------------------------------------------
  
  observeEvent(state$screen, {
    shinyjs::hide("welcome_screen")
    shinyjs::hide("cards_section")
    shinyjs::hide("analysis_footer")
    shinyjs::hide("results_screen")
    shinyjs::hide("analysis_screen")
    
    message(paste("Screen:", state$screen))
    
    if (state$screen == "welcome") {
      shinyjs::show("welcome_screen")
    } else if (state$screen == "analysis") {
      shinyjs::show("cards_section")
      shinyjs::show("analysis_footer")
      shinyjs::show("analysis_screen")
    } else if (state$screen == "results") {
      shinyjs::show("results_screen")
    }
  }, ignoreNULL = FALSE)
  
  
  #-----------------------------------------------------------------------------
  # Screens Outputs
  #-----------------------------------------------------------------------------
  
  output$welcome_ui <- renderUI({ 
    tagList(
      create_topbar_InSetup(),
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
    
    
    
    })
  output$analysis_header_ui <- renderUI({ 
    
    tagList(
      create_topbar_InProject(),
      # Enable Bootstrap tooltips
      tags$script(HTML("
            $(function () {
              $('[data-bs-toggle=\"tooltip\"]').tooltip();
            });
          ")),
      createAnalysisLegend(),
      fluidRow(
        #   column(
        #     width = 12,
        #     lipidDefinitionUI("lipid")
        #   ),
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
              div(class = "flex-distance",
                  div(class = "input-block",numericInput("mz_min", "m/z minimum",  basicFilterValues$mz_min)),
                  div(class = "input-block",numericInput("mz_max", "m/z maximum", basicFilterValues$mz_max)),
                  createPriorityRadioButtons("mz_priority", basicFilterValues$mz_priority, horizontal = FALSE)),
              div(class = "flex-distance",
                  div(class = "input-block",numericInput("rt_min", "RT minimum (in seconds)", basicFilterValues$rt_min)),
                  div(class = "input-block",numericInput("rt_max", "RT maximum (in seconds)", basicFilterValues$rt_max)),
                  createPriorityRadioButtons("rt_priority",basicFilterValues$rt_priority, horizontal = FALSE)),
              div(class = "flex-distance",
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
      )
    )
    
    })
  
  output$rule_coordination <- renderUI({ 
    div(
  titlePanel("MS2 Peak Filter Builder"),
  div(style = "display:flex; gap:10px; align-items:center; margin-bottom:10px;",
      actionButton("add_card", "Add peak filter", icon = icon("plus"), class = "btn-primary"),
      selectInput("card_sort", "Order by",
                  choices = c("Creation time" = "creation", "Mass" = "mass"),
                  selected = "creation", width = "140px"),
      numericInput("tolerance_ppm", "Tolerance (ppm)", value = basicFilterValues$tolerance_ppm, width = "140px"),
      numericInput("minimum_dalton", "Minimum Dalton", value = basicFilterValues$minimum_dalton, width = "140px")
  ),
  hr(),
  div(id="cards_container", class = "cards-container"))
})
  
  output$results_ui <- renderUI({ 
    tagList(
      create_topbar_InProject(),
      verbatimTextOutput("log")
    )
    })
  
  output$featureTable <- DT::renderDataTable({
    message("output$featureTable")
    req(results())
    res <- as.data.frame(results())
    
    DT::datatable(
      res[, c("FEATURE_ID", "precursorMz", "rtime", "score")],
      selection = "single",
      class = "feature-table",
      rownames = FALSE,
      colnames = c("Feature ID", "m/z", "Retention time", "Score"),
      options = list(
        pageLength = 25,
        lengthMenu = c(25, 50, 100),
        order = list(list(3, "desc"))
      )
    ) |>
      DT::formatRound(columns = "precursorMz", digits = 6) |>
      DT::formatRound(columns = "rtime", digits = 4) |>
      DT::formatRound(columns = "score", digits = 3)
    
  }, server = TRUE)
  
  
  # feature details
    output$spectrumDetails <- renderUI({
      if (is.null(selected_feature_id())) {
        return(div(class = "placeholder-text", "Select Feature"))
      }
      req(nrow(results()) > 0)
      res <- results()
      i <- as.numeric(selected_feature_id())
      
      tagList(
        plotOutput("spectrum_plot"),
        div(class = "detail-row",
            div(class = "detail-label", "Feature ID"),
            div(class = "detail-value", as.character(res$FEATURE_ID[i]))
        ),
        div(class = "detail-row",
            div(class = "detail-label", "m/z"),
            div(class = "detail-value", as.character(format(res$precursorMz[i], digits = 6)))
        ),
        div(class = "detail-row",
            div(class = "detail-label", "Retention time"),
            div(class = "detail-value", as.character(format(res$rtime[i], digits = 4)))
        ),
        div(class = "detail-row",
            div(class = "detail-label", "Score"),
            div(class = "detail-value", as.character(format(res$score[i], digits = 3)))
        ),
        div(class = "detail-row",
            div(class = "detail-label", "Adduct"),
            div(class = "detail-value", res$ION[i])
        ),
        div(class = "detail-row",
            div(class = "detail-label", "Polarity"),
            div(class = "detail-value", res$POLARITY[i])
        )
      )
    })



  

  
}

