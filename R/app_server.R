#plan(multisession)
#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#' @import shiny
#' @import shinyFiles
#' @import shinyWidgets
#' @import jsonlite
#' @import zip
#' @import future
#' @import promises
#' @import readr
#' @import shinyTree
#' @import rapiclient
#' @import data.tree
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
    screen = "Setup",
    out_dir = path.expand("~"),
    mode = "pos",
    project_name = "",
    newMS2data = TRUE
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
  # external Analysis tools used by the application.
  #
  # @details
  # The following settings are available:
  # \describe{
  #   \item{LipidIN}{Logical flag indicating whether LipidIN Analysis is enabled.}
  #   \item{LipidIN_ppm1}{Primary mass tolerance (ppm) used by LipidIN.}
  #   \item{LipidIN_ppm2}{Secondary mass tolerance (ppm) used by LipidIN.}
  #   \item{LipidIN_priority}{Execution priority assigned to LipidIN.}
  #   \item{DreaMS}{Logical flag indicating whether DreaMS Analysis is enabled.}
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
  
  # Application results values
  #
  # Stores the spectra object and the analysis results
  #
  # @details
  # The following values are maintained:
  # \describe{
  #   \item{spectra}{MS2 Spectra Object}
  #   \item{results}{dataframe of the passed filters}
  #   \item{running}{state of the analysis}
  #   \item{log}{logs of the analysis}
  # }
  # @keywords internal
  analysis_results <- reactiveValues(
    spectra = NULL,
    results = NULL,
    running = FALSE,
    log = NULL
  )
  
  
  log_state <- reactiveValues(file = NULL)
  log_text <- reactiveVal("")
  
  brkrrba_save <- reactiveVal("")
  head_group <- reactiveVal("")
  
  ms1_data <- reactive({
    read_csv("../../data/example.csv")
    req(input$ms1_file)
    read_csv(input$ms1_file$datapath)
  })
  
  metaboAuth <- reactiveValues(
    authenticated = FALSE, 
    pending_action = NULL
    )
  
  mod_MetaboScape_auth_server("auth", metaboAuth)
  MetaboScapeProjects <- reactiveVal(list())
  MetaboSelectedProjectID <- reactiveVal("")
  
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
  
  observeEvent(input$extract_ms2, {
    metaboAuth$require_auth(function() {
      selectProject_from_metaboscape()   
    })
  })
  
  # show output path
  # observeEvent(input$out_dir, {
  #   req(input$out_dir)
  #   path <- as.character(parseDirPath(roots, input$out_dir))
  #   if (length(path) > 0 && nzchar(path)) {
  #     state$out_dir <- path
  #   }
  # })
  
  # Output MS2 file name
  MS2file_info <- reactive({
    if (is.null(input$ms2_file)) {
      NULL
    } else {
      state$newMS2data <- TRUE
      ext <- tools::file_ext(input$ms2_file$name)
      MS2file_name(switch(ext,
               mgf = input$ms2_file$name,
               validate("Invalid file; Please upload a .mgf file (other file types are following)")
        ))
      input$ms2_file$name
    }
  })
  
  output$MS2file_info <- renderText({
    MS2file_name()
  })
  
  MS2file_name <- reactiveVal("")
  
  # MS2file_name <- reactive({
  #   req(input$ms2_file)
  #   
  #   ext <- tools::file_ext(input$ms2_file$name)
  #   switch(ext,
  #          mgf = input$ms2_file$name,
  #          validate("Invalid file; Please upload a .mgf file (other file types are following)")
  #   )
  # })
  
  # Button Start Analysis
  observeEvent(input$start_analysis, {
    if(input$project_name == "") state$project_name <- paste0("project_",format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
    state$screen <- "Analysis"
  })
  
  # Button Return to Setup
  observeEvent(input$back_to_welcome, {
    state$screen <- "Setup"
  })
  observeEvent(input$back_to_results, {
    state$screen <- "Results"
  })
  observeEvent(input$back_to_analysis, {
    state$screen <- "Analysis"
  })
  #Button run Analysis
  observeEvent(input$run_analysis, {
     state$screen <- "Results"
   })
  
 
  output$export_results <- downloadHandler(
    filename = function() {
      paste0(state$project_name, "_results.csv")
    },
    content = function(file) {
      
      res <- analysis_results$results
      write.csv(res, file, row.names = FALSE)
    }
  )
  
  output$export_brkrrba <- downloadHandler(
    filename = function() {
      paste0(state$project_name, ".brkrrba")
    },
    content = function(file) {
      
      settings_list <- create_settings_list(
        state = reactiveValuesToList(state),
        basicFilterValues = reactiveValuesToList(basicFilterValues),
        toolSettings = reactiveValuesToList(ToolSettings),
        cards   = all_cards_data(), 
        lipid_definition = all_chain_data(),
        headgroup = head_group())
      create_brkrrba_file(settings_list, file)
    }
  )
  
  
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
  # observeEvent(input$save_project, {
  #   setup_and_save_project(create_settings_list(
  #     state = reactiveValuesToList(state),
  #     basicFilterValues = reactiveValuesToList(basicFilterValues),
  #     toolSettings = reactiveValuesToList(ToolSettings),
  #     cards   = all_cards_data(), lipid_definition = all_chain_data(),
  #     headgroup = head_group()),
  #     analysis_results$results,
  #     analysis_results$spectra,
  #     paste0(state$out_dir, "/", state$project_name,".zip")
  #   )
  #   showNotification(paste0("Project saved to ",  paste0(state$out_dir, "/", state$project_name,".zip")))
  # })
  
  output$export_project <- downloadHandler(
    filename = function() {
      paste0(state$project_name, ".zip")
    },
    content = function(file) {
      save <- setup_and_save_project(create_settings_list(
        state = reactiveValuesToList(state),
        basicFilterValues = reactiveValuesToList(basicFilterValues),
        toolSettings = reactiveValuesToList(ToolSettings),
        cards   = all_cards_data(), lipid_definition = all_chain_data(),
        headgroup = head_group()),
        analysis_results$results,
        analysis_results$spectra,
        file)
    }
  )
  

  # Load Project
  observeEvent(input$load_project_file, {
    loaded_project <- load_project(input$load_project_file$datapath)
    settings_list <- loaded_project$settings_list
    project_results <- loaded_project$results
    project_spectra <- loaded_project$spectra
    
    for(id in card_ids()) delete_card_by_id(id)
    
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
    
    head_group(settings_list$chain_info$headgroup)

    for(rule in settings_list$specific_filters){
      add_card_to_ui(rule)
    }
    
    if(!is.null(project_results) & !is.null(project_spectra)){
      analysis_results$results <- project_results
      analysis_results$spectra <- project_spectra
      state$newMS2data <- FALSE
      state$screen <- "Results"
    }
    
    lipid_load_trigger(settings_list$lipid_definition)
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
  
  #-----------------------------------------------------------------------------
  # Lipid Design
  #-----------------------------------------------------------------------------
  
  lipid_load_trigger <- reactiveVal(NULL)
  
  lipidDefinitionServer <- mod_LipidDefinition_server("LipidDefinition_1", lipid_load_trigger, head_group)
  
  all_chain_data <- reactive({
    lipidDefinitionServer()

  })
  
  #-----------------------------------------------------------------------------
  # Results
  #-----------------------------------------------------------------------------
  
  output$log <- renderText(log_text())
  
  output$state_screen <- renderUI({
    h2(class = "app-title", paste0("MS2Search - ", state$screen))
  })
  
  selected_feature_id <- reactive({
    req(nrow(analysis_results$results) > 0)
    req(input$featureTable_rows_selected)
    res <- analysis_results$results
    res$FEATURE_ID[input$featureTable_rows_selected]
  })
  
  output$interactive_spectrum <- renderPlotly({
    req(nrow(analysis_results$results) > 0)
    req(selected_feature_id())
    settings_list <- create_settings_list(
      state = reactiveValuesToList(state),
      basicFilterValues = reactiveValuesToList(basicFilterValues),
      toolSettings = reactiveValuesToList(ToolSettings),
      cards   = all_cards_data(), 
      lipid_definition = all_chain_data(),
      headgroup = head_group())
    plottingSpectra(analysis_results$spectra, selected_feature_id(), settings_list)
    
  })
  
  output$spectrum_plot <- renderPlot({
    req(nrow(analysis_results$results) > 0)
    req(selected_feature_id())
    par(mfrow = c(2, 1))
    plot_spectrum_by_id(analysis_results$spectra, noise = FALSE, selected_feature_id(), id_column = "FEATURE_ID")
  })
  
  #run Analysis with logging
  observeEvent(input$run_analysis, {
    ms2_path <- input$ms2_file$datapath
    SelectedProjectID <- MetaboSelectedProjectID()
    
    if(state$newMS2data && is.null(ms2_path) && SelectedProjectID == "") showNotification("No MS2 data selected", type= "error")
    req(!is.null(ms2_path) || SelectedProjectID != "" || !state$newMS2data)
    
    state$screen <- "Results"
    analysis_results$running <- TRUE
    
    append_log <- function(msg) {
      current <- log_text()
      log_text(paste0(current, format(Sys.time(), "[%H:%M:%S] "), msg, "\n"))
    }
    
    #append_log("Starting Analysis")
    

    settings <- create_settings_list(
      state             = reactiveValuesToList(state),
      basicFilterValues = reactiveValuesToList(basicFilterValues),
      toolSettings      = reactiveValuesToList(ToolSettings),
      cards             = all_cards_data(),
      lipid_definition   = all_chain_data(),
      headgroup         = head_group()
    )
    
    log_file <- tempfile(paste0("analysis_log_", state$project_name), fileext = ".txt")
    writeLines("", log_file)
    log_state$file <- log_file
    
    
    MetaboScapeAPI <- getMetaboScapeConnection()
    future({
      
      write_log <- function(msg) {
        message(msg)
        cat(
          format(Sys.time(), "[%H:%M:%S]"), msg, "\n",
          file = log_file, append = TRUE
        )
      }

      
      tryCatch({
        
        if(state$newMS2data){
          if(!is.null(ms2_path)) spec <- load_spectra(ms2_path, add_log = write_log)
          else if(SelectedProjectID != ""){
            write_log(paste("Loading with MetaboSelectedProjectID", SelectedProjectID))
            featureTable <- createSpectradf(operations = MetaboScapeAPI, featureTableId = SelectedProjectID )
            spec <- Spectra(featureTable)
          } 
          else write_log("NO MS2")
        }
        else spec <- analysis_results$spectra

        res <- run_analysis(settings, spec, add_log = write_log)
        
        list(
          spectra = spec,
          results = res,
          log_file = log_file
        )
        
      }, error = function(e) {
        write_log("ERROR:")
        write_log(conditionMessage(e))
        
        stop(e)
      })
      
    }) %...>% {
      # success handler
      logs <- paste(readLines(.$log_file, warn = FALSE), collapse = "\n")
      state$newMS2data <- FALSE
      
      analysis_results$spectra <- .$spectra
      analysis_results$results <- .$results
      analysis_results$running <- FALSE
      analysis_results$log <- logs
      
    } %...!% {
      # error handler
      analysis_results$running <- FALSE
      
      message("Future failed:")
      message(.)
      
      if (file.exists(log_file)) {
        analysis_results$log <-
          paste(readLines(log_file, warn = FALSE), collapse = "\n")
      }
    }
    
  })
  
  
  observe({
    req(analysis_results$running)
    req(log_state$file)
    
    invalidateLater(500, session)
    
    if (file.exists(log_state$file)) {
      logs <- paste(readLines(log_state$file, warn = FALSE), collapse = "\n")
      log_text(logs)
    }
  })
  
  
  output$has_results <- reactive({
    nrow(analysis_results$results) > 0
  })
  outputOptions(output, "has_results", suspendWhenHidden = FALSE)
  
  
  #-------------------------------------------------------------------------------
  # MetaboScape API
  #-------------------------------------------------------------------------------
  
  selectProject_from_metaboscape <- function(){
    if(length(MetaboScapeProjects()) == 0){
      
      MetaboScapeAPI <- getMetaboScapeConnection()
      res <- MetaboScapeAPI$listAllProjects()
      projects <- fromJSON(httr::content(res, as = "text", encoding = "UTF-8"))
      
      tree <- list()

      for (i in seq_len(nrow(projects))) {
        project_name <- projects$name[i]
        exps <- projects$experiments[[i]]
        exp_nodes <- list()

        for (j in seq_len(nrow(exps))) {
          ft <- exps$featureTables[[j]]

          ft_nodes <- setNames(
            lapply(seq_len(nrow(ft)), function(k) {
              structure("", sticon = "table", stId = ft$id[k], a_attr = list(title = sprintf("Created: %s", ft$creationDate[k])))
            }),
            paste0(ft$name, " (", ft$polarity, ", ", ft$creationDate, ")")
          )

          exp_nodes[[exps$name[j]]] <- structure(
            ft_nodes,
            attr  = list(a_attr = list(title = sprintf(
              "%s%s", exps$creationDate[j],
              if (!is.na(exps$description[j]) && nzchar(exps$description[j])) paste0(" – ", exps$description[j]) else ""
            ))
          ))


        }
        tree[[project_name]] <- exp_nodes
      }

      tree <- set_node_attrs(tree, attr_name = "stdisabled", inner_val = TRUE, leaf_val = FALSE )
      MetaboScapeProjects(tree)
      
    }
    show_metaboProjectSelection_modal()
  }
  
  show_metaboProjectSelection_modal <- function() {
    
    output$metaboProjectTree <- renderTree({
      MetaboScapeProjects()
    })
    
    showModal(modalDialog(
      
      title = "Project Table Selection",
      shinyTree("metaboProjectTree", search = TRUE, theme = "proton"),
      
      footer = tagList(
        modalButton("Cancel"),
        actionButton("SelectProject", "Select")
      ),
      easyClose = FALSE
    ))
  }
  
  observeEvent(input$SelectProject,({
    table <- get_selected(input$metaboProjectTree)
    MetaboSelectedProjectID(attr(table[[1]], "stId"))
    removeModal()
    state$newMS2data <- TRUE
    MS2file_name(table[[1]][1])
  }))

  observeEvent(input$setFlagsToMetaboScape,({
    tableID <- MetaboSelectedProjectID()
    if(tableID == "") {metaboAuth$require_auth(function() {
      selectProject_from_metaboscape()   
    })}
    
    showNotification("Not implemented yet")
    return()
    MetaboScapeAPI <- getMetaboScapeConnection()
    res <- MetaboScapeAPI$addFlagToFeature(featureId = "?????", flag = "BLUE")
    
  }))


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
    
    if (state$screen == "Setup") {
      shinyjs::show("welcome_screen")
    } else if (state$screen == "Analysis") {
      shinyjs::show("cards_section")
      shinyjs::show("analysis_footer")
      shinyjs::show("analysis_screen")
    } else if (state$screen == "Results") {
      shinyjs::show("results_screen")
    }
  }, ignoreNULL = FALSE)
  
  
  #-----------------------------------------------------------------------------
  # Screens Outputs
  #-----------------------------------------------------------------------------
  
  output$welcome_ui <- renderUI({ 
    tagList(
      fluidRow(
        column(6,
               div(class = "card required-card",
                   div(class = "required-badge", "Required"),
                   div(class = "card-title", "Required Inputs"),
                   div(
                     style = "margin-left: auto; min-width: 150px;",
                     fileInput("ms2_file", "MS2 Data", placeholder = ifelse(is.null(MS2file_info()),"No file selected", MS2file_info()),
                               accept = c(".mgf", ".msp", ".mzmine", ".mzTab")),
                     div(
                       style = "text-align: right;",
                       textOutput("MS2file_info", inline = TRUE)
                     ),
                     actionButton("extract_ms2", "MS2 from MetaboScape API")
                   ),
                   
                   selectInput("mode", "Ionization Mode",
                               choices = c("Positive" = "pos", "Negative" = "neg"),
                               selected = state$mode)
                   
               ),
               div(class = "flex-distance",
                   actionButton("start_analysis", "Start Analysis",
                                class = "process-btn"),
                 div(id = "load_project_wrapper",
                     fileInput("load_project_file", label = NULL, accept = ".zip",
                               buttonLabel = tagList(icon("folder-open"), "Load Project"),
                               placeholder = "")
                 )
               )
               
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
                             accept = c(".csv", ".tsv"))
                   

               )
        )
      )
    )     

    })
  
  
  output$saved_brkrrba <- renderText({
    brkrrba_save()
  })
  
  output$analysis_ui <- renderUI({ 
    
    tagList(
      # Enable Bootstrap tooltips
      tags$script(HTML("
            $(function () {
              $('[data-bs-toggle=\"tooltip\"]').tooltip();
            });
          ")),
      createAnalysisLegend(),
      fluidRow(
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
  ))
  })
  
  output$results_ui <- renderUI({ 
    tagList(
      verbatimTextOutput("log"),
      actionButton("setFlagsToMetaboScape", "Set Flags to MetaboScape")
    )
    })
  
  output$score_table <- renderTable({
    req(selected_feature_id())
    i <- as.numeric(selected_feature_id())
    res <- analysis_results$results
    df <- res[res$FEATURE_ID == i, 2:ncol(res)]
    
    df[] <- lapply(df, function(col) {
      if (is.logical(col)) {
        ifelse(col,
               "<span style='color:green;font-weight:bold;'>&#10004;</span>",
               "<span style='color:red;font-weight:bold;'>&#10008;</span>")
      } else {
        col
      }
    })
    colnames(df) <- gsub("[._]", " ", colnames(df))
    t(df)
  }, striped = TRUE, rownames = TRUE, colnames = FALSE, sanitize.text.function = function(x) x)
  
  # output$feature_info_table <- renderTable({
  #   req(nrow(analysis_results$results) > 0)
  #   req(selected_feature_id())
  #   res <- analysis_results$results
  #   featureID <- as.numeric(selected_feature_id())
  #   i <- which(res$FEATURE_ID == featureID)
  #   
  #   data <- spectraData(sp)[sp$FEATURE_ID == featureID]
  #   results <- results[,c("FEATURE_ID")]
  #   
  #   df <- res[res$FEATURE_ID == i, 2:ncol(res)]
  # 
  # 
  # }, striped = TRUE)
  
  
  output$featureTable <- DT::renderDataTable({
    req(analysis_results$results)
    res <- analysis_results$results

    spectrum_info <- as.data.frame(spectraData(analysis_results$spectra)[match(res$FEATURE_ID, spectraData(analysis_results$spectra)$FEATURE_ID), c("FEATURE_ID", "precursorMz", "rtime")])
    spectrum_info$score <- res$score
    
    DT::datatable(
      spectrum_info,
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
  
  output$ms1_plot <- renderPlot({
    req(selected_feature_id())
    req(ms1_data())
    
    plotMS1(ms1_data(), selected_feature_id())
  })
  
  # feature details
    output$spectrumDetails <- renderUI({
      if (is.null(selected_feature_id())) {
        return(div(class = "placeholder-text", "Select Feature"))
      }
      req(nrow(analysis_results$results) > 0)
      featureID <- as.numeric(selected_feature_id())

      spectrum_info <- analysis_results$spectra[analysis_results$spectra$FEATURE_ID == featureID]
      pol_col <- intersect(c("POLARITY", "polarity", "IONMODE", "ionmode"),spectraVariables(analysis_results$spectra))[1]
      adduct_col <- intersect(c("Adduct", "adducts", "ADDUCT", "ION"),spectraVariables(analysis_results$spectra))[1]
  
      tagList(
        #plotOutput("spectrum_plot"),
        plotlyOutput("interactive_spectrum"),
        div(class = "detail-row",
            div(class = "detail-label", "Feature ID"),
            div(class = "detail-value", as.character(spectrum_info$FEATURE_ID))
        ),
        tableOutput("score_table") ,
        div(class = "detail-row",
            div(class = "detail-label", "m/z"),
            div(class = "detail-value", as.character(format(spectrum_info$precursorMz, digits = 6)))
        ),
        div(class = "detail-row",
            div(class = "detail-label", "Retention time"),
            div(class = "detail-value", as.character(format(spectrum_info$rtime, digits = 4)))
        ),
        div(class = "detail-row",
            div(class = "detail-label", "Adduct"),
            div(class = "detail-value", paste(unlist(spectrum_info[[adduct_col]]), collapse = ", "))
        ),
        div(class = "detail-row",
            div(class = "detail-label", "Polarity"),
            div(class = "detail-value", spectrum_info[[pol_col]])
        ),
        plotOutput("ms1_plot")
      )
    })

    

  

  
}

