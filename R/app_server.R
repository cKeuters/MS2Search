#plan(multisession)
#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#' @import shiny
#' @import shinyFiles
#' @import shinyWidgets
#' @importFrom jsonlite fromJSON
#' @import zip
#' @importFrom future future plan
#' @import promises
#' @import readr
#' @import shinyTree
#' @importFrom stats setNames
#' @importFrom utils write.csv read.csv combn
#' @importFrom plotly renderPlotly plotlyOutput
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
    min_filtering_type = "fixed",
    fixed_threshold = 20,
    percent_treshold = 1,
    noiseFiltering_method = "none",
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
    res_settings = NULL,
    priorities = NULL,
    running = FALSE,
    log = NULL
  )
  
  
  log_state <- reactiveValues(file = NULL)
  log_text <- reactiveVal("")
  
  brkrrba_save <- reactiveVal("")

  ms1_data <- reactiveVal(NULL)
  
  #MetaboScape Authentification and Projects, Versions: "2025b" or "2026b"
  metaboAuth <- reactiveValues(
    authenticated = FALSE, 
    MetaboVersion = "2025b",
    pending_action = NULL
    )
  MetaboScapeProjects <- reactiveVal(list())
  MetaboSelectedProjectID <- reactiveVal("")
  MetaboAnnotationData <- reactiveVal(NULL)
  
  MS2file_name <- reactiveVal("")
  
  checked_features <- reactiveVal(NULL) 
  score_threshold <- reactiveVal(0.5)
  
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
  #Noise filtering
  observeEvent(input$fixedThreshold_value, {
    basicFilterValues$fixed_threshold <- input$fixedThreshold_value
  })
  observeEvent(input$percentThreshold_value, {
    basicFilterValues$percent_threshold <- input$percentThreshold_value
  })
  observeEvent(input$radio_minIntensity, {
    basicFilterValues$min_filtering_type <- input$radio_minIntensity
  })
  observeEvent(input$radio_FilteringMethod, {
    basicFilterValues$noiseFiltering_method <- input$radio_FilteringMethod
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

  observeEvent(input$ms1_file,{
    ms1_data(read_csv(input$ms1_file$datapath))
  })
  
  #-----------------------------------------------------------------------------
  # ObserveEvent of Outputs and Buttons
  #-----------------------------------------------------------------------------
  
  observeEvent(input$extract_ms2, {
    metaboAuth$require_auth(function() {
      selectProject_from_metaboscape()   
    })
  })
  
  
  # Output MS2 file name
  MS2file_info <- reactive({
    if (is.null(input$ms2_file)) {
      NULL
    } else {
      state$newMS2data <- TRUE
      MetaboSelectedProjectID("")
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
  
  # Button Start Analysis
  observeEvent(input$start_analysis, {
    if(MS2file_name() == ""){
      showNotification("No MS2 data selected", type="warning")
    }else{
      if(state$project_name == "") state$project_name <- paste0("project_",format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
      state$screen <- "Analysis"
    }

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
      res$selected <- checked_features()
      
      convert_list <- function(x) {
        if (is.list(x)) {
          vals <- sapply(x, function(y) {
            if (length(y) == 0) "" else paste(y, collapse = "/")
          })
        } else {
          x
        }
      }
      res <- lapply(res, convert_list)
      
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
        lipid_definition = lipid_definition())
      create_brkrrba_file(settings_list, file)
    }
  )
  
  #About Page
    mod_about_server("about_1")
    
    observeEvent(
      c(input$back_to_welcome, input$back_to_analysis, input$back_to_results),
      {
        updateTabsetPanel(session, "page", selected = "main")
      },
      ignoreInit = TRUE,
      priority = 10 
    )
    
    
    observeEvent(input$about_open, {
      updateTabsetPanel(session, "page", selected = "about")
      shinyjs::runjs("window.scrollTo(0, 0);")  
    })
    
    observeEvent(input$about_back, {
      updateTabsetPanel(session, "page", selected = "main")
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
  
  
  output$export_project <- downloadHandler(
    filename = function() {
      if(input$project_name == "") state$project_name <- paste0("project_",format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
      paste0(state$project_name, ".zip")
    },
    content = function(file) {
      save <- setup_and_save_project(create_settings_list(
        state = reactiveValuesToList(state),
        basicFilterValues = reactiveValuesToList(basicFilterValues),
        toolSettings = reactiveValuesToList(ToolSettings),
        cards   = all_cards_data(), 
        lipid_definition = lipid_definition()),
        analysis_results$res_settings,
        analysis_results$results,
        analysis_results$spectra,
        file)
    }
  )
  

  # Load Project
  observeEvent(input$load_project_file, {
    loaded_project <- load_project(input$load_project_file$datapath)
    settings_list <- loaded_project$settings_list
    analysis_results$res_settings <- loaded_project$res_settings
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
    basicFilterValues$min_filtering_type <- bf$min_filtering_type
    basicFilterValues$fixed_threshold <- bf$fixed_threshold
    basicFilterValues$percent_treshold <- bf$percent_treshold 
    basicFilterValues$noiseFiltering_method  <- bf$noiseFiltering_method 
    basicFilterValues$tolerance_ppm <- bf$tolerance_ppm
    basicFilterValues$minimum_dalton <- bf$minimum_dalton
    
    for(rule in settings_list$specific_filters){
      id <- next_card_id()
      rule$id <- id
      next_card_id(id + 1)
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
      metabolite_choices = METABOLITES$name
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
  
  lipidDefinitionServer <- mod_LipidDefinition_server("LipidDefinition_1", lipid_load_trigger)
  
  lipid_definition <- reactive({
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
    
    res <- filtered_results()
    res$FEATURE_ID[input$featureTable_rows_selected]
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
    

    settings <- create_settings_list(
      state             = reactiveValuesToList(state),
      basicFilterValues = reactiveValuesToList(basicFilterValues),
      toolSettings      = reactiveValuesToList(ToolSettings),
      cards             = all_cards_data(),
      lipid_definition   = lipid_definition()
      )
    
    analysis_results$res_settings <- settings

    prio <- c(settings$basicFilters$mz_priority, 
              settings$basicFilters$rt_priority, 
              settings$basicFilters$n_peaks_priority)
    prio <- c(prio, unlist(lapply(settings$specific_filters, function(x){x$rule_priority})))
    prio <- c(prio, unlist(lapply(settings$lipid_definition$chains, function(x){x$priority})))
    prio <- c(prio, unlist(lapply(settings$lipid_definition$validations, function(x){x$val_priority})))
    
    
    log_file <- tempfile(paste0("analysis_log_", state$project_name), fileext = ".txt")
    writeLines("", log_file)
    log_state$file <- log_file
    
    
    if(SelectedProjectID != "") MetaboScapeAPI <- getMetaboScapeConnection(metaboAuth$MetaboVersion)
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
          
          if(SelectedProjectID != ""){
            write_log(paste("Loading with MetaboSelectedProjectID", SelectedProjectID))
            featureTable <- createSpectradf(operations = MetaboScapeAPI, featureTableId = SelectedProjectID )
            spec <- Spectra(featureTable)
          } 
          else if(!is.null(ms2_path)) spec <- load_spectra(ms2_path, add_log = write_log)
          else write_log("NO MS2")
        }
        else spec <- analysis_results$spectra

        res <- run_analysis(settings, spec, add_log = write_log)
        
        list(
          spectra = spec,
          results = res,
          priorities = prio,
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
      analysis_results$priorities <- .$priorities
      analysis_results$log <- logs
      
    } %...!% {
      # error handler
      analysis_results$running <- FALSE
      
      message("Future failed:")

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
  mod_MetaboScape_auth_server("auth", metaboAuth)
  
  selectProject_from_metaboscape <- function(){
    if(length(MetaboScapeProjects()) == 0){
      
      MetaboScapeAPI <- getMetaboScapeConnection(metaboAuth$MetaboVersion)
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

          if(metaboAuth$MetaboVersion == "2025b") exps$description[j] <- ""
          exp_nodes[[exps$name[j]]] <- structure(
            ft_nodes,
            attr  = list(a_attr = list(title = sprintf(
              "%s%s", exps$creationDate[j],
              if (!is.na(exps$description[j]) && nzchar(exps$description[j])) paste0(" \u2013 ", exps$description[j]) else ""
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
      
      title = div(class = "modal-title-flex",
      "Project Table Selection",
                  actionButton("Metabo_logout", "Logout")),
      shinyTree("metaboProjectTree", search = TRUE, theme = "proton"),
      
      footer = tagList(
        modalButton("Cancel"),
        actionButton("SelectProject", "Select")
      ),
      easyClose = FALSE
    ))
  }
  
  observeEvent(input$Metabo_logout,{
    metaboAuth$authenticated <- FALSE
    metaboAuth$require_auth(function() {
      selectProject_from_metaboscape()   
    })
  })
  
  observeEvent(input$SelectProject,({
    table <- get_selected(input$metaboProjectTree)
    MetaboSelectedProjectID(attr(table[[1]], "stId"))
    removeModal()
    state$newMS2data <- TRUE
    getMetaboMS1data()
    MetaboAnnotationData(getMetaboPrimaryAnnotation())
    MS2file_name(table[[1]][1])
  }))
  
  getMetaboMS1data <- function(){
    MetaboScapeAPI <- getMetaboScapeConnection(metaboAuth$MetaboVersion)
    res <- MetaboScapeAPI$retrieveIntensityMatrix(MetaboSelectedProjectID())
    IntensityMatrix <- fromJSON(httr::content(res, as = "text", encoding = "UTF-8"))
    res <- MetaboScapeAPI$listAllSamples(MetaboSelectedProjectID())
    samples <- fromJSON(httr::content(res, as = "text", encoding = "UTF-8"))
    
    df <- as.data.frame(IntensityMatrix$intensities)
    colnames(df) <- lapply(samples$analysis, function(x) x$name)[[1]]
    df$FEATURE_ID <- seq(nrow(df))
    ms1_data(df)
  }
  
  # @title Get MetaboScape Primary Annotation  
  #
  # @return Returns the MetaboScape FeatureTable with MetaboScape Feature Id, primary annotation and annotation tool.
  getMetaboPrimaryAnnotation <- function(){
    MetaboScapeAPI <- getMetaboScapeConnection(metaboAuth$MetaboVersion)
    res <- MetaboScapeAPI$retrieveFeatureTable(featureTableId = MetaboSelectedProjectID())
    featureTable <- fromJSON(httr::content(res, as = "text", encoding = "UTF-8"))
    
    data.frame(
      id = featureTable$id,
      primaryAnnotation.name = featureTable$primaryAnnotation$name,
      primaryAnnotation.formula = featureTable$primaryAnnotation$formula,
      primaryAnnotation.tool = featureTable$primaryAnnotation$tool
      
    )

  }
   

  observeEvent(input$setFlagsToMetaboScape,({
    res <- analysis_results$results
    req(res)
    
    if(!"metaboScape_ID" %in% colnames(res)){
      showNotification("The Analysis had to be run on MetaboScape data from a
                       connected Server. Connect to a server in Setup, select 
                       a feature table and repeat the analysis.", 
                       duration = 10, type = "warning")
    }else{
      MetaboScapeAPI <- getMetaboScapeConnection(metaboAuth$MetaboVersion)
      current <- isolate(checked_features())
      
      true_ids <- names(current)[current]
      true_ids <- na.omit(true_ids)
      
      for(feature in true_ids){

        
        metaboID <- res[which(res$FEATURE_ID == feature),]$metaboScape_ID
        MetaboScapeAPI$addFlagToFeature(featureId = metaboID, flag = flag_color()$color)
        

        
      }
    }

      

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
  # Setup Screen
  #-----------------------------------------------------------------------------
  
  output$welcome_ui <- renderUI({ 
    tagList(
      fluidRow(
        column(6,
               div(class = "card required-card",
                   div(class = "required-badge", "Required"),
                   div(class = "card-title", "Required Inputs"),
                   div( class = "ms2Input-row",
                     div(
                       style = "min-width: 150px;",
                       fileInput("ms2_file", "MS2 Data", placeholder = ifelse(is.null(MS2file_info()),"No file selected", MS2file_info()),
                                 accept = c(".mgf", ".msp", ".mzmine", ".mzTab")),
                       
                     ),
                     actionButton("extract_ms2", "MS2 from MetaboScape API")
                   ),
                   div(
                     class = "input-badge",
                     textOutput("MS2file_info", inline = TRUE)
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
  
  #-----------------------------------------------------------------------------
  # Analysis Screen
  #-----------------------------------------------------------------------------
  
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
          width = 6,
          div(
            class = "card analysis-card",
            div(class = "card-title", "Basic filters"),
                div(
                  class = "basic-filter-container",
                  
                  div(
                    class = "filter-group",
                    h5("m/z"),
                    div(class = "input-block",
                        numericInput("mz_min", "Minimum", basicFilterValues$mz_min)),
                    div(class = "input-block",
                        numericInput("mz_max", "Maximum", basicFilterValues$mz_max)),
                    createPriorityRadioButtons("mz_priority",
                                               basicFilterValues$mz_priority,
                                               horizontal = TRUE)
                  ),
                  
                  div(
                    class = "filter-group",
                    h5("Retention Time"),
                    div(class = "input-block",
                        numericInput("rt_min", "Minimum (s)", basicFilterValues$rt_min)),
                    div(class = "input-block",
                        numericInput("rt_max", "Maximum (s)", basicFilterValues$rt_max)),
                    createPriorityRadioButtons("rt_priority",
                                               basicFilterValues$rt_priority,
                                               horizontal = TRUE)
                  ),
                  
                  div(
                    class = "filter-group",
                    h5("Number of peaks"),
                    div(class = "input-block",
                        numericInput("n_peaks_min", "Minimum", basicFilterValues$n_peaks_min)),
                    div(class = "input-block",
                        numericInput("n_peaks_max", "Maximum", basicFilterValues$n_peaks_max)),
                    createPriorityRadioButtons("n_peaks_priority",
                                               basicFilterValues$n_peaks_priority,
                                               horizontal = TRUE)
                  )
                )
              
            )
            

          ),
        column(width = 6,
               div(
                 class = "card analysis-card",
                 
                 div(class = "card-title", 
                     actionButton("button_MS2_intensities", NULL, icon = icon("square-poll-vertical")),
                     "MS2 Noise Filtering"),
                 
                 div(                  
                   class = "basic-filter-container",

                 div(
                   class = "filter-group-radio",
                   
                   div(
                     class = "filter-row",
                     div(class = "radio-group-inline",
                     radioButtons(
                       "radio_minIntensity",
                       label = "Minimum MS2 Intensity",
                       selected = basicFilterValues$min_filtering_type,
                       choiceNames = c("Fixed Threshold",
                                   "Percent from max Peak",
                                   "None"),
                       choiceValues = c("fixed", "percent", "none")
                     )),
                     
                     conditionalPanel(
                       "input.radio_minIntensity == 'fixed'",
                       numericInput(
                         "fixedThreshold_value",
                         "Minimum intensity threshold",
                         basicFilterValues$fixed_threshold
                       )
                     ),
                     
                     conditionalPanel(
                       "input.radio_minIntensity == 'percent'",
                       numericInput(
                         "percentThreshold_value",
                         "Percent from max Peak",
                         basicFilterValues$percent_treshold,
                         min = 0,
                         max = 100
                       )
                     )
                   )
                 ),
                 div(
                   class = "filter-group-radio",
                   div(class = "radio-group-inline",
                   radioButtons("radio_FilteringMethod", label = "MsCoreUtils Filtering Method",
                                choiceNames = c(" MAD (median absolute deviation)", "Friedman\u2019s SuperSmoother", "Weighted Moving Average", "None"),
                                choiceValues = c("MAD", "SuperSmoother", "wma", "none"),
                                selected = basicFilterValues$noiseFiltering_method)
                 )))
                 

               )
      
        )
      )
    )
    
    })
  
  observeEvent(input$button_MS2_intensities, {

    showModal(modalDialog(
      
      title = "Histogram of MS2 intensities",
      plotOutput(outputId = "distPlot")
      
      
      ,easyClose = TRUE
    ))
  })
  
  output$distPlot <- renderPlot({
    ms2_path <- input$ms2_file$datapath
    SelectedProjectID <- MetaboSelectedProjectID()
    
    if(state$newMS2data){
      if(!is.null(ms2_path)) {
        analysis_results$spectra <- load_spectra(ms2_path)
        state$newMS2data <- FALSE
        }
      else if(SelectedProjectID != ""){
        MetaboScapeAPI <- getMetaboScapeConnection(metaboAuth$MetaboVersion)
        
        featureTable <- createSpectradf(operations = MetaboScapeAPI, featureTableId = SelectedProjectID )
        analysis_results$spectra <- Spectra(featureTable)
        state$newMS2data <- FALSE
      } 
      else {
        showNotification("No MS2 data", type = "warning")
        return()
      }
    }
    hist_data <- data.frame(unlist(peaksData(analysis_results$spectra))[,"intensity"])
    colnames(hist_data) <- "intensity"
    
    ggplot(hist_data, aes(intensity)) +
      geom_histogram() +
      scale_x_log10()
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
  
  #-----------------------------------------------------------------------------
  # Result Screen - Top
  #-----------------------------------------------------------------------------
  
  flag_color  <- reactiveVal(list(color = "BLUE", rgb = "color: rgb(0,0,255)"))
  show_picker <- reactiveVal(FALSE)
  
  flag_choices <- list(list(color = "BLUE", rgb = "color: rgb(0,0,255)"),
                    list(color = "GREEN", rgb = "color: rgb(0,255,0)"),
                    list(color = "RED", rgb = "color: rgb(255,0,0)"))
  
  output$results_ui <- renderUI({
    tagList(
      verbatimTextOutput("log"),
      
      div(
        class = "action-toolbar",
        
        div(
          class = "toolbar-left",
          
          div(
            class = "threshold-group",
            
            numericInput(
              "score_threshold",
              NULL,
              value = score_threshold(),
              min = 0,
              max = 1,
              step = 0.1
            ),
            
            actionButton(
              "set_score_threshold",
              "Check",
              class = "toolbar-btn"
            )
          ),
          
          actionButton(
            "uncheckAll",
            "Uncheck all",
            class = "toolbar-btn"
          )
        ),
        
        div(
          class = "toolbar-right",
          
          div(class = "flag-wrapper",
              actionLink("changeFlagColor",
                         label = tags$i(class = "fa-solid fa-flag", style = flag_color()$rgb)),
              uiOutput("flag_picker")
          ),
          
          actionButton(
            "setFlagsToMetaboScape",
            "Set Flags to MetaboScape",
            class = "toolbar-primary"
          )
        )
      )
    )
  })
  
  # Toggle popup when the flag is clicked
  observeEvent(input$changeFlagColor, {
    show_picker(!show_picker())
  })
  
  # Popup content, rendered right below the flag
  output$flag_picker <- renderUI({
    req(show_picker())
    div(class = "flag-picker",
        lapply(flag_choices, function(color) {
          tags$a(
            href = "#", class = "flag-option",
            onclick = sprintf(
              "Shiny.setInputValue('flag_color_pick', '%s', {priority: 'event'}); return false;",
              color$color),
            tags$i(class = "fa-solid fa-flag", style = color$rgb)
          )
        })
    )
  })
  
  # Pick a color: update the flag and close the popup
  observeEvent(input$flag_color_pick, {
    match_idx <- which(sapply(flag_choices, function(x) x$color) == input$flag_color_pick)
    
    flag_color(flag_choices[[match_idx]])
    show_picker(FALSE)
  })
  
  # Close when clicking outside
  observeEvent(input$flag_picker_close, show_picker(FALSE))
  
  
  
  
  #-----------------------------------------------------------------------------
  # Result Screen - Feature Table
  #-----------------------------------------------------------------------------
  
  
  # --- helper: build a checkbox <input> for one row ---
  checkbox_html <- function(id, checked) {
    sprintf('<input type="checkbox" class="feature-checkbox" data-feature-id="%s"%s/>',
            id, ifelse(!is.na(checked) & checked, " checked", ""))
  }
  
  
  # JS that listens for clicks on any checkbox in the table and reports
  # back to Shiny which feature was (un)checked. Delegated on the table
  # body, so it survives DT's server-side redraws.
  checkbox_js <- DT::JS(
    "table.on('change', 'input.feature-checkbox', function() {",
    "  var id = $(this).data('feature-id');",
    "  var checked = $(this).is(':checked');",
    "  Shiny.setInputValue('feature_checkbox_change',",
    "    {id: id, checked: checked}, {priority: 'event'});",
    "});"
  )
  

  observeEvent(analysis_results$results, {
    checked_features(setNames(filtered_results()$score > score_threshold(), as.character(filtered_results()$FEATURE_ID)))
  }, ignoreNULL = TRUE)
  
  observeEvent(input$score_threshold, {
    score_threshold(input$score_threshold)
  })
  
  # Table data without the checkbox column, shared by render + proxy update
  table_data <- reactive({
    res <- filtered_results()

    sd <- spectraData(analysis_results$spectra,
                      columns = c("FEATURE_ID", "FEATURE_KEY", "precursorMz", "rtime"))
    info <- as.data.frame(sd[match(res$FEATURE_ID, sd$FEATURE_KEY),
                             c("FEATURE_KEY", "precursorMz", "rtime")])
    info$score <- res$score
    info
  })

  add_checkboxes <- function(info, checked) {
    ids <- as.character(info$FEATURE_KEY)
    info$selected <- checkbox_html(ids, unname(checked[ids]))
    info
  }
  
  output$featureTable <- DT::renderDataTable({
    req(analysis_results$results)
    spectrum_info <- add_checkboxes(table_data(), isolate(checked_features()))
    
    DT::datatable(
      spectrum_info,
      selection = "single",
      class = "feature-table",
      rownames = FALSE,
      escape = FALSE,
      colnames = c("Feature ID", "m/z", "Retention time", "Score", "Include"),
      options = list(
        pageLength = 25,
        lengthMenu = c(25, 50, 100),
        order = list(list(3, "desc")),
        columnDefs = list(list(orderable = FALSE, targets = 4))
      ),
      callback = checkbox_js
    ) |>
      DT::formatRound("precursorMz", 6) |>
      DT::formatRound("rtime", 4) |>
      DT::formatRound("score", 3)
  }, server = TRUE)
  
  proxy <- DT::dataTableProxy("featureTable")
  
  # Checkbox toggles: only update the stored state, table is untouched
  observeEvent(input$feature_checkbox_change, {
    info <- input$feature_checkbox_change
    current <- checked_features()
    current[[as.character(info$id)]] <- isTRUE(info$checked)
    checked_features(current)
  }, ignoreInit = TRUE)
  
  # Threshold button: update the state AND refresh the table contents
  observeEvent(input$set_score_threshold, {
    res <- filtered_results()
    new_state <- setNames(res$score > score_threshold(), as.character(res$FEATURE_ID))
    checked_features(new_state)
    
    DT::replaceData(
      proxy,
      add_checkboxes(table_data(), new_state),
      resetPaging = FALSE,
      clearSelection = "none",
      rownames = FALSE
    )
  })
  
  observeEvent(input$uncheckAll, {
    res <- filtered_results()
    new_state <- setNames(rep(FALSE, length(res$FEATURE_ID)), as.character(res$FEATURE_ID))
    checked_features(new_state)
    
    DT::replaceData(
      proxy,
      add_checkboxes(table_data(), new_state),
      resetPaging = FALSE,
      clearSelection = "none",
      rownames = FALSE
    )
  })

  
  #filter
  filter_cols <- reactive({
    req(analysis_results$results)
    res <- analysis_results$results

    start_filtercols <- ifelse("metaboScape_ID" %in% colnames(res), 3, 2)
    matches <- res[, start_filtercols:(ncol(res)-1)]
    
    is_match <- function(x) {
      if (is.list(x)) {
        vapply(x, function(el) {
          length(el) > 0 && !all(is.na(el))
        }, logical(1))
      } else {
        !is.na(x) & x
      }
    }
    
    lapply(matches, is_match)
  })
  
  # set checkboxes
  observeEvent(filter_cols(), {
    updateCheckboxGroupInput(session, "filter_select",
                             choices = names(filter_cols()))
  })
  
  #results, but just the ones passing the filters (checkboxes)
  filtered_results <- reactive({
    req(analysis_results$results)
    res <- analysis_results$results
    res_passed <- as.data.frame(filter_cols())
    colnames(res_passed) <- gsub("\\.", "_", colnames(res_passed))
    sel <- gsub(" ", "_",input$filter_select)
    
    if (length(sel) == 0) return(res)
    
    passed_cols <- res_passed[,sel, drop = FALSE]
    keep <- rowSums(passed_cols) == ncol(passed_cols)
    
    res[keep, , drop = FALSE]
  })
  
  
  
  #-----------------------------------------------------------------------------
  # Result Screen - Feature Details
  #-----------------------------------------------------------------------------
  
  set_isotope_sum <- reactiveVal(FALSE)
  
  selected_spectrum <- reactive({
    req(selected_feature_id())
    get_spectrum_by_id(analysis_results$spectra, selected_feature_id())
  })
  
  plot_ctx <- reactive({
    req(analysis_results$res_settings)
    prepare_plot_context(analysis_results$res_settings)
  })
  
  
  observeEvent(input$set_isotope_sum_Button,{
    set_isotope_sum(!set_isotope_sum())
  })
  
  output$filter_results_tableDT <- DT::renderDataTable({
    req(selected_feature_id())
    featureID  <- selected_feature_id()
    res <- analysis_results$results
    
    
    spectrum <- selected_spectrum()
    settings_list <- analysis_results$res_settings
    ctx <- plot_ctx()

    if(set_isotope_sum()) {
      res <- summarize_isotope_peaks(res, spectrum)}
    
    nl <- neutral_losses_from_results(res, featureID, ctx$rules_prepped, spectrum,
                                      tol = settings_list$basicFilters$minimum_dalton)
    
    i_pos <- which(res$FEATURE_ID == featureID)
    res <- res[i_pos, 2:ncol(res), drop = FALSE]
    colnames(res) <- gsub("[\\._]", " ", colnames(res))
    
    tick  <- "<span style='color:green;font-weight:bold;'>&#10004;</span>"
    cross <- "<span style='color:red;font-weight:bold;'>&#10008;</span>"

    df <- as.data.frame(filter_cols())[i_pos,]
    colnames(df) <- gsub("[\\._]", " ", colnames(df))
    
    # Row 1: match yes/no
    row1 <- ifelse(as.logical(df[1, ]), tick, cross)
    df[1,] <- row1
    
    precursorMz <- analysis_results$spectra[analysis_results$spectra$FEATURE_ID == featureID]$precursorMz

    #row 2 found values
    row2 <- vapply(res[,colnames(df)], function(x) {

      if (is.list(x)){
        value_text <- paste(round(as.numeric(unlist(x)), 3), collapse = ", ")
        diff <- lapply(x, function(x2){as.numeric(precursorMz) - as.numeric(x2)})
        diff_text <- paste(round(as.numeric(unlist(diff)), 3), collapse = ", ")

      }else
      {
        value_text <- x
        diff_text <- as.numeric(precursorMz) - as.numeric(x)
      }
      if(!is.logical(x) && !is.na(x)){
        sprintf(
          '<span title="Found Values: %s&#10;Difference to precursor: %s">%s</span>',
          value_text,
          diff_text,
          value_text
        )
      }else ""

      },
      character(1)
    )

    df <- rbind(df, row2) 

    # Row 3: formulas, only for the chain columns
    row3 <- setNames(rep("", ncol(df)), names(df))
    chain_calcs <- ctx$chain_calcs
    chain_cols  <- grep("chain", names(res))
    
    for (k in seq_len(min(length(chain_cols), length(chain_calcs)))) {
      cn       <- names(res)[chain_cols[k]]
      found    <- as.numeric(unlist(res[[cn]]))
      found <- found[!is.na(found)]
      all_data <- chain_calcs[[k]]
      if(length(found > 0)){
        tol <- max(max(found) * basicFilterValues$tolerance_ppm / 1000000, basicFilterValues$minimum_dalton)
        near <- rowSums(abs(outer(all_data$mass, found, "-")) <= tol) > 0
        
        formulas <- all_data$formula[near]
        
        row3[cn] <- paste(unique(formulas), collapse = ", ")
      }
      else row3[cn] <- ""
      
    }
    # Row 3: Sequential neutral loss combinations
    seq_nl_col <- lapply(settings_list$specific_filters, function(filter){
      if(filter$sequential) filter$label else NULL
    })
    
    seq_nl_col <- Filter(Negate(is.null), seq_nl_col)
    seq_nl_col <- unlist(seq_nl_col)
    
    for (col_name in seq_nl_col) {
      filtered_nl <- nl[nl$label == col_name,]
      if (is.null(filtered_nl) || NROW(filtered_nl) == 0) next
      
      annotation_txt <- list()
      for(nl_row in seq_len(nrow(filtered_nl))){
        mass_counts <- table(unlist(filtered_nl$masses[nl_row]))
        txt <- if (sum(mass_counts) > 1) {
          paste0(round(as.numeric(filtered_nl$mz_found[nl_row]), 2),": ",paste(paste0(mass_counts, "x ", names(mass_counts)), collapse = ", "))
        } else {
          names(mass_counts)
        }
        annotation_txt <- c(annotation_txt, txt)
      }

      row3[col_name] <- paste(annotation_txt, collapse = "; ")
      
    }
    df <- rbind(df, row3)    
    df <- t(df)
    
    
    DT::datatable(
      df,
      rownames = TRUE,
      escape = FALSE,
      selection = 'none',
      colnames = c("Filter", "passed", "found values", "formula"),
      options = list(
        paging = FALSE,
        rowCallback = DT::JS(
          rowCallback = DT::JS(
            sprintf(
              "
          function(row, data, index) {
          
            var p = %s[index];
          
            if (p === 'basic') {
              $(row).css(
                'background',
                'linear-gradient(to right, #28a745 0%%, white 20%%, white 100%%)'
              );
              $('td', row).css('background', 'transparent');
            }
          
            else if (p === 'important') {
              $(row).css(
                'background',
                'linear-gradient(to right, #fd7e14 0%%, white 20%%, white 100%%)'
              );
              $('td', row).css('background', 'transparent');
            }
          
            else if (p === 'required') {
              $(row).css(
                'background',
                'linear-gradient(to right,  #dc3545 0%%, white 20%%, white 100%%)'
              );
              $('td', row).css('background', 'transparent');
            }
          
          }
          ",
    jsonlite::toJSON(analysis_results$priorities)
            )
          )
        )
      )
    )
    
 
  }, server = TRUE)
  
  
  output$interactive_spectrum <- renderPlotly({
    req(nrow(analysis_results$results) > 0)
    req(selected_feature_id())
    settings_list <- analysis_results$res_settings
    results <- analysis_results$results
    results <- results[results$FEATURE_ID == selected_feature_id(),]
    plottingSpectra(selected_spectrum(), set_isotope_sum(), settings_list, results)
    
  })
  
  mod_spectrum_server(
    "spectrum_plot_mod",
    spectrum_r        = selected_spectrum,
    settings_list_r   = reactive(analysis_results$res_settings),
    results_r         = reactive({
      req(nrow(analysis_results$results) > 0, selected_feature_id())
      analysis_results$results[analysis_results$results$FEATURE_ID == selected_feature_id(), ]
    }),
    set_isotope_sum_r = set_isotope_sum,
    ctx_r             = plot_ctx
  )
  
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
      featureID <- selected_feature_id()

      spectrum_info <- analysis_results$spectra[analysis_results$spectra$FEATURE_KEY == featureID]
      pol_col <- intersect(c("POLARITY", "IONMODE", "ionmode", "polarity"),spectraVariables(analysis_results$spectra))[1]
      adduct_col <- intersect(c("Adduct", "adducts", "ADDUCT", "ION"),spectraVariables(analysis_results$spectra))[1]
      
      if(!is.null(MetaboAnnotationData())){
        annotationDetails <- MetaboAnnotationData()
        annotationDetails <- annotationDetails[featureID,]
      }else{
        calculated_formulas <- calculateFormulaFromMass(spectrum_info$precursorMz,unlist(spectrum_info[[adduct_col]])[1])
        
        annotationDetails <- data.frame(id = featureID, 
                                        primaryAnnotation.formula = calculated_formulas$formula[1], 
                                        primaryAnnotation.name = "-", 
                                        primaryAnnotation.tool = "Decomposed Mass")
      }

      tagList(
        div(
          class = "flex-distance",
          div(class = "info-pill", strong("Formula: "),
              annotationDetails$primaryAnnotation.formula
          ),
          div(class = "info-pill",strong("Name: "),
              annotationDetails$primaryAnnotation.name
          ),
          div(class = "info-pill",strong("Tool: "),
              annotationDetails$primaryAnnotation.tool
          )
        ),
        div(class = "detail-row",
          div(
              div(class = "detail-label", "Feature ID"),
              div(class = "detail-value", as.character(spectrum_info$FEATURE_ID))
          ),
          div(
              div(class = "detail-label", "m/z"),
              div(class = "detail-value", as.character(format(spectrum_info$precursorMz, digits = 6)))
          ),
          div(
              div(class = "detail-label", "Retention time"),
              div(class = "detail-value", as.character(format(spectrum_info$rtime, digits = 4)))
          ),
          div(
              div(class = "detail-label", "Adduct"),
              div(class = "detail-value", paste(unlist(spectrum_info[[adduct_col]]), collapse = ", "))
          ),
          div(
              div(class = "detail-label", "Polarity"),
              div(class = "detail-value", spectrum_info[[pol_col]])
          ) 
        ),
        #plotlyOutput("interactive_spectrum"),
        mod_spectrum_ui("spectrum_plot_mod"),
        DT::dataTableOutput("filter_results_tableDT"),
        actionButton("set_isotope_sum_Button", "Summarize Peaks"),
        plotOutput("ms1_plot")
      )
    })

}

mod_spectrum_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # shinyWidgets::checkboxGroupButtons() gives the same input as real buttons
    checkboxGroupInput(ns("overlays"), NULL, inline = TRUE,
                       choices = c("Neutral losses" = "nl", "Chains" = "chain", "Validations" = "val")),
    plotly::plotlyOutput(ns("spectrum_plot"))
  )
}

mod_spectrum_server <- function(id, spectrum_r, settings_list_r, results_r,
                                set_isotope_sum_r, ctx_r) {
  moduleServer(id, function(input, output, session) {
    ctx <- ctx_r      
    
    built <- reactive({
      req(spectrum_r(), nrow(results_r()) > 0)
      plottingSpectra(spectrum_r(), set_isotope_sum_r(), settings_list_r(), results_r(),
                      ctx(), show = isolate(input$overlays))
    })
    
    output$spectrum_plot <- plotly::renderPlotly(built()$fig)
    
    proxy <- plotly::plotlyProxy("spectrum_plot", session)
    
    observeEvent(input$overlays, ignoreNULL = FALSE, ignoreInit = TRUE, {
      idx <- built()$idx
      for (grp in names(idx)) {
        if (length(idx[[grp]]) == 0) next
        plotly::plotlyProxyInvoke(proxy, "restyle",
                                  list(visible = grp %in% input$overlays),
                                  as.list(idx[[grp]]))
      }
    })
  })
}