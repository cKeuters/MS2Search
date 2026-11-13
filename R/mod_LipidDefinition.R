#' @title Lipid definition UI Function
#'
#' @description The Shiny UI module for the chain and chain validation enviroment.
#' 
#' @details
#' A panel with the headgroup and with one container for the chains and one 
#' for the validations.
#'
#' @param id Shiny module id
#'
#' @importFrom shiny NS tagList textInput actionButton
#' @noRd 
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
        div(
          class = "headgroup-box",
          textInput(ns("shorthand"), label = "Shorthand", placeholder = "z.B. CER", width = "160px")
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



#' @title Lipid definition Server Function
#'
#' @description The Shiny Server module the lipid definition coordination.
#' 
#' @details
#' Updates the headgroup on input. Adds chains to the UI, creates their Servers 
#' and stores them in a reactiveValue. Does the same with the validations.
#' 
#' @param id Shiny module id
#' @param load_trigger Holds the chains when new data is loaded.
#' @param input,output,session Internal parameters for {shiny}.
#' 
#' @importFrom shinyvalidate InputValidator sv_regex sv_optional
#' 
#' @return Returns a list with a chain list and a validation list.
#' 
#' @noRd 
mod_LipidDefinition_server <- function(id, load_trigger){
  moduleServer(id, function(input, output, session){
    ns <- session$ns
    
    main_iv <- shinyvalidate::InputValidator$new()
    main_iv$add_rule("headgroup", shinyvalidate::sv_regex("^$|^[=]?([A-Z][a-z]?\\d*)([=]?[A-Z][a-z]?\\d*)*$",
                                                          "Invalid formula"))
    main_iv$enable()

    head_group <- reactiveVal("")
    
    observeEvent(input$headgroup, {
      head_group(input$headgroup)
    })

    shorthand <- reactiveVal("")
    
    observeEvent(input$shorthand, {
      shorthand(input$shorthand)
    })
    
    #-----------------------------------------------------------------------------
    # Chains
    #-----------------------------------------------------------------------------
    
    chain_modules <- reactiveValues()
    chain_ids <- reactiveVal(character(0))

    remove_chain <- function(cid) {
      cid <- as.character(cid)
      removeUI(selector = paste0("#", ns(paste0("chain_wrap_", cid))))
      chain_modules[[cid]] <- NULL

      chain_ids(setdiff(chain_ids(), cid))
    }
    
    observeEvent(input$add_chain, {
      chain <- empty_chain
      chain$chain_id <- create_unique_ids()
        
      add_chain_to_ui(init_chain = chain)
    }, ignoreInit = TRUE)
    
    add_chain_to_ui <- function(init_chain){
      
      cid <- init_chain$chain_id
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
      main_iv$add_validator(chain_modules[[cid]]$validator)
      
    }

    
    #-----------------------------------------------------------------------------
    # Validations
    #-----------------------------------------------------------------------------
    
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
        ui = div(id = ns(paste0("validation_wrap_", vid)), 
                 validationUI(ns(child_id), init_validation, fixed)),
        immediate = TRUE
      )
      validation_modules[[vid]] <<- validationServer(
        child_id, init_validation, chain_ids, fixed, on_remove = remove_validation(vid))
      main_iv$add_validator(validation_modules[[vid]]$validator)
      
      next_validation_id(as.numeric(vid) + 1)
    }

    
    # always-present precursor validation
    session$onFlushed(function() {
      withReactiveDomain(session, {
        isolate(add_validation_to_ui(
          list(name = "precursor validation", adjustments = "",
               val_priority = "basic", components = list()),
          fixed = TRUE
        ))
      })
    }, once = TRUE)
    
    observeEvent(input$add_validation, {
      add_validation_to_ui(empty_validation)
    })
    
    
    #---- if data was loaded, adds the chains and validations to the UI -----
    observeEvent(load_trigger(), {
      req(load_trigger())
      
      for(id in chain_ids()){
        remove_chain(id)
      }
      
      
      chains_to_load <- load_trigger()$chains
      
      if(!is.null(chains_to_load)){
        isolate({
          for (chain in chains_to_load) {
            tryCatch({
              add_chain_to_ui(chain)
            }, error = function(e) {
              cat("ERROR in add_chain_to_ui:", conditionMessage(e), "\n")
            })
          }
        })
      }
      vals_to_load <- load_trigger()$validations

      if(!is.null(vals_to_load)){
        for(val_id in validation_ids()) remove_validation(val_id)
        
        isolate({
          for (val in vals_to_load) {
            tryCatch({
              add_validation_to_ui(val)
            }, error = function(e) {
              cat("ERROR in add_validation_to_ui:", conditionMessage(e), "\n")
            })
          }
        })
      }
      
      head_group(load_trigger()$head_group)
      updateTextInput(session = session,"headgroup",value = head_group())
      shorthand(load_trigger()$shorthand)
      updateTextInput(session = session,"shorthand",value = shorthand())
    }, ignoreNULL = TRUE, ignoreInit = TRUE)
    
    
    #Return of the function
    reactive({list(
      chains = ({
        ids <- chain_ids()
        if (length(ids) == 0) return(list())
        lapply(ids, function(cid) chain_modules[[cid]]$chain_data())
      }),
      validations = ({
        if (length(validation_ids()) == 0) return(list())
        lapply(validation_ids(), function(v){
          validation_modules[[v]]$val_data()
        })
      }),
      head_group = head_group(),
      shorthand = shorthand(),
      is_valid = main_iv$is_valid()
      
    )
    })
    
  })
}


#' @title Chain Block UI Function
#'
#' @description The Shiny UI module for one chain.
#' 
#' @details
#' UI for setting the parameters of one chain. Controls modifications, the svg 
#' of the chain, and the mz information of possible chain combinations.
#'
#' @param id Shiny module id
#' @param init_chain Initial parameters
#'
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
        tags$label("Even Steps"),
        tags$label(
          class = "switch",
          tags$input(
            id = ns("even_steps"),
            type = "checkbox",
            checked = init_chain$even_steps
          ),
          tags$span(class = "slider round")
        )
      ),
      div(
        class = "control-group",
        tags$label("C-Atoms"),
        div(
          class = "control-inputs",
          numericInput(ns("c_min"), "min", value = init_chain$C_min,
                       min = 0, max = 100,
                       step = if (isTRUE(init_chain$even_steps)) 2 else 1,
                       width = "62px"),
          numericInput(ns("c_max"), "max", value = init_chain$C_max,
                       min = 0, max = 100,
                       step = if (isTRUE(init_chain$even_steps)) 2 else 1,
                       width = "62px")
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
    )  )
}

#' @title Chain Block Server Function
#'
#' @description The Shiny Server module for one chain.
#' 
#' @details
#' It is updating the parameters from the UI. It coordinates the modifications 
#' modules.
#' 
#' @param server_id Shiny module id
#' @param id ID of the chain, used for removing the chain.
#' @param init_chain Parameters to initialize the chain.
#' @param on_remove Function to delete the chain.
#' @param input,output,session Internal parameters for {shiny}.
#' 
#' @return Returns a list with C_min, C_max, DB_min, DB_max, priority, modifications
#' 
#' @noRd 
chain_block_server <- function(id, module_id, init_chain, on_remove ) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    chain_iv <- shinyvalidate::InputValidator$new()
    
    # ReactiveValues of the Chain parameters
    #
    # @details
    # Variables and values:
    # \describe{
    #   \item{C_min}{Minimum number of carbon atoms for the chain. (initial 10)}
    #   \item{C_max}{Maximum number of carbon atoms for the chain. (initial 20)}
    #   \item{even_steps}{Determines if the chain only has C-values divided by 2. (initial TRUE)}
    #   \item{DB_min}{Minimum number of double bonds in the chain (initial 0)}
    #   \item{DB_max}{Maximum number of double bonds in the chain (initial 1)}
    #   \item{priority}{Which scoring priority: basic, important, required (initial basic)}
    #   \item{modifications}{List of modications added to the chain.}
    # }
    chain_state <- reactiveValues(
      chain_id = init_chain$chain_id,
      C_min = init_chain$C_min,
      C_max = init_chain$C_max,
      even_steps = init_chain$even_steps,
      DB_min = init_chain$DB_min,
      DB_max = init_chain$DB_max,
      priority = init_chain$priority,
      modifications = NULL
    )
    mod_ids <- reactiveVal(1)

    
    chain_mz_values <- reactive({
      req(chain_iv$is_valid())
      calc_lipid_chain_masses(reactiveValuesToList(chain_state))
    })
    
    output$chain_priority <- renderUI({
      createPriorityRadioButtons(ns("chain_priority"), chain_state$priority)
    })
    observeEvent(input$chain_priority,{
      chain_state$priority <- input$chain_priority
    })
    
    observeEvent(input$remove_chain, on_remove(module_id), ignoreInit = TRUE)
    
    #--------------------------------------------------------------------------
    # Validating NumericInputs
    #--------------------------------------------------------------------------
    
    # helper: only send an update if the value really differs
    update_if_changed <- function(id, current, new, ...) {
      if (is.null(current) || is.na(current) || !isTRUE(all.equal(current, new))) {
        updateNumericInput(session, id, value = new, ...)
      }
    }
    
    observeEvent(input$even_steps, {
      chain_state$even_steps <- isTRUE(input$even_steps)
      step_size <- if (isTRUE(input$even_steps)) 2 else 1
      updateNumericInput(session, "c_min", step = step_size)
      updateNumericInput(session, "c_max", step = step_size)
    }, ignoreInit = FALSE)
    
    # --- C atoms ---------------------------------------------------------------
    observeEvent(input$c_min, {
      req(input$c_min)                       # ignore NA while the user is typing
      v <- min(max(input$c_min, 0), 100)
      if (v != input$c_min) updateNumericInput(session, "c_min", value = v)
      
      chain_state$C_min <- v
      updateNumericInput(session, "c_max", min = v)          # only the limit
      if (!is.null(input$c_max) && !is.na(input$c_max) && input$c_max < v) {
        updateNumericInput(session, "c_max", value = v)      # only if violated
      }
    }, ignoreInit = TRUE)
    
    observeEvent(input$c_max, {
      req(input$c_max)
      v <- min(max(input$c_max, chain_state$C_min), 100)
      if (v != input$c_max) updateNumericInput(session, "c_max", value = v)
      
      chain_state$C_max <- v
      test_maxDB()
    }, ignoreInit = TRUE)
    
    # --- Double bonds ----------------------------------------------------------
    observeEvent(input$db_min, {
      req(input$db_min)
      max_db <- floor(chain_state$C_max / 2)
      v <- min(max(input$db_min, 0), max_db)
      if (v != input$db_min) updateNumericInput(session, "db_min", value = v)
      
      chain_state$DB_min <- v
      updateNumericInput(session, "db_max", min = v)
      if (!is.null(input$db_max) && !is.na(input$db_max) && input$db_max < v) {
        updateNumericInput(session, "db_max", value = v)
      }
    }, ignoreInit = TRUE)
    
    observeEvent(input$db_max, {
      req(input$db_max)
      print(input$db_max)
      max_db <- floor(chain_state$C_max / 2)
      v <- min(max(input$db_max, chain_state$DB_min), max_db)
      if (v != input$db_max) updateNumericInput(session, "db_max", value = v)
      
      chain_state$DB_max <- v
    }, ignoreInit = TRUE)
    
    # --- Only adjust limits, and only change values that violate them ----------
    test_maxDB <- function() {
      max_db <- floor(chain_state$C_max / 2)
      updateNumericInput(session, "db_max", max = max_db)
      updateNumericInput(session, "db_min", max = max_db)
      
      if (!is.null(input$db_max) && !is.na(input$db_max) && input$db_max > max_db) {
        updateNumericInput(session, "db_max", value = max_db)
      }
      if (!is.null(input$db_min) && !is.na(input$db_min) && input$db_min > max_db) {
        updateNumericInput(session, "db_min", value = max_db)
      }
    }
    #--------------------------------------------------------------------------
    # Modification coordination
    #--------------------------------------------------------------------------
    
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
      }, ignoreInit = TRUE, once = TRUE)
    }
    
    register_mod_validation <- function(mid) {
      field <- paste0(mid, "_formula")
      
      chain_iv$add_rule(field, function(value) {
        if (is.null(chain_state$modifications[[mid]])) return(NULL)
        if (length(value) > 0 && !grepl("^(=?[A-Z][a-z]?|\\d)*$", value)) {
          return("Invalid formula")
        }        
        NULL
      })
    }
    
    add_mod_to_ui <- function(init_mod) {
      n_mods <- mod_ids()
      mid <- paste0("mods_", n_mods)
      target <- paste0("mods_above")
      
      
      insertUI(
        selector = paste0("#", ns(target)),
        where = "beforeEnd",
        ui = mod_chip_ui(ns, mid, init_mod)
      )

      register_mod_observers(mid, init_mod)
      register_mod_validation(mid)
      mod_ids(mod_ids() + 1)
    }
    
    isolate({
      init_mods <- init_chain$modifications
      for (init_mod in init_mods) {
        add_mod_to_ui(init_mod)
        }
    })
    
    #--------------------------------------------------------------------------
    # Chain information visualization
    #--------------------------------------------------------------------------
    
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
    
    # box with min and max mz values of chain
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
    
    chain_data <- reactive({
      list(
        chain_id = chain_state$chain_id,
        C_min = chain_state$C_min,
        C_max = chain_state$C_max,
        even_steps = chain_state$even_steps,
        DB_min = chain_state$DB_min,
        DB_max = chain_state$DB_max,
        priority = chain_state$priority,
        modifications = chain_state$modifications
      )
    })
    
    list(
      chain_data = chain_data,
      validator = chain_iv
    )
  })
}


#' @title Modification UI Function
#'
#' @description The Shiny UI module for one modification of one chain.
#' 
#' @details
#' UI for the input of one modification and if it is optional to the chain.
#'
#' @param ns Function for building Shiny module id
#' @param mid ID of the modification used for naming.
#' @param init_mod Parameters to initialize modification.
#'
#' @noRd 
mod_chip_ui <- function(ns, mid, init_mod) {
  div(
    id = ns(paste0("mod_chip_", mid)),
    class = "mod-chip",
    div(
      class = "formula-container",
      textInput(
        ns(paste0(mid, "_formula")),
        NULL,
        placeholder = "e.g. OH",
        value = init_mod$formula
      )
    ),
    checkboxInput(ns(paste0(mid,"_mod_optional")), "Optional", value = init_mod$optional),
    actionButton(ns(paste0(mid, "_remove")), "\u00d7", class = "btn-remove-mod")
  )
}

#' @title Modification Server Function
#'
#' @description The Shiny Server module for one modification of one chain.
#' 
#' @details
#' Its just there to give back the input values from the UI mod.
#' 
#' @param id Shiny module id
#' @param input,output,session Internal parameters for {shiny}.
#' 
#' @return Returns a list with the formula and the optionality.
#' 
#' @noRd 
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

#' @title Chain-Information Pop-Up UI Function
#'
#' @description The Shiny UI module for the Chain-Information.
#' 
#' @details
#' Its just a button with a i in it. The Pop-Up is done by the server when clicking
#' the button.
#'
#' @param id Shiny module id
#'
#' @noRd 
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

#' @title Chain-Information Pop-Up Server Function
#'
#' @description The Shiny Server module for the Chain-Information.
#' 
#' @details
#' It shows a dialog when the button in the UI is pressed and renders the
#' datatable in it.
#' 
#' @param id Shiny module id
#' @param df dataframe of the chain information containing the formula, optional
#' modifications and the masses.
#' @param title Title of the dialog (detault = NULL)
#' @param input,output,session Internal parameters for {shiny}.
#' 
#' @noRd 
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

#' @title Validation UI Function
#'
#' @description The Shiny UI module for one validation.
#' 
#' @details
#' The validation UI for a fixed validation adds its name and all chains and
#' backbone. Otherwise the components can be changed and the validation can 
#' be deleted.
#'
#' @param id Shiny module id
#' @param validation List of parameters to Initialze the UI.
#' @param fixed Defines whether  the user can change the parameters.
#'
#' @noRd 
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
                value = validation$adjustments,
                placeholder = "individual adjustment")
    ),
    createPriorityRadioButtons(ns("val_priority"),
                               validation$val_priority,
                               horizontal = TRUE),
    actionButton(ns("remove_validation"), label = NULL, icon = icon("trash"), class = "btn-remove-validation")
  )
}

#' @title Validation Server Function
#'
#' @description The Shiny Server module for one validation.
#' 
#' @details
#' Updates the validations from the UI. Can add components to the UI. If it
#' is fixed, gives all possible components back.
#' 
#' @param id Shiny module id
#' @param init_validation List of parameters to initialze the Server
#' @param chain_ids Chain IDs to choose from when selecting the components.
#' @param fixed Boolean to tell whether  the user can change the values.
#' @param on_remove Function to delete the validation.
#' @param input,output,session Internal parameters for {shiny}.
#' 
#' @return Returns a list with name, adjustments, val_priority and components
#' 
#' @noRd 
validationServer <- function(id, init_validation, chain_ids, fixed = FALSE, on_remove = NULL) {
  moduleServer(id, function(input, output, session) {
    
    val_iv <- shinyvalidate::InputValidator$new()

    val_iv$add_rule("adjustment", function(value) {
        if (is.null(input$adjustment))return(NULL)
        if (length(value) > 0 && !grepl("^$|^[-=]?([A-Z][a-z]?\\d*)([-=]?[A-Z][a-z]?\\d*)*$", value)) {
          return("Invalid adjustment")
        }        
        NULL
      })
    
    ns <- session$ns
    selection_ids  <- reactiveVal(list())
    selection_count <- reactiveVal(0)

    add_component_to_ui <- function(selection = NULL) {
      chain <- isolate(chain_ids())
      selection_count(selection_count() + 1)
      i <- selection_count()
      selection_ids(c(selection_ids(), i ))

      insertUI(
        selector = paste0("#", ns("selector_container")),
        where = "beforeEnd",
        ui = selectInput(
          ns(paste0("sel_", i)),
          label = NULL,
          choices = setNames(c("Backbone", chain),c("Backbone", seq_along(chain))),
          selected = selection
        ),
        immediate = TRUE
      )
      
    }
    isolate(for (com in init_validation$components) add_component_to_ui(com))
    
    observeEvent(input$add_component,{
      add_component_to_ui()
    })
    
    # keep the choices in sync when chains are added/removed
    observeEvent(chain_ids(), {
      
      req(length(selection_ids()) > 0)
      
      choices <- setNames(
        c("Backbone", chain_ids()),
        c("Backbone", seq_along(chain_ids()))
      )
      
      for (i in selection_ids()) {
        
        id <- paste0("sel_", i)
        cur <- input[[id]]

        updateSelectInput(
          session,
          id,
          choices = choices,
          selected = cur
        )
      }
    }, ignoreInit = TRUE)
    
    
    components <- reactive({
      if (fixed) {
        c("Backbone", chain_ids())
      } else {
        unlist(lapply(selection_ids(), function(i) input[[paste0("sel_", i)]]))
      }
    })
    
    observeEvent(input$remove_validation,{
      if(is.function(on_remove)) on_remove()
    }, ignoreInit = TRUE, once = TRUE)
    
    #Return of the function
    val_data <- reactive(list(
      name         = input$val_name %||% init_validation$name,
      adjustments  = input$adjustment,
      val_priority = input$val_priority %||% init_validation$val_priority,
      components   = components()
    ))
    list(
      val_data = val_data,
      validator = val_iv
    )
  })
}



#' Building the svg for a acrylchain
#'
#' @param c_min,c_max max defines the length of the chain, min 
#' defines when the chain gets grey (variabel length between them).
#' @param db_min,db_max db_max double bonds are distributed over the chain.
#' 
#' @noRd
build_chain_svg <- function(c_min, c_max, db_min = 0, db_max = 0) {
  c_max  <- max(c_max, c_min, 1)
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

#' Generate unique identifier
#'
#' @details
#' https://samfirke.com/2018/08/22/generating-unique-ids-using-r/
#' 
#' @noRd
create_unique_ids <- function(n = 1, char_len = 5){
  
  pool <- c(letters, LETTERS, 0:9)
  
  res <- character(n) # pre-allocating vector is much faster than growing it
  for(i in seq(n)){
    this_res <- paste0(sample(pool, char_len, replace = TRUE), collapse = "")
    while(this_res %in% res){ # if there was a duplicate, redo
      this_res <- paste0(sample(pool, char_len, replace = TRUE), collapse = "")
    }
    res[i] <- this_res
  }
  res
}


#' @title Default chain values
#'
#' @details
#' Variables and values:
#' \describe{
#'   \item{C_min}{Minimum number of carbon atoms for the chain. (initial 10)}
#'   \item{C_max}{Maximum number of carbon atoms for the chain. (initial 20)}
#'   \item{even_steps}{Determines if the chain only has C-values divided by 2. (initial TRUE)}
#'   \item{DB_min}{Minimum number of double bonds in the chain (initial 0)}
#'   \item{DB_max}{Maximum number of double bonds in the chain (initial 1)}
#'   \item{priority}{Which scoring priority: basic, important, required (initial basic)}
#'   \item{modifications}{List of modications added to the chain.}
#' }
#' @noRd
empty_chain <- list(
  C_min = 10,
  C_max = 20,
  even_steps = TRUE,
  DB_min = 0,
  DB_max = 1,
  priority = "basic",
  modifications = list()
)

#' @title Default Validation values
#'
#' @details
#' Variables and values:
#' \describe{
#'   \item{name}{Name of the validation}
#'   \item{adjustments}{String of what to change on the validation}
#'   \item{val_priority}{Which scoring priority: basic, important, required (initial basic)}
#'   \item{selectedComponents}{Selected Backbone/chains}
#' }
#' @noRd
empty_validation <- list(
  name = "enter name",
  adjustments = "",
  val_priority = "basic",
  selectedComponents = list()
)
