#' Rules UI Function
#'
#' @description A shiny Module.
#'
#' @param id, card, metabolite_choices, input,output,session Internal parameters for {shiny}.
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_Rules_ui <- function(id, card, metabolite_choices) {
  ns <- NS(id)
  card_class <- "peak-card"

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
              div(class = "flex-distance",
                  div(class = "toggle-switch",
                   tags$input(type = "checkbox", id = ns("peakToggle"), checked = if (isTRUE(card$is_mzPeak)) "checked" else NULL),
                   div(class = "toggle-slider"),
                   div(class = "toggle-labels",
                       span(class = "toggle-label-left", "Neutral loss"),
                       span(class = "toggle-label-right", "m/z peak")
                   )
              )),
              
              # --- Mass entries ---
              selectInput(ns("input_type"), "Input type",
                          choices = c("Exact mass" = "mass", "Metabolite name" = "metabolite"),
                          selected = card$input_type),
              uiOutput(ns("entries_ui")),
              
              actionButton(ns("add_entry"), "Add value", icon = icon("plus"), class = "btn-sm"),
              uiOutput(ns("mass_preview")),
              
              # --- UI for "Neutral loss" (unchecked) ---
              conditionalPanel(
                condition = "!input.peakToggle",
                ns = ns,
                div(
                  checkboxInput(ns("sequential"), "Sequential neutral losses", value = card$sequential),
                  uiOutput(ns("sequential_slider_ui")),
                  checkboxInput(ns("precursor"), "Just after precursor", value = card$just_after_precursor)
                )
              ),
 
              checkboxInput(ns("sumOfSmall"), "Validation by sum of smaller peaks", value = card$sumOfSmall),
              div(style = "text-align: center;", createPriorityRadioButtons(ns("priority"), "basic"))
              
          )
      )
  )
}
    
#' Rules Server Functions
#'
#' @noRd 
mod_Rules_server <- function(id, initial, on_delete, on_duplicate, metabolite_choices, parent_session ) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    rule_state <- reactiveValues(
      id = initial$id, 
      label = initial$label, 
      is_mzPeak = initial$is_mzPeak,
      input_type = initial$input_type,
      min_repeat = initial$min_repeat, 
      max_repeat = initial$max_repeat,
      just_after_precursor = initial$just_after_precursor, 
      sumOfSmall = initial$sumOfSmall, 
      sequential = initial$sequential,
      next_card_id = initial$next_card_id
    )
    
    
    #--------------------------------------------------------------------------
    # Entry modules
    #--------------------------------------------------------------------------
    
    entry_modules <- reactiveValues()
    entry_ids <- reactiveVal(character(0))
    next_entry_id <- reactiveVal(initial$next_entry_id)
    
    output$entries_ui <- renderUI({
      ids <- entry_ids()
      if (length(ids) == 0) return(NULL)
      
      lapply(ids, function(eid) {
        init <- isolate(entry_modules[[eid]]())
        entryModuleUI(ns(paste0("entry_", eid)), init, rule_state$input_type, metabolite_choices)
      })
    })
    
    remove_entry <- function(eid) {
      eid <- as.character(eid)
      if (length(entry_ids()) <= 1) return(invisible())  # mind. 1 Zeile behalten
      removeUI(selector = paste0("#", ns(paste0("entry_wrap_", eid))))
      entry_modules[[eid]] <- NULL
      entry_ids(setdiff(entry_ids(), eid))
    }
    
    add_entry <- function(entry_init) {
      message("add_entry called with: ", entry_init$id, " / ", entry_init$value)
      eid <- as.character(entry_init$id)
      entry_modules[[eid]] <- entryModuleServer(paste0("entry_", eid), entry_init, remove_entry)
      entry_ids(c(entry_ids(), eid))
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
    
        current_entries <- reactive({
      lapply(entry_ids(), function(eid) entry_modules[[eid]]())
    })
    
    #--------------------------------------------------------------------------
    # Toggle
    #--------------------------------------------------------------------------
    
    output$sequential_slider_ui <- renderUI({
      req(input$sequential)
      noUiSliderInput(
        inputId = ns("slider"),
        label = "Number of losses",
        min = 0,
        max = 10,
        value = c(rule_state$min_repeat, rule_state$max_repeat),
        step = 1,
        format = wNumbFormat(decimals = 0),
        tooltips = TRUE,
        pips = list(mode = "values", values = list(0, 10), density = 10),
        color = "#337ab7"
      )
    })
    
    observeEvent(input$peakToggle, {
      rule_state$is_mzPeak <- isTRUE(input$peakToggle)
    }, ignoreNULL = FALSE, ignoreInit = TRUE)
  

    
    
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
    
    
    output$mass_preview <- renderUI({
      masses <- vapply(current_entries(), function(e) {
        resolve_mass(list(input_type = rule_state$input_type, value = e$value))
      }, numeric(1))
      masses <- abs(as.numeric(masses[!is.na(masses)]))
      if (length(masses) == 0) return(NULL)
      if(!rule_state$is_mzPeak)masses <- - masses
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
      list(id = rule_state$id, 
           label = rule_state$label, 
           is_mzPeak = rule_state$is_mzPeak,
           input_type = rule_state$input_type,
           entries = current_entries(), 
           next_entry_id = next_entry_id(),
           min_repeat = rule_state$min_repeat, 
           max_repeat = rule_state$max_repeat,
           just_after_precursor = rule_state$just_after_precursor, 
           sumOfSmall = rule_state$sumOfSmall, 
           sequential = rule_state$sequential, 
           next_card_id = rule_state$next_card_id)
    })
  }, session = parent_session)
}

#' Metabolite Entry UI Functions
#'
#' @noRd 
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

#' Metabolite Entry Server Functions
#'
#' @noRd 
entryModuleServer <- function(id, initial, on_remove) {
  moduleServer(id, function(input, output, session) {
    
    state <- reactiveValues(
      id = initial$id, 
      value = initial$value)
    
    observeEvent(input$value, {
      message(paste("observeEvent(input$value",input$value))
      state$value <- input$value
    }, ignoreInit = TRUE)
    
    observeEvent(input$remove, on_remove(initial$id), ignoreInit = TRUE)
    reactive(list(id = state$id, value = state$value))
  })
}

#' @description Calculates the mass of metabolite names
#' @noRd
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


METABOLITES <- data.frame(
  name = c(
    "Glucose Loss", "Fructose", "Sucrose", "Lactate", "Pyruvate",
    "Citrate", "Malate", "ATP", "ADP", "AMP"
  ),
  mass = c(
    -180.0634, 180.0634, 342.1162, 89.0244, 87.0082,
    192.0270, 133.0137, 506.9957, 427.0294, 347.0631
  ),
  stringsAsFactors = FALSE
)

empty_entry <- function(id) list(id = id, value = "")

empty_card <- function(id) {
  list(
    id = id,
    label = paste("Peak filter", id),
    is_mzPeak = FALSE,
    input_type = "mass",
    entries = list(empty_entry(1)),
    next_entry_id = 2,
    min_repeat = 0,
    max_repeat = 1,
    just_after_precursor = FALSE,
    sumOfSmall = FALSE,
    sequential = FALSE,
    priority = 1,
    next_card_id = NA
    )
}


