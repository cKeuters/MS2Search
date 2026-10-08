#' @title Rules UI Function
#'
#' @description The Shiny UI module for one Rule card.
#' 
#' @details
#' Takes the values from the card parameter and builds the UI to make them 
#' adjustable.
#'
#' @param id Shiny module id
#' @param card Holds the informations of one rule card from the UI.
#' @param metabolite_choices Dataframe of Metabolites and the possible masses.
#'
#' @importFrom shiny NS tagList 
#' @noRd 
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
              div(id = ns("entries_ui")),
              
              actionButton(ns("add_entry"), "Add value", icon = icon("plus"), class = "btn-sm"),
              uiOutput(ns("mass_preview")),
              
              # --- UI for "Neutral loss" (unchecked) ---
              conditionalPanel(
                condition = "!input.peakToggle",
                ns = ns,
                div(
                  checkboxInput(ns("sequential"), "Sequential neutral losses", value = card$sequential),
                  uiOutput(ns("sequential_slider_ui")),
                  checkboxInput(ns("precursor"), "Just after precursor", value = card$just_after_precursor),
                  div(class = "flex-distance card-label",
                      numericInput(ns("rule_mz_min"), "Minimum m/z value", value = card$rule_mz_min),
                      numericInput(ns("rule_mz_max"), "Maximum m/z value", value = card$rule_mz_max))
                )
              ),
 
              checkboxInput(ns("filter_noise"), "Remove noise from MS2 spectrum before search", value = card$filter_noise),
              div(style = "text-align: center;", createPriorityRadioButtons(ns("rule_priority"), card$rule_priority))
              
          )
      )
  )
}
    
#' @title Rules Server Function
#'
#' @description The Shiny Server module for one Rule card.
#' 
#' @details
#' The Rule Server has a rule_state reactive with the rule parameters. It observes
#' the value changes in the UI server to change the rule_state.
#' In addition to that it coordinates the modification modules with the ability 
#' to add and remove them from the UI. It also renders the output of the mass
#' preview.
#' 
#' @param id Shiny module id
#' @param initial Holds the informations of one rule.
#' @param polarity Polarity of the dataset.
#' @param on_delete Function to delete the card.
#' @param on_duplicate Function to duplicate the card.
#' @param metabolite_choices Dataframe of Metabolites and the possible masses.
#' @param input,output,session Internal parameters for {shiny}.
#' 
#' @return Returns a list with label, is_mzPeak, input_type, mass entries, 
#' next_entry_id, min_repeat, max_repeat, just_after_precursor, filter_noise, 
#' sequential, next_card_id, rule_priority, rule_mz_min, rule_mz_max.
#' 
#' @noRd 
mod_Rules_server <- function(id, initial, polarity, on_delete, on_duplicate, metabolite_choices ) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # ReactiveValues of the Rule parameters
    #
    # @details
    # Variables and values:
    # \describe{
    #   \item{id}{The id of this rule for deleting purposes.}
    #   \item{label}{Name of the rule (initial paste("Peak filter", id))}
    #   \item{is_mzPeak}{Boolean if it is a neutral loss or a peak. (initial FALSE)}
    #   \item{input_type}{Input as masses or as names (initial mass)}
    #   \item{min_repeat}{Minumum number of sequential repeats (initial 2)}
    #   \item{max_repeat}{Maximum number of sequential repeats (initial 3)}
    #   \item{just_after_precursor}{Defines if neutral losses are only searched after the precursor. (initial FALSE)}
    #   \item{filter_noise}{Boolean defines if noise is filtered in ms2 (initial TRUE)}
    #   \item{sequential}{If sequential neutral losses are searched. (initial FALSE)}
    #   \item{next_card_id}{id of the next card (initial 2)}
    #   \item{rule_priority}{Which scoring priority: basic, important, required (initial basic)}
    #   \item{rule_mz_min}{Minimal mz value in the ms2 to search for neutral losses (initial NULL)}
    #   \item{rule_mz_max}{Maximum mz value in the ms2 to search for neutral losses (initial NULL)}
    # }
    rule_state <- reactiveValues(
      id = initial$id, 
      label = initial$label, 
      is_mzPeak = initial$is_mzPeak,
      input_type = initial$input_type,
      min_repeat = initial$min_repeat, 
      max_repeat = initial$max_repeat,
      just_after_precursor = initial$just_after_precursor, 
      filter_noise = initial$filter_noise, 
      sequential = initial$sequential,
      next_card_id = initial$next_card_id,
      rule_priority = initial$rule_priority,
      rule_mz_min = initial$rule_mz_min,
      rule_mz_max = initial$rule_mz_max
    )
    
    
    #--------------------------------------------------------------------------
    # Entry modules (mass/molecule name)
    #--------------------------------------------------------------------------
    
    entry_modules <- reactiveValues()
    entry_ids <- reactiveVal(character(0))
    next_entry_id <- reactiveVal(initial$next_entry_id)
    
    remove_entry <- function(eid) {
      eid <- as.character(eid)
      if (length(entry_ids()) <= 1) return(invisible()) 
      removeUI(selector = paste0("#", ns(paste0("entry_wrap_", eid))))
      entry_modules[[eid]] <- NULL
      entry_ids(setdiff(entry_ids(), eid))
    }
    
    add_entry <- function(entry_init, defer_ui = FALSE) {
      eid <- as.character(entry_init$id)
      entry_modules[[eid]] <- entryModuleServer(paste0("entry_", eid), entry_init, remove_entry)
      entry_ids(c(entry_ids(), eid))
      
      insert_entry_ui <- function() {
        insertUI(
          selector = paste0("#", ns("entries_ui")),
          where = "beforeEnd",
          ui = div(
            id = ns(paste0("entry_wrap_", eid)),
            entryModuleUI(ns(paste0("entry_", eid)), entry_init, isolate(rule_state$input_type), metabolite_choices)
          )
        )
      }
      
      if (defer_ui) {
        session$onFlushed(insert_entry_ui, once = TRUE)
      } else {
        insert_entry_ui()
      }
    }
    
    #adds all entries which come from the saved data
    isolate({ for (e in initial$entries) add_entry(e, defer_ui = TRUE) })
    
   #gets the current mass entries value lists from their servers
    current_entries <- reactive({
      lapply(entry_ids(), function(eid) entry_modules[[eid]]())
    })
    
    #--------------------------------------------------------------------------
    # ObserveEvents & Toggle UI
    #--------------------------------------------------------------------------
    
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
    observeEvent(input$rule_priority, rule_state$rule_priority <- input$rule_priority, ignoreInit = TRUE)
    observeEvent(input$precursor, rule_state$just_after_precursor <- isTRUE(input$precursor), ignoreInit = TRUE)
    observeEvent(input$filter_noise, rule_state$filter_noise <- isTRUE(input$filter_noise), ignoreInit = TRUE)
    observeEvent(input$sequential, rule_state$sequential <- isTRUE(input$sequential), ignoreInit = TRUE)
    observeEvent(input$rule_mz_min, rule_state$rule_mz_min <- input$rule_mz_min)
    observeEvent(input$rule_mz_max, rule_state$rule_mz_max <- input$rule_mz_max)
    
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
    
    #Takes the masses from the entries and builds UI badges for them
    output$mass_preview <- renderUI({
      masses <- vapply(current_entries(), function(e) {
        resolve_mass(list(input_type = rule_state$input_type, value = e$value, is_mzPeak = rule_state$is_mzPeak, polarity = polarity))
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
           filter_noise = rule_state$filter_noise, 
           sequential = rule_state$sequential, 
           next_card_id = rule_state$next_card_id,
           rule_priority = rule_state$rule_priority,
           rule_mz_min = rule_state$rule_mz_min,
           rule_mz_max = rule_state$rule_mz_max)
    })
  })
}

#' @title Entry UI Function
#'
#' @description The Shiny UI module for one entry of mass or metabolite name.
#' 
#' @details
#' Is a textinput or a selectinput based on the inputtype to adjust the 
#' entry mass or metabolite name.
#'
#' @param id Shiny module id
#' @param entry The mass or the metabolite name.
#' @param input_type Defines weather it is masses or metabolite name.
#' @param metabolite_choices Dataframe of Metabolites and the possible masses.
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

#' @title Entry Server Function
#'
#' @description The Shiny Server module for one entry of mass or metabolite name.
#' 
#' @details
#' The Rule Server has the id and observes the value of the entry. It can delete
#' by the id.
#' 
#' @param id Shiny module id
#' @param initial Holds the informations of one entry
#' @param on_remove Function to delete the entry
#' @param input,output,session Internal parameters for {shiny}.
#' 
#' @return Returns a list with the id and the value.
#' 
#' @noRd 
entryModuleServer <- function(id, initial, on_remove) {
  moduleServer(id, function(input, output, session) {
    
    state <- reactiveValues(
      id = initial$id, 
      value = initial$value)
    
    observeEvent(input$value, {
      state$value <- input$value
    }, ignoreInit = TRUE)
    
    observeEvent(input$remove, on_remove(initial$id), ignoreInit = TRUE)
    
    reactive(list(id = state$id, value = state$value))
  })
}




#' @title Resolve metabolite name masses 
#'
#' @description Calculates the mass of metabolite names
#' 
#' @details
#' Uses the dataframe METABOLITES to get the mass from a neutral loss or the 
#' mz peak of the chosen entries.
#' 
#' @param card_info Rule card parameters.
#' @param polarity Polarity of the dataset.
#'
#' @return Returns the masses of the metabolite names.
#'
#' @noRd
resolve_mass <- function(card_info, polarity) {
  #adducts have to be added (at the moment only H+ and H-)
  
  if (card_info$input_type == "mass") {
    m <- suppressWarnings(as.numeric(card_info$value))
    if (!is.na(m)) return(round(m, 4))
    return(NA)
  }
  idx <- which(METABOLITES$name == card_info$value)
  if(card_info$is_mzPeak){
    if (length(idx) == 1){
      monoMass <- METABOLITES$monoMass[idx]
      ifelse(card_info$polarity == "pos",monoMass <- monoMass + 1.008, monoMass <- monoMass - 1.008 )
      return(monoMass)
    }
  }
  else{
    if (length(idx) == 1) return(METABOLITES$NeutralLoss[idx])
  }
  
  NA
}

#' Metabolite names and masses
#'
#' A dataframe containing common metabolites with their neutral loss and
#' monoisotopic masses.
#'
#' @format A data frame with 9 obs. of  3 variables:
#' \describe{
#'   \item{name}{metabolite name}
#'   \item{NeutralLoss}{value of the loss}
#'   \item{monoMass}{monoisotopic mass (Da)}
#' }
METABOLITES <- data.frame(
  name = c(
    "Glucose", "Fructose", "Sucrose", "Lactate", "Pyruvate",
    "Citrate", "Malate", "Water", "Phosphate"
  ),
  NeutralLoss = c(
    162.0528, 162.0528, 162.0528, 72.0211, 43.9898,
    174.0164, 116.0110, 18.0106, 97.9769
  ),
  # monoisotopic mass (Da)
  monoMass = c(
    180.0634, 180.0634, 342.1162, 90.0317, 88.0160,
    192.0270, 134.0215, 18.0106, 97.9769
  ),
  stringsAsFactors = FALSE
)

#' @title Initial values of a new entry
#' 
#' @details
#' Variables and values:
#' \describe{
#'   \item{id}{The id of this entry for deleting purposes.}
#'   \item{value}{mass or metabolite name}
#' }
#' @param id id of the entry
empty_entry <- function(id) list(id = id, value = "")



#' @title Initial values of a new rule
#' 
#' @details
#' Variables and values:
#' \describe{
#'   \item{id}{The id of this rule for deleting purposes.}
#'   \item{label}{Name of the rule (initial paste("Peak filter", id))}
#'   \item{is_mzPeak}{Boolean if it is a neutral loss or a peak. (initial FALSE)}
#'   \item{input_type}{Input as masses or as names (initial mass)}
#'   \item{min_repeat}{Minumum number of sequential repeats (initial 2)}
#'   \item{max_repeat}{Maximum number of sequential repeats (initial 3)}
#'   \item{just_after_precursor}{Defines if neutral losses are only searched after the precursor. (initial FALSE)}
#'   \item{filter_noise}{Boolean defines if noise is filtered in ms2 (initial TRUE)}
#'   \item{sequential}{If sequential neutral losses are searched. (initial FALSE)}
#'   \item{next_card_id}{id of the next card (initial 2)}
#'   \item{rule_priority}{Which scoring priority: basic, important, required (initial basic)}
#'   \item{rule_mz_min}{Minimal mz value in the ms2 to search for neutral losses (initial NULL)}
#'   \item{rule_mz_max}{Maximum mz value in the ms2 to search for neutral losses (initial NULL)}
#' }
#' 
#' @param id id of the Rule card
empty_card <- function(id) {
  list(
    id = id,
    label = paste("Peak filter", id),
    is_mzPeak = FALSE,
    input_type = "mass",
    entries = list(empty_entry(1)),
    next_entry_id = 2,
    min_repeat = 2,
    max_repeat = 3,
    just_after_precursor = FALSE,
    sequential = FALSE,
    rule_priority = "basic",
    rule_mz_min = NULL,
    rule_mz_max = NULL,
    filter_noise = TRUE
    )
}


