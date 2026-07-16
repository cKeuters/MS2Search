#-------------------------------------------------------------------------------
# functions for analysis screen
#-------------------------------------------------------------------------------
createPriorityRadioButtons <- function(id, horizontal = TRUE){
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
    next_entry_id = 2,
    min_repeat = 0,
    max_repeat = 1,
    just_after_precursor = FALSE,
    sumOfSmall = FALSE,
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

#-------------------------------------------------------------------------------
# Save and Load parameters
#-------------------------------------------------------------------------------

save_parameters <- function(state, cards, ppm1, ppm2, project_path){
  config <- list(
    state = state,
    rv = rv,
    ppm1  = ppm1,
    ppm2  = ppm2
  )
  
  qsave(config, paste0(project_path,"parameters.qs"))
  
}

load_parameters <- function(project_path){
  config <- qread( paste0(project_path,"parameters.qs"))
  
  for (nm in names(config$state)) {
    state[[nm]] <- config$state[[nm]]
  }
  for (nm in names(config$rv)) {
    rv[[nm]] <- config$rv[[nm]]
  }
  updateNumericInput(session, "ppm1", value = config$ppm1)
  updateNumericInput(session, "ppm2", value = config$ppm2)
}



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

