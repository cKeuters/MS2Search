#' ui 
#'
#' @description A utils function
#'
#' @noRd

#' @description Title and buttons for Saving Project and Return to Setup in the topbar
#' @noRd
create_topbar_InProject <- function(){
  div(class = "top-bar",
      #h2(class = "app-title","MS2Search - Analysis"),
      uiOutput("state_screen"),
      div(class= "top-bar-nav",
        actionButton("back_to_welcome", "Setup",
                     icon = icon("house"), class = "btn-default"),
        actionButton("back_to_analysis", "Analysis",
                     icon = icon("sliders"), class = "btn-default"),
        conditionalPanel(
          condition = "output.has_results",
          actionButton("back_to_results", "Results",
                     icon = icon("square-poll-horizontal"), class = "btn-default"))
      ),
      actionButton("save_project", "Save Project",
                   icon = icon("folder-open"), class = "btn-default")
  )
}

#' @description Title and button for Loading Project in the topbar
#' @noRd
create_topbar_InSetup <- function(){
  div(class = "top-bar",
      h2(class = "app-title", "Cera Search Project - Setup"),
      div(id = "load_project_wrapper",
          fileInput("load_project_file", label = NULL, accept = ".zip",
                    buttonLabel = tagList(icon("folder-open"), "Load Project"),
                    placeholder = "")
      )
  )
}

#' @description Legend of Priorities
#' @noRd
createAnalysisLegend <- function(){
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
  )
}

#' @description Priority Radio Buttons
#' @noRd
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

#' @description Tooltips
#' @noRd
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

