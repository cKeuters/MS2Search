#' About UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_about_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h2("About"),
    p("Short description of what the app does."),
    tags$dl(
      tags$dt("Version"),
      tags$dd(textOutput(ns("version"), inline = TRUE)),
      tags$dt("Authors"),
      tags$dd("Carina Keuters, IPB"),
      tags$dt("License"),
      tags$dd("MIT"),
      tags$dt("Source code"),
      tags$dd(tags$a(href = "https://github.com/cKeuters/MS2Search",
                     "GitHub: https://github.com/cKeuters/MS2Search", target = "_blank"))
    )
  )
}


    
#' About Server Functions
#'
#' @importFrom utils packageVersion
#' @noRd 
mod_about_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    output$version <- renderText({
      as.character(utils::packageVersion(golem::get_golem_options("pkg_name") %||%
                                           "MS2Search"))
    })
  })
}