#' Rules UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_Rules_ui <- function(id) {
  ns <- NS(id)
  tagList(
 
  )
}
    
#' Rules Server Functions
#'
#' @noRd 
mod_Rules_server <- function(id){
  moduleServer(id, function(input, output, session){
    ns <- session$ns
 
  })
}
    
## To be copied in the UI
# mod_Rules_ui("Rules_1")
    
## To be copied in the server
# mod_Rules_server("Rules_1")
