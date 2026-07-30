
library(shiny)
library(shinyFiles)


fluidPage(
  shinyjs::useShinyjs(),
    includeCSS("custom.css"),
    uiOutput("main_screen")
    
  
)
