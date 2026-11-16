#' DreaMS_Network UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_DreaMS_Network_ui <- function(id) {
  ns <- NS(id)
  tagList(
    actionButton(ns("getEmbeddings"), "get Embeddings")
 
  )
}
    
#' DreaMS_Network Server Functions
#'
#' @noRd 
mod_DreaMS_Network_server <- function(id, ms2_file_path){
  moduleServer(id, function(input, output, session){
    ns <- session$ns
 
    observeEvent(input$getEmbeddings,{
      print("Button clicked")
      mgf_path <- ms2_file_path()
      if(mgf_path != "" && !is.null(mgf_path))get_dreams_embeddings(mgf_path)
      
    })
  })
}
    
#' @title Getting DreaMS Embeddings 
#'
#' @description Requesting the dreams embeddings from the dreams container.
#' 
#' 
#' @param mgf_path Holds the filepath of a mgf-file.
#'
#' @return Returns a matrix of the embeddings.
#'
#' @importFrom httr2 request req_url_path req_timeout req_perform resp_header
#' @importFrom curl form_file
#' @noRd
get_dreams_embeddings <- function(mgf_path) {
  print("Getting mgf path")
  resp <- httr2::request(Sys.getenv("DREAMS_URL", "http://localhost:8000")) |>
    httr2::req_url_path("/embed") |>
    httr2::req_body_multipart(file = curl::form_file(mgf_path)) |>
    httr2::req_timeout(600) |>
    httr2::req_perform()
  
  n <- as.integer(httr2::resp_header(resp, "X-Rows"))
  print(paste("n",n))
  
  p <- as.integer(httr2::resp_header(resp, "X-Cols"))
  print(paste("p",p))
  
  raw <- httr2::resp_body_raw(resp)
  matrix(readBin(raw, "numeric", n = n * p, size = 4), nrow = n, ncol = p, byrow = TRUE)
}
