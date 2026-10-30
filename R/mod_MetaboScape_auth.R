
#' MetaboScape_auth Server Functions
#' @importFrom sodium password_verify

#' @noRd 
mod_MetaboScape_auth_server <- function(id, r){
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    show_login_modal <- function() {
      showModal(modalDialog(
        title = "Institute login for Bruker MetaboScape",
        div(class = "modal-flex",
            passwordInput(ns("password"), "Access code"),
            div(class = "small-selection", 
                selectInput(ns("metaboVersionSelection"), NULL,choices = c("2025b", "2026b"))) 
            ),
        tags$script(HTML(sprintf("
      $(document).off('keyup.modalLogin');
      $(document).on('keyup.modalLogin', function(e) {
        if (e.key === 'Enter') {
          $('#%s').click();
        }
      });
    ", ns("submit")))),
        tags$details(
          tags$summary(
            class = "validation-summary",
            textOutput(ns("versionText"))
          ),
          div(class = "card",
              passwordInput(ns("APIkey"), "API key"),
              passwordInput(ns("newPassword"), "New Password for API key"),
              textInput(ns("newServerAdress"), "New Server Adress", placeholder= "e.g. 10.20.0.123"),
              actionButton(ns("SetAPIkey"), "Set API Key")
          )
        
          ),
        footer = tagList(
          modalButton("Cancel"),
          actionButton(ns("submit"), "Log in")
        ),
        easyClose = FALSE
      ))
    }
    
    output$versionText <- renderText({
      paste("Set new API key for", r$MetaboVersion)
    })
    
    observeEvent(input$metaboVersionSelection,{
      r$MetaboVersion <- input$metaboVersionSelection
    })
    
    observeEvent(input$submit, {
      if (identical(Sys.getenv(paste0("APP_ACCESS_PASSWORD_HASH_", r$MetaboVersion)), "")){
        showNotification("APP_ACCESS_PASSWORD_HASH not set on this server", type = "error")
        return()
      }
      ok <- isTRUE(sodium::password_verify(
        Sys.getenv(paste0("APP_ACCESS_PASSWORD_HASH_", r$MetaboVersion)),
        input$password
      ))
      if (ok) {
        r$authenticated <- TRUE
        removeModal()
        if (!is.null(r$pending_action)) {
          action <- r$pending_action
          r$pending_action <- NULL
          action()
        }
      } else {
        showNotification("Wrong access code", type = "error")
      }
    })
    
    observeEvent(input$SetAPIkey, ({
      req(input$APIkey)
      version <- r$MetaboVersion
      do.call(Sys.setenv, as.list(setNames(input$APIkey, paste0("METABOSCAPE_API_KEY_", version))))
      do.call(Sys.setenv, as.list(setNames(sodium::password_store(input$newPassword), 
                          paste0("APP_ACCESS_PASSWORD_HASH_", version))))
      do.call(Sys.setenv, as.list(setNames(input$newServerAdress, paste0("METABOSCAPE_SERVER_ADRESS_", version))))
      
    }))
    
    # exposed to every other module via the shared `r` bus
    r$require_auth <- function(on_success) {
      if (isTRUE(r$authenticated)) {
        on_success()
      } else {
        r$pending_action <- on_success
        show_login_modal()
      }
    }
  })
}


#' @title Setting up MetaboScape connection
#' 
#' @description Uses a System variable to get the API key and connects to the
#' API to create rapiclient operations.
#'
#' @return The operation functions to interact with the MetaboScape API.
#'
#' @import rapiclient
#' @noRd
getMetaboScapeConnection <- function(version){
  print("Setting up API connection...")
  api_key <- Sys.getenv(paste0("METABOSCAPE_API_KEY_", version))
  if (identical(api_key, "")){
    showNotification("METABOSCAPE_API_KEY not set on this server", type = "error")
    return()
  }
  server_address <- Sys.getenv(paste0("METABOSCAPE_SERVER_ADRESS_", version))
  
  bruker_api <- rapiclient::get_api(url = paste0("http://", server_address, "/cxf/metaboscape/openapi.json"))
  bruker_api$basePath <- "/cxf/metaboscape"
  operations <- rapiclient::get_operations(bruker_api, .headers = c("api-key" = api_key))
  
}
    

