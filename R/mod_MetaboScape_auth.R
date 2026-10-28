
#' MetaboScape_auth Server Functions
#' @importFrom sodium password_verify

#' @noRd 
mod_MetaboScape_auth_server <- function(id, r){
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    show_login_modal <- function() {
      showModal(modalDialog(
        title = "Institute login",
        passwordInput(ns("password"), "Access code"),
        tags$script(HTML(sprintf("
      $(document).off('keyup.modalLogin');
      $(document).on('keyup.modalLogin', function(e) {
        if (e.key === 'Enter') {
          $('#%s').click();
        }
      });
    ", ns("submit")))),
        
        div(class = "card",
          passwordInput(ns("APIkey"), "API key"),
          passwordInput(ns("newPassword"), "New Password for API key"),
          actionButton(ns("SetAPIkey"), "Set API Key")
        ),
        footer = tagList(
          modalButton("Cancel"),
          actionButton(ns("submit"), "Log in")
        ),
        easyClose = FALSE
      ))
    }
    
    observeEvent(input$submit, {
      if (identical(Sys.getenv("APP_ACCESS_PASSWORD_HASH"), "")){
        showNotification("APP_ACCESS_PASSWORD_HASH not set on this server", type = "error")
        return()
      }
      ok <- isTRUE(sodium::password_verify(
        Sys.getenv("APP_ACCESS_PASSWORD_HASH"),
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
      Sys.setenv(METABOSCAPE_API_KEY = input$APIkey)
      Sys.setenv(APP_ACCESS_PASSWORD_HASH = sodium::password_store(input$newPassword))
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
getMetaboScapeConnection <- function(){
  print("Setting up API connection...")
  api_key <- Sys.getenv("METABOSCAPE_API_KEY")
  if (identical(api_key, "")){
    showNotification("METABOSCAPE_API_KEY not set on this server", type = "error")
    return()
  }
  bruker_api <- rapiclient::get_api(url = "http://ServerAddress/cxf/metaboscape/openapi.json")
  print("Got Bruker API")
  bruker_api$basePath <- "/cxf/metaboscape"
  print("Bruker API basePath set")
  operations <- rapiclient::get_operations(bruker_api, .headers = c("api-key" = api_key))
  
}
    

