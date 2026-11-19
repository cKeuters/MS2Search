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
    actionButton(ns("getEmbeddings"), "get Embeddings"),
    textOutput(ns("dreams_online")),
    textOutput(ns("dreams_status"))
 
  )
}
    
#' DreaMS_Network Server Functions
#'
#' @noRd 
mod_DreaMS_Network_server <- function(id, ms2_file_path){
  moduleServer(id, function(input, output, session){
    ns <- session$ns


    job_id      <- reactiveVal(NULL)
    job_state   <- reactiveVal(list(status = "idle"))
    embeddings  <- reactiveVal(NULL) 
    fails       <- reactiveVal(0)

    dreams_up <- reactive({ invalidateLater(15000, session); dreams_available() })
    output$dreams_online <- renderText(if (dreams_up()) "DreaMS: available" else "DreaMS: not reachable")

    observeEvent(input$getEmbeddings, {
      if (!dreams_available()) {
        showNotification("DreaMS server is not reachable right now. Please try again later.", type = "error")
        return()
      }
      id <- tryCatch(dreams_submit(ms2_file_path()), 
                    error = function(e) {
                      showNotification(paste("Could not start DreaMS job:", conditionMessage(e)), type = "error")
                      NULL
                    })
      req(id)
      fails(0)
      job_id(id)
      job_state(list(status = "queued", ahead = 0))
    })

    observe({
      id <- job_id()
      req(id)
      invalidateLater(5000, session)

      st <- tryCatch(dreams_status(id), error = function(e) list(status = "unreachable"))
      job_state(st)

      if (st$status == "unreachable") {
        fails(isolate(fails()) + 1)
        if (isolate(fails()) >= 12) {   # about one minute without an answer
          showNotification("Lost connection to DreaMS. Please start again.", type = "error", duration = NULL)
          job_id(NULL)
        }
        return()
      }
      fails(0)

      if (st$status == "done") {
        embeddings(tryCatch(dreams_result(id), error = function(e) {
          showNotification(conditionMessage(e), type = "error"); NULL
        }))
        job_id(NULL)
      } else if (st$status %in% c("error", "lost")) {
        showNotification("DreaMS job failed or was lost (server restarted?). Please start again.",
                        type = "error", duration = NULL)
        job_id(NULL)
      }
    })

    output$dreams_status <- renderText({
      st <- job_state()
      switch(st$status,
        queued      = paste0("Waiting for DreaMS (", st$ahead, " job(s) ahead) ..."),
        running     = paste0("DreaMS is computing ... ", st$elapsed %/% 60, " min elapsed"),
        done        = "DreaMS embeddings ready.",
        error       = paste("DreaMS failed:", st$error),
        lost        = "DreaMS was restarted, job lost. Please start again.",
        unreachable = "DreaMS temporarily unreachable, retrying ...",
        ""
      )
    })
  })
}
    

dreams_req <- function(path, timeout = 5) {
  httr2::request(Sys.getenv("DREAMS_URL", "http://localhost:8000")) |>
    httr2::req_url_path(path) |>
    httr2::req_timeout(timeout) |>
    httr2::req_error(is_error = function(resp) FALSE)  # we check the status ourselves
}

dreams_available <- function() {
  tryCatch(
    httr2::resp_status(httr2::req_perform(dreams_req("/health", 2))) == 200,
    error = function(e) FALSE
  )
}

dreams_submit <- function(mgf_path) {
  resp <- dreams_req("/jobs", timeout = 120) |>
    httr2::req_body_multipart(file = curl::form_file(mgf_path)) |>
    httr2::req_perform()
  if (httr2::resp_status(resp) != 200) stop("DreaMS rejected the job (HTTP ", httr2::resp_status(resp), ")")
  httr2::resp_body_json(resp)$job_id
}

dreams_status <- function(job_id) {
  resp <- httr2::req_perform(dreams_req(paste0("/jobs/", job_id), 3))
  if (httr2::resp_status(resp) == 404) return(list(status = "lost"))
  httr2::resp_body_json(resp)
}

dreams_result <- function(job_id) {
  resp <- httr2::req_perform(dreams_req(paste0("/jobs/", job_id, "/result"), 120))
  if (httr2::resp_status(resp) != 200) stop("Could not fetch result (HTTP ", httr2::resp_status(resp), ")")
  n <- as.integer(httr2::resp_header(resp, "X-Rows"))
  p <- as.integer(httr2::resp_header(resp, "X-Cols"))
  matrix(readBin(httr2::resp_body_raw(resp), "numeric", n = n * p, size = 4),
         nrow = n, ncol = p, byrow = TRUE)
}