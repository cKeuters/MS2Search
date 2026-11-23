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
    
    fluidRow(
      column(
        width = 12,
        div(
          class = "dreams-status-bar",
          uiOutput(ns("dreams_online")),
          textOutput(ns("dreams_online_msg")),
          span(
            class = "dreams-status-text",
            textOutput(ns("dreams_status"), inline = TRUE)
          )
        )
      )
    ),
    uiOutput(ns("embedding_ui"))
  )
}
    
#' DreaMS_Network Server Functions
#'
#' @importFrom plotly add_markers layout event_register event_data
#' @importFrom uwot umap 
#'
#' @noRd 
mod_DreaMS_Network_server <- function(id, Feature_keys, ms2_file_path_r, selected_feature_id_r, results_r){
  moduleServer(id, function(input, output, session){
    ns <- session$ns


    job_id      <- reactiveVal(NULL)
    job_state   <- reactiveVal(list(status = "idle"))
    embeddings  <- reactiveVal(NULL) 
    fails       <- reactiveVal(0)

    dreams_up <- reactive({ invalidateLater(15000, session); dreams_available() })
    output$dreams_online_msg <- renderText(if (dreams_up()) "DreaMS: available" else "DreaMS: not reachable")

    observeEvent(input$getEmbeddings, {
      if (!dreams_available()) {
        showNotification("DreaMS server is not reachable right now. Please try again later.", type = "error")
        return()
      }
      id <- tryCatch(dreams_submit(ms2_file_path_r()), 
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
    
    has_umap <- reactive({
      !is.null(embeddings()) && nrow(embeddings()) > 0
    })
    
    output$embedding_ui <- renderUI({

      if(has_umap()){
        div(
          div(
            class = "dreams-controls",
            radioButtons(
              ns("color_mode"),
              "Colour by",
              choices = c(
                "Spectra only" = "none",
                "Score" = "score",
                "Filter" = "filter"
              ),
              inline = TRUE
            ),
            conditionalPanel(
              condition = "input.color_mode == 'filter'",
              ns = ns,
              selectInput(
                ns("filters"),
                "Filters",
                choices = c("None"),
                multiple = TRUE
              )
            )
          ),
          
          
          plotly::plotlyOutput(ns("umap"), height = "600px"),
          
          div(
            style = "margin-top:15px;",
            actionButton(ns("getEmbeddings"), "Get Embeddings", class = "dreams-btn")
          )
        )
      }
      else{
        div(
          style = "text-align:center; margin-top:150px;",
          actionButton(ns("getEmbeddings"), "Get Embeddings", class = "dreams-btn")
        )
      }
    })
    
    output$dreams_online <- renderUI({
      
      online <- dreams_up()
      cls <- if (online)"dreams-dot dreams-online"
      else "dreams-dot dreams-offline"
      
      txt <- if (online)"DreaMS: online"
      else "DreaMS: offline"
      
      tags$span(
        class = cls,
        title = txt
      )
    })
    
    #---------------------------------------------------------------------------
    # Visualization of DreaMS
    #---------------------------------------------------------------------------
    
    # data for the plot
    plot_df <- reactive({
      req(embeddings())
      
      coords <- uwot::umap(
        embeddings(),
        n_neighbors = 5,
        min_dist    = 0.1,
        metric      = "cosine",
        n_threads   = 4
      )
      embed_df <- data.frame(id = Feature_keys, x = coords[, 1], y = coords[, 2])
      
      res <- results_r()
      if (!is.null(res) && nrow(res) > 0) {
        embed_df <- dplyr::left_join(
          embed_df,
          dplyr::select(res, id = FEATURE_ID, score),
          by = "id"
        )
      }

      
      embed_df
    })
    
    output$umap <- plotly::renderPlotly({
      req(plot_df())
      
      df <- plot_df()
      df$hover <- paste0(df$id, if ("score" %in% names(df)) paste0("<br>Score: ", round(df$score, 3)) else "")
      
      p <- switch(
        input$color_mode,
        none = plotly::plot_ly(
          df, x = ~x, y = ~y, key = ~id, text = ~hover, hoverinfo = "text",
          source = "umap", type = "scattergl", mode = "markers",
          marker = list(color = "grey60", size = 6, opacity = 0.7)
        ),
        score = plotly::plot_ly(
          df, x = ~x, y = ~y, key = ~id, text = ~hover, hoverinfo = "text",
          source = "umap", type = "scattergl", mode = "markers",
          color = ~score, colors = "viridis",
          marker = list(size = 6, opacity = 0.8)
        ),
        filter = {
            # no filter chosen yet -> neutral
            plotly::plot_ly(df, x = ~x, y = ~y, key = ~id, text = ~hover,
                            hoverinfo = "text", source = "umap",
                            type = "scattergl", mode = "markers",
                            marker = list(color = "grey60", size = 6))
        }
      )
      
      # highlight the currently selected spectrum
      sel <- selected_feature_id_r()
      if (!is.null(sel) && sel %in% df$id) {
        s <- df[df$id == sel, ]
        p <- plotly::add_markers(
          p, data = s, x = ~x, y = ~y, key = NULL, inherit = FALSE,
          showlegend = FALSE, hoverinfo = "skip",
          marker = list(size = 14, color = "rgba(0,0,0,0)",
                        line = list(color = "black", width = 2))
        )
      }
      
      p |>
        plotly::layout(
          xaxis = list(title = "UMAP 1", zeroline = FALSE),
          yaxis = list(title = "UMAP 2", zeroline = FALSE),
          uirevision = "keep"   # keeps zoom/pan when colour mode or selection changes
        ) |>
        plotly::event_register("plotly_click")
      
      
      # click on a dot -> set selection
      observeEvent(plotly::event_data("plotly_click", source = "umap"), {
        ed <- plotly::event_data("plotly_click", source = "umap")
        if (!is.null(ed$key)) selected_feature_id_r(ed$key[[1]])
      })
      
      p
    })
    
  
    
  })
}

#-------------------------------------------------------------------------------
# DreaMS Container communication
#-------------------------------------------------------------------------------

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