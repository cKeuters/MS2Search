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
      
      res <- results_r()
      if(is.null(res)){
        filters <- NULL
      }else{
        start_filtercols <- ifelse("metaboScape_ID" %in% colnames(res), 3, 2)
        filters <- colnames(res)[start_filtercols:(ncol(res)-1)]
      }
      

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
              selectizeInput(
                inputId = ns("filter_cols"),
                label = "Filters",
                choices  = filters,
                selected = NULL,
                multiple = TRUE,                 
                options  = list(
                  placeholder = "No filters selected",
                  plugins     = list("remove_button")
                )
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
    
    circle_idx <- reactiveVal(1L)
    
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
      if(is.null(Feature_keys())) keys <- as.character(1:nrow(embeddings())) else keys <- Feature_keys()
      embed_df <- data.frame(id = keys, x = coords[, 1], y = coords[, 2])
      
      res <- results_r()
      if (!is.null(res) && nrow(res) > 0) {
        start_filtercols <- ifelse("metaboScape_ID" %in% colnames(res), 3, 2)
        filter_res <- res[, start_filtercols:(ncol(res) - 1), drop = FALSE]
        
        is_match <- function(x) {
          if (is.list(x)) {
            vapply(x, function(el) length(el) > 0 && !all(is.na(el)), logical(1))
          } else {
            !is.na(x) & x
          }
        }
        
        flags <- as.data.frame(lapply(filter_res, is_match), check.names = FALSE)
        names(flags) <- paste0("f__", names(filter_res))
        flags$id <- res$FEATURE_ID
        
        embed_df <- dplyr::left_join(
          embed_df, dplyr::select(res, id = FEATURE_ID, score), by = "id"
        )
        embed_df <- dplyr::left_join(embed_df, flags, by = "id")
      }
      embed_df
    })
    
      
    output$umap <- plotly::renderPlotly({
      req(plot_df())
      df <- plot_df()

      base_marker <- list(size = 9, line = list(color = "white", width = 1))
      
      # selected filters that exist as flag columns
      fcols <- character(0)
      if (input$color_mode == "filter") {
        sel_filters <- input$filter_cols

        sel_filters <- sel_filters[nzchar(sel_filters)] 
        fcols <- intersect(paste0("f__", sel_filters), names(df))
      }
      
      # hover text
      df$hover <- paste0("<b>", df$id, "</b>")
      if ("score" %in% names(df)) {
        df$hover <- paste0(df$hover, "<br>Score: ", round(as.numeric(df$score), 3))
      }
      if (length(fcols) > 0) {
        F <- as.matrix(df[, fcols, drop = FALSE])
        F[is.na(F)] <- FALSE
        df$hover <- paste0(
          df$hover, "<br>Filters: ",
          vapply(seq_len(nrow(F)), function(i) {
            m <- sub("^f__", "", fcols[F[i, ]])
            if (length(m)) paste(m, collapse=", ") else "-"
          }, character(1))
        )
      }
      
      p <- if (input$color_mode == "score") {
        plotly::plot_ly(
          df, x = ~x, y = ~y, key = ~id, text = ~hover, hoverinfo = "text",
          source = "umap", type = "scattergl", mode = "markers",
          marker = c(base_marker, list(
            color = ~score, colorscale = "Viridis", opacity = 0.85,
            colorbar = list(title = "Score", thickness = 10, len = 0.5, outlinewidth = 0)
          ))
        )
      } else {
        # grey base trace (carries hover + click)
        plotly::plot_ly(
          df, x = ~x, y = ~y, key = ~id, text = ~hover, hoverinfo = "text",
          source = "umap", type = "scattergl", mode = "markers", showlegend = FALSE,
          marker = c(base_marker, list(color = "#c5cad6", opacity = 0.8))
        )
      }
      
      n_traces <- 1L
      
      if (length(fcols) > 0) {
        pal <- stats::setNames(scales::hue_pal()(length(fcols)), fcols)
        rank <- matrix(0L, nrow = nrow(F), ncol = ncol(F))
        for (j in seq_len(ncol(F))) {
          rank[, j] <- if (j == 1) F[, 1] else rowSums(F[, 1:j, drop = FALSE])
        }
        
        rank <- rank * F
        m <- rowSums(F)    # number of filters per dot
        
        if (length(fcols) == 1) { rank <- matrix(rank, ncol = 1); }
        
        for (i in seq_along(fcols)) {
          rows <- which(F[, i])
          if (length(rows) == 0) next
          sizes <- 9 + 5 * (m[rows] - rank[rows, i])   # first filter = outer ring
          p <- plotly::add_trace(
            p, data = df[rows, ], x = ~x, y = ~y, key = NULL, inherit = FALSE,
            type = "scattergl", mode = "markers",
            name = sub("^f__", "", fcols[i]), showlegend = TRUE, hoverinfo = "skip",
            marker = list(size = sizes, color = pal[[fcols[i]]],
                          line = list(color = "white", width = 1))
          )
          n_traces <- n_traces + 1L
        }
      }
      
      # selection circle = always the last trace
      sel <- isolate(selected_feature_id_r())
      idx <- match(as.character(sel), as.character(df$id))
      has_sel <- !is.na(idx)
      p <- plotly::add_trace(
        p,
        x = if (has_sel) df$x[idx] else df$x[1],
        y = if (has_sel) df$y[idx] else df$y[1],
        key = NULL, inherit = FALSE, type = "scattergl", mode = "markers",
        visible = has_sel, showlegend = FALSE, hoverinfo = "skip",
        marker = list(size = 30, color = "rgba(255,99,71,0.18)",
                      line = list(color = "#e8553d", width = 2))
      )
      circle_idx(n_traces)        
      
      umap_layout(p) |>
        plotly::layout(legend = list(orientation = "h", y = -0.05, x = 0.5, xanchor = "center"))
    })
   
    
    # click -> set selection
    observeEvent(has_umap(), {
      req(has_umap())
      observeEvent(plotly::event_data("plotly_click", source = "umap"), {
        ed <- plotly::event_data("plotly_click", source = "umap")
        req(ed$key, Feature_keys())
        selected_feature_id_r(ed$key[[1]])
      })
      
    }, once = TRUE)

    # selection changed (click or main server) -> only move the circle
    observe({
      sel <- selected_feature_id_r()
      df  <- plot_df()
      ci  <- circle_idx()
      req(df)
      
      proxy <- plotly::plotlyProxy("umap", session)
      idx   <- match(as.character(sel), as.character(df$id))
      
      if (is.na(idx)) {
        plotly::plotlyProxyInvoke(proxy, "restyle", list(visible = FALSE), list(ci))
      } else {
        plotly::plotlyProxyInvoke(
          proxy, "restyle",
          list(x = list(list(df$x[idx])), y = list(list(df$y[idx])), visible = TRUE),
          list(ci)
        )
      }
    })
  
    
    
  })
}


umap_layout <- function(p) {
  axis_style <- list(
    title = list(font = list(size = 12, color = "#8a8f98")),
    showgrid = FALSE, zeroline = FALSE, showline = FALSE,
    showticklabels = FALSE, ticks = ""
  )
  p |>
    plotly::layout(
      xaxis = c(axis_style, list(title = list(text = "UMAP 1", font = axis_style$title$font))),
      yaxis = c(axis_style, list(title = list(text = "UMAP 2", font = axis_style$title$font))),
      paper_bgcolor = "rgba(0,0,0,0)",
      plot_bgcolor  = "white",
      margin = list(l = 40, r = 10, t = 10, b = 40),
      hoverlabel = list(
        bgcolor = "white", bordercolor = "#d0d4dc",
        font = list(size = 14, color = "#2b2f38")
      ),
      uirevision = "keep"
    ) |>
    plotly::config(
      displaylogo = FALSE,
      displayModeBar = "hover",
      modeBarButtonsToRemove = c("select2d", "lasso2d", "autoScale2d",
                                 "hoverClosestCartesian", "hoverCompareCartesian",
                                 "toggleSpikelines")
    ) |>
    plotly::event_register("plotly_click")
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