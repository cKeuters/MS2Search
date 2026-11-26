#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @import shinyFiles
# @import shinyWidgets
#' @import splines
#' @noRd
app_ui <- function(request) {
  tagList(
    golem_add_external_resources(),
    shinyjs::useShinyjs(),
    
    fluidPage(
      create_topbar_InProject(),
      
      tabsetPanel(
        id = "page",
        type = "hidden",
        
        # ---- Main app ----
        tabPanelBody(
          "main",
          
          # Welcome screen
          shinyjs::hidden(
            div(id = "welcome_screen",
                tagList(
                  fluidRow(
                    column(6,
                           uiOutput("requiredInputUI"),
                           div(class = "flex-distance",
                               actionButton("start_analysis", "Start Analysis",
                                            class = "process-btn"),
                               div(id = "load_project_wrapper",
                                   fileInput("load_project_file", label = NULL, accept = ".zip",
                                             buttonLabel = tagList(icon("folder-open"), "Load Project"),
                                             placeholder = "")
                               )
                           )),
                    column(6,
                           uiOutput("optionalInputUI"))
                  )
                )
            )
          ),
          
          # Analysis screen
          shinyjs::hidden(
            div(id = "analysis_screen",
                uiOutput("analysis_ui"),
                mod_LipidDefinition_ui("LipidDefinition_1"),
                
                shinyjs::hidden(
                  div(id = "cards_section",
                      uiOutput("rule_coordination"),
                      uiOutput("cards_ui", class = "cards-container"),
                      div(id = "cards_container", class = "cards-container"))
                ),
                
                shinyjs::hidden(
                  div(id = "analysis_footer",
                      fluidRow(
                        column(
                          12,
                          div(
                            style = "display:flex; justify-content:center; gap:10px;",
                            actionButton("run_analysis", "Run Analysis",
                                         class = "btn-primary btn-lg start-btn"),
                            div(
                              tags$span(
                                id = "download_wrapper",
                                style = "display: inline-block;",
                                title = "No data available yet. Please load data first.",
                                shinyjs::disabled(
                                  downloadButton("export_brkrrba", "Export brkrrba",
                                                 class = "btn-primary btn-lg start-btn")
                                )
                              ),
                              tags$head(tags$style(HTML("
                              #download_wrapper .disabled,
                              #download_wrapper [disabled] { pointer-events: none; }
                            ")))
                            )
                          )
                        )
                      ))
                )
            )
          ),
          
          # Results screen
          shinyjs::hidden(
            div(id = "results_screen",
                uiOutput("results_ui"),
                fluidRow(
                  column(
                    width = 6,
                    tabsetPanel(
                      tabPanel(
                        "Feature Table",
                        div(
                          tags$details(
                            tags$summary(
                              class = "validation-summary",
                              icon("filter"), "Filtering Table for Filters"
                            ),
                            checkboxGroupInput("filter_select",
                                               "Select Filters that have to be passed:",
                                               choices = NULL)
                          ),
                          uiOutput("featureTableUI")
                          
                        )
                      ),
                      tabPanel(
                        "Network",
                        mod_DreaMS_Network_ui("mod_DreaMS_Network")
                      )
                    ),
                  ),
                  column(
                    width = 6,
                    div(id = "detailPanel", class = "detail-panel",
                        h4("Spectrum Details"),
                        uiOutput("spectrumDetails"))
                  )
                ))
          )
        ),
        
        # ---- About page ----
        tabPanelBody(
          "about",
          tags$div(
            class = "about-back-wrap",
            actionLink("about_back", "\u2190 Back to app", class = "about-back")
          ),
          mod_about_ui("about_1")
        )
      )
    )
  )
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
	add_resource_path(
		"www",
		app_sys("app/www")
	)

	tags$head(
		favicon(),
		bundle_resources(
			path = app_sys("app/www"),
			app_title = "MS2Search"
		),
		tags$link(rel = "stylesheet", type = "text/css", href = "www/custom.css")
		
		# Add here other external resources
		# for example, you can add shinyalert::useShinyalert()
	)
}
