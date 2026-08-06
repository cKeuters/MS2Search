#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @import shinyFiles
#' @import shinyWidgets
#' @import splines
#' @noRd
app_ui <- function(request) {
	tagList(
		# Leave this function for adding external resources
		golem_add_external_resources(),
		shinyjs::useShinyjs(),
		fluidPage(
		  tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
		    
		    # --- Welcome Screen ---
		    shinyjs::hidden(
		      div(id = "welcome_screen", uiOutput("welcome_ui"))
		    ),
		    
		    # --- Analysis Screen ---
		  shinyjs::hidden(div(id = "analysis_screen",
		        uiOutput("analysis_header_ui"), 
		        shinyjs::hidden(
		          div(id = "cards_section",
		              uiOutput("rule_coordination"),
		              uiOutput("cards_ui", class = "cards-container")
		        ),
		        shinyjs::hidden(
		          div(id = "analysis_footer",
		              fluidRow(column(12, align = "center",
		                              actionButton("run_analysis", "Run Analysis", class = "btn-primary btn-lg start-btn")))
		          )
		       ) )
		      )
		    ),
		  # --- Results Screen ---
		  shinyjs::hidden(
		    div(id = "results_screen", uiOutput("results_ui"))
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
			app_title = "MS2SearchGolem"
		)
		# Add here other external resources
		# for example, you can add shinyalert::useShinyalert()
	)
}
