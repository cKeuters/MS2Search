#' Run the Shiny Application
#'
#' @param ... arguments to pass to golem_opts.
#' See `?golem::get_golem_options` for more details.
#' @inheritParams shiny::shinyApp
#'
#' @export
#' @importFrom shiny shinyApp
#' @importFrom golem with_golem_options
run_app <- function(
	onStart = NULL,
	options = list(),
	enableBookmarking = NULL,
	uiPattern = "/",
	...
) {
  options(shiny.maxRequestSize = 500 * 1024^2)
  #seed for chain id generation
  set.seed(Sys.time())
  
  #cores for Parallelity 
  n_cores   <- parallelly::availableCores()
  print(paste("Available cores:", n_cores))
  n_workers <- max(2, min(n_cores - 1, 4))  # leave one core for Shiny, cap at 4
  future::plan(future::multisession, workers = n_workers)
  
	golem::with_golem_options(
		app = shinyApp(
			ui = app_ui,
			server = app_server,
			onStart = onStart,
			options = options,
			enableBookmarking = enableBookmarking,
			uiPattern = uiPattern
		),
		golem_opts = list(...)
	)
}
