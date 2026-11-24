options(shiny.autoreload = TRUE)

# Set options here
options(golem.app.prod = FALSE) # TRUE = production mode, FALSE = development mode

devtools::install(quick = TRUE, upgrade = FALSE)

# Detach all loaded packages and clean your environment
golem::detach_all_attached()
# rm(list=ls(all.names = TRUE))

# Document and reload your package
golem::document_and_reload()

#Only Reload
#devtools::load_all()
pkgload::load_all()

# Run the application
run_app()
