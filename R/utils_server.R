
#-------------------------------------------------------------------------------
# Save and Load parameters
#-------------------------------------------------------------------------------

#' server 
#'
#' @description A utils function
#'
#' @return The return value, if any, from executing the utility.
#'
#' @noRd
create_settings_list <- function(state, basicFilterValues, toolSettings, cards, lipid_definition){
  
  settings_list <- list(
    project_info = list(
      project_name = state$project_name,
      date = as.character(Sys.Date()),
      mode = state$mode,
      out_dir = state$out_dir,
      screen = state$screen
    ),
    LipidIN_settings = list(
      active = toolSettings$LipidIN,
      priority = toolSettings$LipidIN_priority,
      ppm1 = toolSettings$LipidIN_ppm1,
      ppm2 = toolSettings$LipidIN_ppm2
    ),
    DreaMS_settings = list(
      active = toolSettings$DreaMS,
      priority = toolSettings$DreaMS_priority
    ),
    basicFilters = list(
      mz_min = basicFilterValues$mz_min ,
      mz_max = basicFilterValues$mz_max ,
      mz_priority = basicFilterValues$mz_priority,
      rt_min = basicFilterValues$rt_min ,
      rt_max = basicFilterValues$rt_max ,
      rt_priority = basicFilterValues$rt_priority,
      n_peaks_min = basicFilterValues$n_peaks_min ,
      n_peaks_max = basicFilterValues$n_peaks_max,
      n_peaks_priority = basicFilterValues$n_peaks_priority,
      tolerance_ppm = basicFilterValues$tolerance_ppm,
      minimum_dalton = basicFilterValues$minimum_dalton
    )
  )
  settings_list <- c(
    settings_list,
    list(specific_filters = cards)
  )
  settings_list <- c(settings_list, list(lipid_definition = lipid_definition))
  
  return(settings_list)
}

#' @import arrow
#' @noRd
setup_and_save_project <- function(settings_list, results, spectra, output_zip_path_file) {
  
  tmp_dir <- tempfile("CeraSearch_save_dir")
  dir.create(tmp_dir)
  
  if(!is.null(spectra)){
    con <- gzfile(file.path(tmp_dir, "spectra.mgf.gz"), "w")
    export(spectra, backend = MsBackendMgf(), file = con)
    close(con)    
    #write_parquet(spectra, file.path(tmp_dir, "spectra.parquet"))
  }

  jsonlite::write_json(settings_list, file.path(tmp_dir, "settings.json"), auto_unbox = TRUE)
  
  if(!is.null(results))write.csv(as.data.frame(results), file.path(tmp_dir, "results.csv"))

  #zip files
  files_to_zip <- list.files(tmp_dir, full.names = TRUE)
  zip::zip(zipfile = output_zip_path_file, files = files_to_zip, mode = "cherry-pick")
  
  print(paste("File saved at", output_zip_path_file))
  
  unlink(tmp_dir, recursive = TRUE)
}

#' @noRd
load_project <- function(zip_path) {
  tmp_dir <- tempfile("CeraSearch_load_dir")
  dir.create(tmp_dir)
  if(!file.exists(zip_path))stop("Error: Project file does not exist")
  zip::unzip(zip_path, exdir = tmp_dir)
  
  settings_list <- jsonlite::fromJSON(
    file.path(tmp_dir, "settings.json"),
    simplifyVector = FALSE  
  )
  
  results <- NULL
  if (("results.csv" %in% list.files(tmp_dir))) {
    results <- read.csv(file.path(tmp_dir, "results.csv"))
  }
  spectra <- NULL
  if (("spectra.mgf.gz" %in% list.files(tmp_dir))) {
    spectra <- Spectra::Spectra(file.path(tmp_dir, "spectra.mgf.gz"), source = MsBackendMgf::MsBackendMgf())
  }
  message("File loaded")
  unlink(tmp_dir, recursive = TRUE)
  return(list(settings_list = settings_list, results = results, spectra = spectra))
}


#' Plot a spectrum selected by feature ID
#'
#' @param spectra A `Spectra` object.
#' @param feature_id Character or numeric. The ID to match against.
#' @param id_column Character. spectraVariable holding the feature ID.
#' @param ... Additional arguments passed to `Spectra::plotSpectra()`
#'   (e.g. `main`, `xlim`, `ylim`).
#'
#' @return Invisibly returns the subset Spectra object (useful if you
#'   want to reuse it, e.g. for peak annotation on top of the plot).
#'
#' @importFrom Spectra plotSpectra
plot_spectrum_by_id <- function(spectra, noise, feature_id, id_column = "FeatureID", ...) {
  
  spec_subset <- get_spectrum_by_id(spectra, feature_id, id_column, noise)
  
  Spectra::plotSpectra(
    spec_subset,
    main = sprintf("%s: %s", id_column, feature_id),
    ...
  )
  
  invisible(spec_subset)
}


#' Get a single spectrum from a Spectra object by feature ID
#'
#' @param spectra A `Spectra` object containing multiple spectra.
#' @param feature_id Character or numeric. The ID to match against.
#' @param id_column Character. Name of the spectraVariable that holds
#'   the feature ID (e.g. "spectrumId", "feature_id", "acquisitionNum").
#'   Defaults to "spectrumId".
#'
#' @return A `Spectra` object subset to the matching spectrum/spectra.
#'   Throws an informative error if no match is found.
#'
#' @importFrom Spectra spectraData
get_spectrum_by_id <- function(spectra, feature_id, id_column = "FeatureID", noise) {
  
  spectraData <- Spectra::spectraData(spectra)
  
  if (!id_column %in% colnames(spectraData)) {
    stop(sprintf(
      "Column '%s' not found in spectraData. Available columns: %s",
      id_column, paste(colnames(spectraData), collapse = ", ")
    ))
  }
  
  match_idx <- which(spectraData[[id_column]] == feature_id)
  
  if (length(match_idx) == 0) {
    stop(sprintf("No spectrum found for %s = '%s'", id_column, feature_id))
  }
  
  subsp <- spectra[match_idx]
  
  # if(noise){
  #   peaks <- unlist(peaksData(subsp))
  #   peak_noise_vector <- compute_Noise(peaks)
  #   subsp <- addProcessing(subsp, compute_Noise)
  # }
  # print(unlist(peaksData(subsp)))
  if(noise)subsp <- cleaning_ref(subsp)
  subsp
  
}


#-------------------------------------------------------------------------------
# create brkrrba file
#-------------------------------------------------------------------------------

create_brkrrba_file <- function(settings_list){
  
  neutralPrecursorFormula <- "{backbone}+C{carbonsInAllChains}+H({carbonsInAllChains}*2-{DBEinAllChains}*2)"
  
  searchedSpecies <- "CER"
  annotationName <- sprintf("%s {carbonsInAllChains}:{DBEinAllChains}", searchedSpecies )
  
  
  #chains to query
  if(is.null(settings_list$lipid_definition) || length(settings_list$lipid_definition) == 0) sublevel_chain_definition <- list()
  else{
  sublevel_chain_definition <- setNames(
    lapply(seq_along(settings_list$lipid_definition), function(i) {
      x <- settings_list$lipid_definition[[i]]
      
    formula = paste0(
      "C{carbonsPerChain", i, "}+H({carbonsPerChain", i,"}*2-{DBEPerChain", i, "}*2)",
      if (length(x$modifications) > 0)
        paste0(" + ", paste(sapply(x$modifications, `[[`, "formula"), collapse = " + "))
      else ""
      )
    
      list(
        level = "Chain Level",
        annotationName = "Unknown Annotation Name",
        descriptors = list(list(
          label = paste0("Chain ", i),
          query = sprintf(
            "QUERY scaninfo(MS2DATA) WHERE MS2NL=formula(%s):TOLERANCEMZ=0.01",
            formula
          ),
          signalInfoLabel = paste0("{carbonsPerChain", i, "}:{DBEPerChain", i, "}")
          #requiredSignalsForAnnotation = 1
          #validationRule = ""
        ))
      )
      
    
    }),
    paste0("Chain", seq_along(settings_list$lipid_definition))
  )}
  
  if(is.null(settings_list$specific_filters) || length(settings_list$specific_filters) == 0) sublevel_rule_definition <- list()
  else{
  
  sublevel_rule_definition <- setNames(
    lapply(seq_along(settings_list$specific_filters), function(i) {
      x <- settings_list$specific_filters[[i]]
      masses <- vapply(x$entries, function(e) {
        resolve_mass(list(input_type = x$input_type, value = e$value, is_mzPeak = x$is_mzPeak),settings_list$project_info$mode)
      }, numeric(1))
      
      NL_or_peak <- ifelse(x$is_mzPeak, "MS2PROD", "MS2NL")
      query <- paste0("QUERY scaninfo(MS2DATA) WHERE ",NL_or_peak ,"=(",
                      paste(masses, collapse = " OR "),
                      "):TOLERANCEMZ=0.01")
      
      signalInfoLabel <- paste0(NL_or_peak, ": ", paste(sapply(x$entries, `[[`, "value"), collapse = ", "))
      
      list(
        level = "Rule Level",
        annotationName = "Unknown Annotation Name",
        descriptors = list(list(
          label = x$label,
          query = query,
          signalInfoLabel = signalInfoLabel
        )),
        requiredSignalsForAnnotation = 1
      )
    }),
    paste0("Rule", seq_along(settings_list$specific_filters))
  )
  }
  
  #full structure
  lipid_def <- list(
    variables = c(list(
      list(
        name = "backbone",
        valueType = "formula",
        value = settings_list$head_group
      ),
      list(name = "numberOfChains",
           valueType = "int",
           value = length(settings_list$lipid_definition)
           ),
      list(
        name = "carbonsInAllChains",
        valueType = "intRange",
        value = c(sum(unlist(lapply(settings_list$lipid_definition, function(x){x$C_min}))),
                  sum(unlist(lapply(settings_list$lipid_definition, function(x){x$C_max}))))
      ),
      list(
        name = "DBEinAllChains",
        valueType = "intRange",
        value = c(sum(unlist(lapply(settings_list$lipid_definition, function(x){x$DB_min}))),
                  sum(unlist(lapply(settings_list$lipid_definition, function(x){x$DB_max}))))
      )
    ),
    unname(Map(function(x, i) {
      list(
        name = paste0("carbonsPerChain", i),
        valueType = "intRange",
        value = c(x$C_min, x$C_max)
      )
    },
    settings_list$lipid_definition,
    seq_along(settings_list$lipid_definition)))
    ,
    unname(Map(function(x, i) {
      list(
        name = paste0("DBEPerChain", i),
        valueType = "intRange",
        value = c(x$DB_min, x$DB_max)
      )
    },
    settings_list$lipid_definition,
    seq_along(settings_list$lipid_definition)))
    ),
    annotation = 
      list(list(
        polarityType = ifelse(settings_list$project_info$mode == "pos", "positive", "negative"),
        annotationProfile = 
          list(list(
            ionNotation = ifelse(settings_list$project_info$mode == "pos", "[M+H]+", "[M-H]-"),
            neutralPrecursorFormula = neutralPrecursorFormula,
            #regressionParametersCCS = list(),
            annotationConfiguration =
              c(
                level = "Species Level",
                annotationName = annotationName,
                #descriptors = list(),
                sublevel_chain_definition,
                sublevel_rule_definition
              )
          ))
      )
      )
  )
  
  json_file <- paste0(settings_list$project_info$out_dir,"/lipid_def_",settings_list$project_info$project_name,".brkrrba")
  write(toJSON(lipid_def, auto_unbox = TRUE, pretty = TRUE), json_file)
  
  txt <- readLines(json_file)
  txt <- gsub('"([^"]+)":', '\\1:', txt)
  writeLines(txt, json_file)
  
}








