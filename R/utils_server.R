
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
create_settings_list <- function(state, basicFilterValues, toolSettings, cards, lipid_definition, headgroup = NULL){
  
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
      min_filtering_type  = basicFilterValues$min_filtering_type,
      fixed_threshold = basicFilterValues$fixed_threshold ,
      percent_trehsold = basicFilterValues$percent_trehsold ,
      noiseFiltering_method = basicFilterValues$noiseFiltering_method,
      tolerance_ppm = basicFilterValues$tolerance_ppm,
      minimum_dalton = basicFilterValues$minimum_dalton
    ),
    chain_info = list(
      headgroup = headgroup
    )
  )
  settings_list <- c(
    settings_list,
    list(specific_filters = cards)
  )
  settings_list <- c(settings_list, list(lipid_definition = lipid_definition))
  
  return(settings_list)
}

#' @noRd
setup_and_save_project <- function(settings_list, results, spectra, output_zip_path_file) {
  print(output_zip_path_file)
  tmp_dir <- tempfile("CeraSearch_save_dir")
  dir.create(tmp_dir)
  
  if(!is.null(spectra)){
    con <- gzfile(file.path(tmp_dir, "spectra.mgf.gz"), "w")
    Spectra::export(spectra, backend = MsBackendMgf(), file = con)
    close(con)    
  }

  jsonlite::write_json(settings_list, file.path(tmp_dir, "settings.json"), auto_unbox = TRUE)
  

  if(!is.null(results)){
    convert_list <- function(x) {
      if (is.list(x)) {
        vals <- sapply(x, function(y) {
          if (length(y) == 0) "" else paste(y, collapse = "/")
        })
      } else {
        x
      }
    }
    results <- lapply(results, convert_list)
    write.csv(as.data.frame(results), file.path(tmp_dir, "results.csv"))
  }
  


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
    
    results[] <- lapply(results, function(col) {
      if (is.character(col)) {
        lapply(col, function(x) {
          if (grepl("/", x, fixed = TRUE)) {
            as.numeric(strsplit(x, "/", fixed = TRUE)[[1]])
          } else {
            x
          }
        })
      } else {
        col
      }
    })
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
plot_spectrum_by_id <- function(spectra, feature_id, id_column = "FeatureID", ...) {
  
  spec_subset <- get_spectrum_by_id(spectra, feature_id, id_column)
  
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
get_spectrum_by_id <- function(spectra, feature_id, id_column = "FeatureID") {
  
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
  
}


#-------------------------------------------------------------------------------
# create brkrrba file
#-------------------------------------------------------------------------------

#' @import dplyr
#' @noRd
create_brkrrba_file <- function(settings_list, file_path){
  
  searchedSpecies <- "CER"
  annotationName <- sprintf("%s {carbonsInAllChains}:{DBEinAllChains}", searchedSpecies )
  
  chain_opt_mods <- unname(lapply(seq_along(settings_list$lipid_definition), function(i) {
    x <- settings_list$lipid_definition[[i]]
    if(length(x$modifications) == 0)return(NULL)
    
    mods <- x$modifications
    formulas <- lapply(mods, function(m){if(m$optional) m$formula})

    counts <- tibble(formula = unlist(formulas)) %>% count(formula)
    lapply(seq_len(nrow(counts)), function(j) {
      list(
        name = mod_name <- gsub( "=", "DB_",paste0(counts$formula[j], "_", i)),
        valueType = "intRange",
        value = c(0, counts$n[j])
      )
    })
    
  }))

  chain_formulas <- lapply(seq_along(settings_list$lipid_definition), function(i) {
    x <- settings_list$lipid_definition[[i]]

    opt_mods <- lapply(x$modifications, function(m){if(m$optional) m$formula})
    opt_mods <- Filter(Negate(is.null), opt_mods)
    opt_mods <- unique(opt_mods)
    
    mod_DB_count <- 0
    req_mods <- lapply(x$modifications, function(m){
      if(!m$optional){
        if(grepl("=",m$formula, fixed = TRUE)) mod_DB_count <<- mod_DB_count + 1
        gsub( "=", "",m$formula)
      }})
    req_mods <- Filter(Negate(is.null), req_mods)
    

    formula <- paste0(
      "C{carbonsPerChain", i, "} + H({carbonsPerChain", i, "}*2) - H({DBEPerChain", i, "}*2)",
      if (length(req_mods) > 0)
        paste0(paste0(" + ", paste(req_mods, collapse = " + ")), " - H",length(req_mods))
      else "",
      if (length(opt_mods) > 0)
        paste0(
          " + ",
          paste(
            sapply(opt_mods, function(m) {
              mod_name <- gsub( "=", "DB_",paste0(m, "_", i))
              elements <- unlist(strsplit(m, "")[[1]])
              elements <- elements[elements != "="]
              paste0(
              paste(
                sapply(elements, function(e) {
                  paste0(e, "({", mod_name, "})")
                }),
                collapse = " + "
              ), "- H({", mod_name, "})")
            }),
            collapse = " + "
          )
        )
      else ""
    )
    
    #removing one H for each double bound in the required modifications
    print(mod_DB_count)
    if(mod_DB_count > 0){
      formula <- paste0(formula, " - H", mod_DB_count)
    }
    
    #removing one H for each double bound in the optional modifications
    matches <- gregexpr("\\{([^}]+)\\}", formula, perl = TRUE)
    vars <- unique(regmatches(formula, matches)[[1]])
    if(sum(grepl("DB_", vars)) != 0){
      formula <- paste0(formula, " - ", paste(lapply(vars[grepl("DB_", vars)], function(x){
        paste0(" - H(",x,")")
      })))
    }

    
    formula
    
  })
  
  chain_descriptors <- lapply(seq_along(chain_formulas),function(i){
    list(
      label = paste0("Chain ", i),
      query = sprintf(
        "QUERY scaninfo(MS2DATA) WHERE MS2PROD=formula(%s):TOLERANCEMZ=0.01",
        chain_formulas[i]
      ),
      signalInfoLabel = paste0("{carbonsPerChain", i, "}:{DBEPerChain", i, "}")
    )
  })
  
  neutralPrecursorFormula <- paste0("{backbone} + ", paste(chain_formulas, collapse = " + "))
  
  #Building a sublevel for Rules
  if(length(settings_list$specific_filters)!= 0) {subLevel <- list(
    level = "Rule Level",
    annotationName = "Unknown Annotation Name",
    descriptors = lapply(settings_list$specific_filters, function(rule){
      
      masses <- vapply(rule$entries, function(e) {
        resolve_mass(list(input_type = rule$input_type, value = e$value, is_mzPeak = rule$is_mzPeak),settings_list$project_info$mode)
      }, numeric(1))
      
      NL_or_peak <- ifelse(rule$is_mzPeak, "MS2PROD", "MS2NL")
      query <- paste0("QUERY scaninfo(MS2DATA) WHERE ",NL_or_peak ,"=(",
                      paste(masses, collapse = " OR "),
                      "):TOLERANCEMZ=0.01")

      list(
        label = rule$label,
        query = query,
        signalInfoLabel = paste0(NL_or_peak, ": ", paste(sapply(rule$entries, `[[`, "value"), collapse = ", "))
      )
    }),
    requiredSignalsForAnnotation = 0
  )}
  else subLevel <- NULL
  
  
  #full structure
  lipid_def <- list(
    variables = c(list(
      list(
        name = "backbone",
        valueType = "formula",
        value = ifelse(is.null(settings_list$chain_info$headgroup),"", settings_list$chain_info$headgroup)
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
    seq_along(settings_list$lipid_definition))),
    do.call(c, chain_opt_mods)
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
              list(
                level = "Species Level",
                annotationName = annotationName,
                descriptors = chain_descriptors,
                subLevel = subLevel
              )
            
          ))
      )
      )
  )
  

  #json_file <- paste0(settings_list$project_info$out_dir,"/lipid_def_",settings_list$project_info$project_name,".brkrrba")
  json_file <- file_path
  write(toJSON(lipid_def, auto_unbox = TRUE, pretty = TRUE), json_file)
  
  txt <- readLines(json_file)
  txt <- gsub('"([^"]+)":', '\\1:', txt)
  writeLines(txt, json_file)
  
}


#-------------------------------------------------------------------------------
# Plotting with plotly
#-------------------------------------------------------------------------------

#' @import plotly
#' @noRd
plottingSpectra <- function(sp, feature_id, settings_list, results){
  rules <- settings_list$specific_filters
  spectrum <- get_spectrum_by_id(sp, feature_id, id_column = "FEATURE_ID") 
  mz <- unlist(spectrum$mz)
  intensities <- unlist(spectrum$intensity)
  filtered_spectrum <- cleaning_ref(spectrum, settings_list$basicFilters)
  mz_filtered <- unlist(filtered_spectrum$mz)

  filter_results <- useRulesOnSpectrum(spectrum, settings_list)
  filter_results[["Intensity Filter"]] <- setdiff(mz, mz_filtered)

  filter_priorities <- unlist(lapply(rules, function(x){x$priority}))
  if(!is.null(filter_results[["Intensity Filter"]]))filter_priorities <- c(filter_priorities, "low")
  
  category <- resolve_filter_category(mz = mz,
                                      filter_results = filter_results,
                                      filter_priority = filter_priorities,
                                      tol = 1e-4)
  levs <- sort(unique(unlist(lapply(rules, function(x){x$label}))))
  colors <- setNames(rainbow(length(levs)), levs)
  colors["none"] <- "#404040"
  colors["Intensity Filter"] <- "#D3D3D3"
  
  library(plotly)
  
  # 1. Initialize an empty plotly base canvas
  fig <- plot_ly() %>%
    # 2. Draw perfectly straight vertical spikes as vector segments
    add_segments(
      x = unlist(spectrum$mz), 
      xend = unlist(spectrum$mz),
      y = 0, 
      yend = unlist(spectrum$intensity),
      color = category,
      colors = colors,
      text = paste0("m/z: ", unlist(spectrum$mz), "<br>Intensity: ", unlist(spectrum$intensity)),
      hoverinfo = "text",
      showlegend = TRUE
    ) %>%
    # 3. Optimize the layout interaction mechanics
    layout(
      hovermode = "x",          # Pinpoints the nearest line on the X axis instantly
      xaxis = list(title = "m/z"),
      yaxis = list(title = "Intensity")
    )%>%
    add_segments(
      x = spectrum$precursorMz, xend = spectrum$precursorMz,
      y = 0, yend = max(unlist(spectrum$intensity)),
      line = list(color = "red", dash = "dash", width = 1),
      name = "Precursor",
      inherit = FALSE
      )
  
  fig <- fig %>% layout(title = paste("Spectum", feature_id),
                        xaxis = list(title = "m/z"),
                        yaxis = list(title = "intensity"))
  fig
      }


useRulesOnSpectrum <- function(spectrum, settings_list){
  rules <- settings_list$specific_filters
  filter_results <- list()
  denoised_spectrum <- cleaning_ref(spectrum,settings_list$basicFilters )
  for(rule in rules){
    
    if(!rule$filter_noise) peaks_data <- peaksData(spectrum)[[1]]
    else peaks_data <- peaksData(denoised_spectrum)[[1]]
        
    masses <- vapply(rule$entries, function(e) {
      resolve_mass(list(input_type = rule$input_type, value = e$value, is_mzPeak = rule$is_mzPeak),settings_list$project_info$mode)
    }, numeric(1))
    masses <- abs(as.numeric(masses[!is.na(masses)]))
    if(!rule$is_mzPeak){
      
      #searching for one or multiple neutral losses
      if(!rule$sequential){
        rule$min_repeat = 1
        rule$max_repeat = 1
      }
      
      mz_values <-find_ladder(
            peaks = peaks_data,
            precursor = spectrum$precursorMz[[1]],
            values = masses,
            losses_min = rule$min_repeat,
            losses_max = rule$max_repeat,
            rule_min_mz = rule$rule_mz_min,
            rule_max_mz = rule$rule_mz_max,
            tol_ppm = settings_list$basicFilters$tolerance_ppm,
            tol_mz = settings_list$basicFilters$minimum_dalton,
            just_after_precursor = rule$just_after_precursor
          )
      

    }else{
      #Searching for specific mz values in the ms2 data
      retained_value <- findMZpeak(peaksData(spectrum)[[1]],
                                    values = masses, 
                                    tol_ppm = settings_list$basicFilters$tolerance_ppm, 
                                    tol_mz = settings_list$basicFilters$minimum_dalton)
      if(retained_value == TRUE) mz_values <- masses
      else mz_values <- NULL
      
    }
    filter_results[[rule$label]] <- mz_values
  }
  filter_results
}

resolve_filter_category <- function(mz, filter_results, filter_priority, tol = 1e-4) {
  priority_levels = c( "low","basic", "important", "required")
  priority_tuples <- data.frame(
    name = names(filter_results),
    priority = filter_priority
  )
  
  sapply(mz, function(m) {
    hits <- names(filter_results)[
      vapply(filter_results, function(vals) any(abs(vals - m) < tol), logical(1))
    ]
    if (length(hits) == 0) return("none")
    
    priority_tuples <- priority_tuples[priority_tuples$name %in% hits,]
    priority_tuples$name[
      which.max(match(priority_tuples$priority, priority_levels))
    ]

  })
}

#' @import ggplot2
#' @noRd
plotMS1 <- function(ms1_data, selectedFeature){
  if(! selectedFeature %in% ms1_data$FEATURE_ID ) return()
  
  samples <- as.data.frame(ms1_data[, 12:ncol(ms1_data)])
  rownames(samples) <- ms1_data$FEATURE_ID
  
  intensities <- samples[as.character(selectedFeature),]
  
  df <- data.frame(
    Sample = colnames(intensities),
    Intensity = as.numeric(intensities)
  )
  
  ggplot(df, aes(x = Sample, y = Intensity)) +
    geom_bar(stat = "identity",fill = "#337ab7")+
    coord_flip()+
    theme_minimal()
  
}

