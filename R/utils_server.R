
#-------------------------------------------------------------------------------
# Save and Load parameters
#-------------------------------------------------------------------------------

#' @title Combining all parameters  
#'
#' @description Combines all parameters in one list.
#' 
#' @param state project parameters
#' @param basicFilterValues mz, rt and count parameters, tolerances and noise filtering.
#' @param toolSettings external tool settings
#' @param cards parameters of Rule cards
#' @param lipid_definition parameters of chains and validations
#' @param headgroup backbone / headgroup of chains
#'
#' @return Returns the list of all parameters
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

#' @title Saving Project  
#'
#' @description Saves the parameters, results and spectra in a zip file.
#' 
#' @details
#' Creates a tmp directionary. Then exports the Spectra object as mgf.gz, 
#' writes a json file from the settings parameters and writes the results 
#' in a csv file while seperating lists with /. At the end the tmp dir is
#' zipped.
#' 
#' @param settings_list settings parameters
#' @param results dataframe of analysis results of the filters and chains.
#' @param spectra Spectra object of the ms2 spectra
#' @param output_zip_path_file path where to save
#'
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
    write.csv(as.data.frame(results), row.names = FALSE, file.path(tmp_dir, "results.csv"))
  }
  


  #zip files
  files_to_zip <- list.files(tmp_dir, full.names = TRUE)
  zip::zip(zipfile = output_zip_path_file, files = files_to_zip, mode = "cherry-pick")
  
  print(paste("File saved at", output_zip_path_file))
  
  unlink(tmp_dir, recursive = TRUE)
}


#' @title Loading Project  
#'
#' @description Loads the parameters, results and spectra from a zip file.
#' 
#' @details
#' Unzips the zip file. The setting parameters are loaded from json, the results
#' from the csv file and the Spectra are loaded from the mgf.gz.
#' 
#' @param zip_path path to zip file
#'
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

#' @title Creating and saving brkrrba-file
#'
#' @description Builds a json format from the chains and then saves it as brkrrba-file.
#' 
#' @param settings_list A list of setting parameters.
#' @param file_path path to save file
#'
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
#' @title Getting the neutral loss locations from the results
#'
#' @description Takes the results, which include the found peaks, and looks 
#' between which peaks the neutral loss is located.
#' 
#' @param results A list of setting parameters.
#' @param feature_id id of the selected spectrum
#' @param rules list of rule parameters
#' @param spectrum A Spectra object with only one spectrum
#' @param tol tolerance in Dalton
#' @param include_precursor Boolean if precusor is added to the spectrum mz values
#' 
#' @return The found peaks with partner
#' 
#' @noRd
neutral_losses_from_results <- function(results, feature_id, rules, spectrum,
                                        tol = 0.01, include_precursor = TRUE) {
  mz  <- unlist(spectrum$mz)
  int <- unlist(spectrum$intensity)
  
  # candidate partner peaks (optionally the precursor as well)
  cand_mz  <- mz
  cand_int <- int
  if (include_precursor) {
    cand_mz  <- c(cand_mz,  spectrum$precursorMz)
    cand_int <- c(cand_int, max(int))
  }
  
  row <- results[results$FEATURE_ID == feature_id, , drop = FALSE]
  colnames(row) <- gsub("\\.", " ", colnames(row))
  if (nrow(row) == 0) return(NULL)

  nm <- gsub('\\.'," ",names(results))
 
  first <- which(nm == "peakCount_filter") + 1
  last  <- which(nm == "score") - 1
  if (length(first) == 0 || length(last) == 0 || first > last) return(NULL)
  filter_cols <- nm[first:last]
  
  rule_names <- gsub('\\.'," ",make.names(vapply(rules, function(r) as.character(r$label)[1], character(1))))
  
  out <- list()
  for (col in filter_cols) {
    cell <- row[[col]][1]

    val  <- paste(unlist(cell), collapse = " ")
    
    if (is.na(val) || !nzchar(val) || val == "NA") next

    r_idx <- which(rule_names == col)

    if (length(r_idx) == 0) r_idx <- which(rule_names == sub("\\.[0-9]+$", "", col))
    if (length(r_idx) == 0) next 
    rule <- rules[[r_idx]]

    label    <- rule$label
    
    found    <- as.numeric(regmatches(val, gregexpr("[0-9]+\\.?[0-9]*", val))[[1]])
    masses <- vapply(rule$entries, function(e) {
      resolve_mass(list(input_type = rule$input_type, value = e$value, is_mzPeak = rule$is_mzPeak),settings_list$project_info$mode)
    }, numeric(1))
    searched <- abs(as.numeric(masses[!is.na(masses)]))

    for (f in found) {
      j <- which.min(abs(mz - f))            # snap to the real found peak
      if (abs(mz[j] - f) > tol) next
      
      for (m in searched) {
        for (target in c(mz[j] - m, mz[j] + m)) {   # partner below or above
          d <- abs(cand_mz - target)
          p <- which.min(d)
          if (d[p] > tol) next
          if (abs(cand_mz[p] - mz[j]) < tol) next    # skip self
          
          out[[length(out) + 1]] <- data.frame(
            mz_found = mz[j],      int_found = int[j],
            mz_partner = cand_mz[p], int_partner = cand_int[p],
            mass = m, delta = abs(cand_mz[p] - mz[j]), label = label
          )
        }
      }
    }
  }
  if (length(out) == 0) return(NULL)
  res <- do.call(rbind, out)
  # remove duplicates (same pair found from both ends)
  key <- paste(pmin(res$mz_found, res$mz_partner),
               pmax(res$mz_found, res$mz_partner), res$mass)
  res[!duplicated(key), , drop = FALSE]
}

#' @title Plotting one spectrum
#'
#' @description Creating a plotly plot with colors for the rule filters and 
#' neutral loss visualization.
#' 
#' @param sp A Spectra object with all spectra.
#' @param feature_id id of the selected spectrum
#' @param settings_list list of all setting parameters.
#' @param results The results of the analysis.
#' 
#' @return The plotly figure.
#' 
#' @import plotly
#' @importFrom grDevices rainbow
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
  if(!is.null(filter_results[["Intensity Filter"]])) filter_priorities <- c(filter_priorities, "low")
  
  category <- resolve_filter_category(mz = mz,
                                      filter_results = filter_results,
                                      filter_priority = filter_priorities,
                                      tol = 1e-4)
  levs <- sort(unique(unlist(lapply(rules, function(x){x$label}))))
  colors <- setNames(rainbow(length(levs)), levs)
  colors["none"] <- "#404040"
  colors["Intensity Filter"] <- "#D3D3D3"
  
  # ---- neutral losses -> shapes + annotations (hidden by default) ----
  nl <- neutral_losses_from_results(results, feature_id, rules, spectrum, tol = 0.1)

  shapes <- list(); annots <- list()

  if(!is.null(nl)){
    
    for (i in seq_len(nrow(nl))) {
      y <- max(nl$int_found[i], nl$int_partner[i])  
      col <- unname(colors[nl$label[i]])            
      if (is.na(col)) col <- "black"
      
      shapes[[i]] <- list(
        type = "line", xref = "x", yref = "y",
        x0 = nl$mz_found[i],   y0 = y,
        x1 = nl$mz_partner[i], y1 = y,
        color = category, colors = colors,
        line = list(color = col, width = 1.5, dash = "dot"),
        visible = FALSE
      )
      annots[[i]] <- list(
        x = (nl$mz_found[i] + nl$mz_partner[i]) / 2,
        y = y, yshift = 8,
        xref = "x", yref = "y",
        text = nl$mass[i],
        showarrow = FALSE, font = list(size = 10, color = col),
        visible = FALSE
      )
    }
    n <- length(shapes)
    set_vis <- function(v) {
      out <- list()
      for (i in seq_len(n)) {
        out[[sprintf("shapes[%d].visible", i - 1)]]      <- v
        out[[sprintf("annotations[%d].visible", i - 1)]] <- v
      }
      out
    }
    
  }
  

  fig <- plot_ly() %>%
    add_segments(
      x = mz, xend = mz, y = 0, yend = intensities,
      color = category, colors = colors,
      text = paste0("m/z: ", mz, "<br>Intensity: ", intensities),
      hoverinfo = "text", showlegend = TRUE
    ) %>%
    add_segments(
      x = spectrum$precursorMz, xend = spectrum$precursorMz,
      y = 0, yend = max(intensities),
      line = list(color = "red", dash = "dash", width = 1),
      name = "Precursor", inherit = FALSE
    ) %>%
    layout(
      hovermode = "x",
      xaxis = list(title = "m/z"),
      yaxis = list(title = "Intensity"),
      shapes = shapes,
      annotations = annots,
      margin = list(b = 110)
    )
  
  if (!is.null(nl)) {
    fig <- fig %>% layout(
      updatemenus = list(list(
        type = "buttons", direction = "right",
        x = 0, y = -0.2, xanchor = "left", yanchor = "top",
        active = -1,
        buttons = list(
          list(
            label  = "Show neutral losses",
            method = "relayout",
            args   = list(c(set_vis(TRUE),  list("updatemenus[0].buttons[0].label" = "Hide neutral losses"))),
            args2  = list(c(set_vis(FALSE), list("updatemenus[0].buttons[0].label" = "Show neutral losses")))
          )
        )
      ))
    )
  }
  fig
}


#' @title Using the Rules on one spectrum
#'
#' @description Results are the same of the analysis results - function 
#' obsolete.
#' 
#' @param spectrum A Spectra object with one spectrum.
#' @param settings_list list of all setting parameters.
#' 
#' @return The results of the analysis.
#' 
#' @noRd
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
      
      mz_values <-find_neutralLosses(
            peaks = peaks_data[,"mz"],
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

#' @title Ordering the filter results by priority
#'
#' @description For each mz peak in the spectrum the results are compared and
#' from the fitting ones the filter name with the highest priority is taken.
#' 
#' @param mz The mz values of the spectrum
#' @param filter_results results of the analysis
#' @param filter_priority priority of the result filters
#' @param tol tolerance in dalton (default 1e-4)
#' 
#' @return The names of the priotized rules for each mz value.
#' 
#' @noRd
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

#' @title Plotting ms1 for one feature
#'
#' @description Selects the row in the ms1 data with the selected feature id and
#' plots it with ggplot2.
#' 
#' @param ms1_data The dataframe of the ms1 data.
#' @param selectedFeature The feauture id of the selected feature.
#' 
#' @return The ggplot plot
#' 
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
  
  ggplot(df, aes(x = .data$Sample, y = .data$Intensity)) +
    geom_bar(stat = "identity",fill = "#337ab7")+
    coord_flip()+
    theme_minimal()
  
}

