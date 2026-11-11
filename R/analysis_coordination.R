
#' @title Loading MS2 data 
#'
#' @description Creating a Spectra Object from a .mgf file
#' 
#' @details
#' Loading in a .mgf file as Spectra object by MsBackendMgf with a file path. 
#' Setting the retention time from possible variables in th .mgf file.
#' 
#' @param ms2file Holds the filepath of a mgf-file.
#' @param add_log A function to write logs.
#'
#' @return Returns a spectra object.
#'
#' @importFrom Spectra Spectra 
#' @importFrom MsBackendMgf MsBackendMgf
#' @noRd
load_spectra <- function(ms2_file, add_log = NULL){
  
  #if(is.null(ms2_file))ms2_file <- "path/data/Test_dataset.mgf"
  print(ms2_file)
  
  sp <- Spectra::Spectra(ms2_file, source = MsBackendMgf::MsBackendMgf())
  sp$FEATURE_ID <- as.numeric(sp$FEATURE_ID)
  sp <- getRT(sp)
  
  FEATURE_ID <- spectraData(sp)$FEATURE_ID
  if(length(FEATURE_ID) != length(unique(FEATURE_ID))){
    n_peaks <- lengths(peaksData(sp))
    
    key <- paste(
      sprintf("%.5f", spectraData(sp)$precursorMz),
      sprintf("%.2f", spectraData(sp)$rtime),
      n_peaks,
      sep = "|"
    )
    sp$FEATURE_KEY <- key
  }else sp$FEATURE_KEY <- FEATURE_ID
  
  if(!is.null(add_log)) add_log(paste("Loaded MS2-file with ", length(sp)," spectra."))

  sp
}


#' @title Get the MetaboScape Feature Table as dataframe 
#'
#' @description Gets the feature table from the MetaboScape API and combines all
#' values in a dataframe ready to create a Spectra Object.
#' 
#' @details
#' Gets the feature table from the API and selects it in the project to find out its 
#' polarity. MS-level is set to 2. Retention time, precursor mz, Feature_ID, 
#' metaboScape feature ID, adducts and polarity values are combined in a dataframe.
#' For each feature the ms2 values are pulled from the API and mz and intensity 
#' of the ms2 spectra is added to the dataframe.
#' 
#' @param operations The MetaboScape API functions from rapiclient.
#' @param featureTableId MetaboScape ID of the chosen feature table.
#'
#' @return Returns a dataframe containing the spectra information.
#'
#' @importFrom S4Vectors DataFrame
#' @importFrom future.apply future_lapply
#' @importFrom purrr map list_rbind
#' @import dplyr
#' @noRd
createSpectradf <- function(operations, featureTableId){
  res <- operations$retrieveFeatureTable(featureTableId = featureTableId)
  featureTable <- fromJSON(httr::content(res, as = "text", encoding = "UTF-8"))
  
  res <- operations$listAllProjects()
  projects <- fromJSON(httr::content(res, as = "text", encoding = "UTF-8"))

  all_feature_tables <- projects$experiments |>
    map(\(exp) {
      map(exp$featureTables, identity) |> 
        list_rbind()
    }) |>
    list_rbind()
  
  polarity <- all_feature_tables |>
    filter(id == featureTableId) |>
    pull(polarity)
  if(polarity == "POSITIVE") pol <- 1 else pol <- 0

  adducts <- vapply(
    featureTable$featureIons,
    function(x) paste(x$ionNotation, collapse = ";"),
    character(1)
  )

  # plan(multisession, workers = 4)
  msms <- future_lapply(featureTable$id, function(id){
    res <- operations$retrieveFeatureMsMsSpectra(id)
    fromJSON(httr::content(res, as = "text", encoding = "UTF-8"))
  })
  mz <- lapply(msms, function(feature){feature$spectra[[1]]$signals[[1]]$mz})
  intensity <- lapply(msms, function(feature){feature$spectra[[1]]$signals[[1]]$intensity})
  

  spectra_df <- S4Vectors::DataFrame(
    msLevel = as.integer(2),
    rtime = as.numeric(featureTable$rtInSeconds), 
    precursorMz = as.numeric(featureTable$mass),  
    FEATURE_ID = as.character(rownames(featureTable)),
    FEATURE_KEY = as.character(rownames(featureTable)),
    scanIndex  = as.integer(rownames(featureTable)),
    metaboScape_ID = featureTable$id,
    adducts = adducts,
    polarity = as.integer(pol),
    mz = IRanges::NumericList(mz, compress = FALSE),
    intensity = IRanges::NumericList(intensity, compress = FALSE)
  )

  spectra_df
}



#' @title Coodinates the analysis 
#'
#' @description Uses the filters, chains and rules to process the ms2 data.
#' 
#' @details
#' First there are the basic filters (ms-level, polarity, mz, rt, peak count).
#' Then the rules are applied. The highest peak from a neutral loss is saved. 
#' For the chains peaks near the approximated masses are searched. For the 
#' validation the masses of the chains are combined and searched.
#' All filter with the priority "required" remove spectra which did not pass.
#' The score is calculated as 1 for basic and 2 for important priority if the 
#' filter got passed for that spectrum. It is divided by the sum of all priorities. 
#' 
#' @param settings_list A list of lists (json) holding all information of the 
#' analysis settings.
#' @param sp A Spectra object containing the ms2 spectra information.
#' @param add_log A function to write logs.
#'
#' @return Returns a dataframe with the results of each filter, chain and rule. 
#'
#' @importFrom Spectra filterMsLevel filterEmptySpectra spectraVariables spectraData rtime peaksData precursorMz
#'  
#' @importFrom MetaboCoreUtils calculateMass 
#' @noRd
run_analysis <- function(settings_list, sp, add_log){
  #settings_list <- load_project("path/data/28.9.26_meeting2.zip")$settings_list

  #-----------------------------------------------------------------------------
  # Basic filters and checks
  #-----------------------------------------------------------------------------

  #check MS-Level 2 & empty spectra
  keep <- !vapply(peaksData(sp), is.null, logical(1))
  sp <- sp[keep]
  sp <- sp |>
    filterMsLevel(2) |>
    filterEmptySpectra() 
  add_log(paste("Non-empty spectra with MS-level 2: ", length(sp)))
  
  results <- data.frame(FEATURE_ID = spectraData(sp)$FEATURE_ID)
  if(length(results$FEATURE_ID) != length(unique(results$FEATURE_ID))){
    results$FEATURE_ID <- spectraData(sp)$FEATURE_KEY
  }
  
  if("metaboScape_ID" %in% spectraVariables(sp)) results$metaboScape_ID <- spectraData(sp)$metaboScape_ID

  #check mode
  pol_col <- intersect(c("POLARITY", "IONMODE", "ionmode", "polarity"),spectraVariables(sp))[1]
  if(settings_list$project_info$mode == "pos") mode <- c("POS","Positive", "POSITIVE", "positive", 1) else mode <- c("NEG", "NEGATIVE", "negative", 0)

  if (length(pol_col) > 0) {
    retained_vector <- sp[[pol_col[1]]] %in% mode
    results <- results[retained_vector, , drop = FALSE]
    sp <- sp[retained_vector,]
    }
  add_log(paste0("Spectra in positive mode: ", length(sp)))
  
  #check mz
    all_mz <- precursorMz(sp)
    retained_vector <- all_mz >= settings_list$basicFilters$mz_min & all_mz <= settings_list$basicFilters$mz_max
    retained_vector[is.na(retained_vector)] <- FALSE
    add_log(paste0("Spectra in mz-Range (", settings_list$basicFilters$mz_min,", ", settings_list$basicFilters$mz_max, "): ", sum(retained_vector)))
    
    results$mz_filter <- retained_vector
    if(settings_list$basicFilters$mz_priority == "required"){
      results <- results[retained_vector,]
      sp <- sp[retained_vector,]
    }

  #check RT
    all_RT <- rtime(sp)
    retained_vector <- all_RT >= settings_list$basicFilters$rt_min & all_RT <= settings_list$basicFilters$rt_max
    retained_vector[is.na(retained_vector)] <- FALSE
    add_log(paste0("Spectra in RT-Range (", settings_list$basicFilters$rt_min,", ", settings_list$basicFilters$rt_max, "): ", sum(retained_vector)))
    results$rt_filter <- retained_vector
    if(settings_list$basicFilters$rt_priority == "required"){
      results <- results[retained_vector,]
      sp <- sp[retained_vector,]
    }
    
  sp_cleaned <- cleaning_ref(sp, settings_list$basicFilters)
  
    
  #check peak count
    all_PeakCount <- lengths(sp_cleaned)
    retained_vector <- all_PeakCount >= settings_list$basicFilters$n_peaks_min & all_PeakCount <= settings_list$basicFilters$n_peaks_max
    retained_vector[is.na(retained_vector)] <- FALSE
    add_log(paste0("Spectra in peak-count-range (", settings_list$basicFilters$n_peaks_min,", ", settings_list$basicFilters$n_peaks_max, "): ", sum(retained_vector)))
    results$peakCount_filter <- retained_vector
    if(settings_list$basicFilters$n_peaks_priority == "required"){
      results <- results[retained_vector,]
      sp <- sp[retained_vector,]
      sp_cleaned <- sp_cleaned[retained_vector,]
    }
    
    
    
  #-----------------------------------------------------------------------------
  # Rule filters
  #-----------------------------------------------------------------------------

  for(rule in settings_list$specific_filters){
    
    if(!rule$filter_noise)  {peaks_data <- peaksData(sp)
    }else  peaks_data <- peaksData(sp_cleaned)
    
    masses <- vapply(rule$entries, function(e) {
      resolve_mass(list(input_type = rule$input_type, value = e$value, is_mzPeak = rule$is_mzPeak),settings_list$project_info$mode)
    }, numeric(1))
    masses <- abs(as.numeric(masses[!is.na(masses)]))

    if(!rule$is_mzPeak){
      
      #searching for multiple neutral losses
      if(!rule$sequential){
        rule$min_repeat = 1
        rule$max_repeat = 1
      }
      ladders <- Map(
        function(peaks, precursor) {
          find_neutralLosses(
            peaks = peaks,
            precursor = precursor,
            values = masses,
            losses_min = rule$min_repeat,
            losses_max = rule$max_repeat,
            rule_min_mz = rule$rule_mz_min,
            rule_max_mz = rule$rule_mz_max,
            tol_ppm = settings_list$basicFilters$tolerance_ppm,
            tol_mz = settings_list$basicFilters$minimum_dalton,
            just_after_precursor = rule$just_after_precursor
          )
        },
        peaks_data[,"mz"],
        sp$precursorMz
      )

      retained_vector <- lengths(ladders) > 0
      add_log(paste0("Spectra after ",rule$label,": ", sum(retained_vector)))
      results[[rule$label]] <- ladders
      if(rule$rule_priority == "required"){
        results <- results[retained_vector,]
        sp <- sp[retained_vector,]
        sp_cleaned <- sp_cleaned[retained_vector,]
        
      }
      
    }else{
      #Searching for specific mz values in the ms2 data
      found_mz <- lapply(peaks_data, findMZpeak,
                values = masses, 
                tol_ppm = settings_list$basicFilters$tolerance_ppm, 
                tol_mz = settings_list$basicFilters$minimum_dalton)
      retained_vector <- lengths(found_mz) > 0
      add_log(paste0("Spectra after ",rule$label,": ", sum(retained_vector)))
      results[[rule$label]] <- found_mz
      if(rule$rule_priority == "required"){
        
        results <- results[retained_vector,]
        sp <- sp[retained_vector,]
        sp_cleaned <- sp_cleaned[retained_vector,]
      }

    }

  }
    
    #-----------------------------------------------------------------------------
    # Chain filters
    #-----------------------------------------------------------------------------
    peaks_data <- peaksData(sp_cleaned)

    chain_id <- 1
    for(ld in settings_list$lipid_definition$chains){
      masses <- calc_lipid_chain_masses(ld)$mass
      
      found_mz <- lapply(peaks_data, findMZpeak,
                                       values = masses, 
                                       tol_ppm = settings_list$basicFilters$tolerance_ppm, 
                                       tol_mz = settings_list$basicFilters$minimum_dalton)
      
      retained_vector <- !is.na(found_mz)
      add_log(paste0("Spectra after Chain ",chain_id,": ", sum(retained_vector)))
      
      results[[paste("chain",chain_id)]] <- found_mz
      
      chain_id <- chain_id + 1

      if(ld$priority == "required"){
        results <- results[retained_vector,]
        sp <- sp[retained_vector,]
        sp_cleaned <- sp_cleaned[retained_vector,]
        peaks_data <- peaksData(sp_cleaned)
      }
      
    }
    
    #-----------------------------------------------------------------------------
    # validation filters
    #-----------------------------------------------------------------------------
    
    if(length(settings_list$lipid_definition$chains) == 0){
      settings_list$lipid_definition$validations <- NULL
    }
    else{
      headgroupMass <- if(is.null(settings_list$lipid_definition$head_group)||settings_list$lipid_definition$head_group == "") 0 
      else {calculate_adjustment_mass(settings_list$lipid_definition$head_group)}

      for(val in settings_list$lipid_definition$validations){
        chain_ids <- val$components
        chain_ids <- setdiff(chain_ids, "Backbone")
        
        components <- lapply(chain_ids, function(id){
          chain <- settings_list$lipid_definition$chains[
            sapply(settings_list$lipid_definition$chains, `[[`, "chain_id") == id
          ][[1]]
          calc_lipid_chain_masses(chain)$mass}
          )

        combos <- do.call(expand.grid, components)

        adjustmentMass <- calculate_adjustment_mass(val$adjustments)
        
        masses <- rowSums(combos) + adjustmentMass + if("Backbone" %in% val$components) headgroupMass else 0 
        
        found_mz <- lapply(peaks_data, findMZpeak,
                           values = masses, 
                           tol_ppm = settings_list$basicFilters$tolerance_ppm, 
                           tol_mz = settings_list$basicFilters$minimum_dalton)
        retained_vector <- !is.na(found_mz)
        
        add_log(paste0("Spectra after ",val$name,": ", sum(retained_vector)))
        
        results[[val$name]] <- found_mz
        
        if(val$val_priority == "required"){
          
          results <- results[retained_vector,]
          sp <- sp[retained_vector,]
          sp_cleaned <- sp_cleaned[retained_vector,]
        }
      }
    }

  #------------------------- creating dataframe  ----------------------
    priority_weights <- c(
      basic = 1,
      important = 2,
      required = 0
    )
    
    priorities <- c(
      priority_weights[settings_list$basicFilters$mz_priority],
      priority_weights[settings_list$basicFilters$rt_priority],
      priority_weights[settings_list$basicFilters$n_peaks_priority],
      sapply(settings_list$specific_filters, \(x) priority_weights[x$rule_priority]),
      sapply(settings_list$lipid_definition$chains, \(x) priority_weights[x$priority]),
      sapply(settings_list$lipid_definition$validations, \(x) priority_weights[x$val_priority])
    )
    priorities <- unname(sapply(priorities, function(x) if (length(x) == 0) 0 else x))
    
    start_filtercols <- ifelse("metaboScape_ID" %in% colnames(results), 3, 2)
    matches <- results[, start_filtercols:ncol(results)]
    
    is_match <- function(x) {
      if (is.list(x)) {
        vapply(x, function(el) {
          length(el) > 0 && !all(is.na(el))
        }, logical(1))
      } else {
        !is.na(x) & x
      }
    }
    
    matches <- lapply(matches, is_match)
    matches <- as.data.frame(matches)
    matches[] <- lapply(matches, as.integer)
    matches <- as.matrix(matches)
    
    results$score <- rowSums((sweep(matches, MARGIN=2,priorities, `*`))) / sum(priorities)
  
    as.data.frame(results)
}



#' @title Adds Retention time to a Spectra 
#'
#' @description Converts various writings of retention time to rtime in a Spectra
#' Object.
#' 
#' @details
#' Takes RTINSECONDS, RETENTIONTIME, RETENTION_TIME, RT, RTINMINUTES from the
#' spectraVariables and adds the variable rtime in seconds to the Spectra Object.
#' 
#' @param sp A Spectra object containing the ms2 spectra information.
#'
#' @return Returns a Spectra object containing the ms2 spectra information including
#' retentiontime as rtime.
#' 
#' @importFrom Spectra spectraVariables
#'
#' @noRd
getRT <- function(sp){
  if (all(is.na(rtime(sp)))) {
    
    vars <- spectraVariables(sp)
    
    rt_map <- c(
      RTINSECONDS  = "seconds",
      RETENTIONTIME = "seconds",
      RETENTION_TIME = "seconds",
      RT = "seconds",
      RTINMINUTES = "minutes"
    )
    
    match_idx <- match(names(rt_map), toupper(vars))
    
    if (any(!is.na(match_idx))) {
      
      field <- vars[match_idx[which(!is.na(match_idx))[1]]]
      
      rt_values <- as.numeric(sp[[field]])
      
      if (rt_map[[toupper(field)]] == "minutes")
        rt_values <- rt_values * 60
      
      sp$rtime <- rt_values
    }
  }
  sp
}

#' @title Calculates the mass of the adjustment or backbone string.
#' 
#' @param adjustment A string of elements and "-".
#'
#' @return Returns the mass of the elements.
#' 
#' @importFrom MetaboCoreUtils calculateMass
#'
#' @noRd
calculate_adjustment_mass <- function(adjustment){
  
  adjustmentMass <- 0
  
  if(!is.null(adjustment) && adjustment != ""){
    adjustment <- str_trim(unlist(strsplit(split = "-", adjustment)))
    adjustment_masses <- unlist(lapply(adjustment, function(adjustment)
      {if(adjustment != "")MetaboCoreUtils::calculateMass(adjustment) else 0}))
    adjustmentMass <- ifelse(is.na(adjustment_masses[1]), 0, adjustment_masses[1]) - sum(adjustment_masses[2:length(adjustment_masses)])
  }
  
  adjustmentMass
}

#' @title Finding neutral loss pattern 
#'
#' @description Searching one spectrum for specific neutral losses between peaks.
#' 
#' @details
#' A possible precursor is added to the mz values of the spectrum (peaks). After
#' that mz values which are not in the min and max mz range are removed.  
#' A dataframe is build with the anchor as index of the mz value and the offset as 
#' one searched value. In a loop the candidates are calculated as mz value plus
#' the offset plus the offset from the previous run. For each combination of 
#' peak and possible offsets a row exists. Each loop run candidates that are not 
#' found are removed from the dataframe and the remaining candidates get another
#' offset. 
#' 
#' 
#' @param peaks mz values of one spectrum.
#' @param precursor mz value of the precursor (default NULL).
#' @param values The values which are to be searched in the spectrum.
#' @param losses_min The minimum number of sequential neutral losses (default 1).
#' @param losses_max The maximum number of sequential neutral losses (default 1).
#' @param tol_ppm Tolerance in ppm (default 10).
#' @param tol_mz Tolerance in Dalton (default 0.1).
#' @param rule_min_mz Minimum mz where to be searched (default NULL). 
#' @param rule_max_mz Maximum mz where to be searched (default NULL). 
#' @param just_after_precursor If only searched after the precursor (default FALSE).
#'
#' @return returns the values of the mz values in range of the tolerance to the 
#' searched values. For sequences the highest value.
#'
#' @importFrom MsCoreUtils closest 
#' @export
find_neutralLosses <- function(peaks, precursor = NULL, values, losses_min = 1, 
                               losses_max = 1, tol_ppm = 10, tol_mz = 0.1, 
                               rule_min_mz = NULL, rule_max_mz = NULL, 
                               just_after_precursor = FALSE) {

  mz <- peaks
  if(!is.null(precursor)) mz <- unique(c(mz, precursor))

  if(!is.na(rule_min_mz) && !is.null(rule_min_mz))mz <- mz[mz > rule_min_mz]
  if(!is.na(rule_max_mz) && !is.null(rule_max_mz))mz <- mz[mz < rule_max_mz]
  mz <- sort(mz[!is.na(mz)])
  
  values <- sort(as.numeric(values))
  n_peaks <- length(mz)
  
  if (n_peaks < losses_min) return(numeric(0))
  
  max_check <- losses_max + 1L
  live <- data.frame(anchor = seq_len(n_peaks), offset = 0)

  if(just_after_precursor) live <- live[live$anchor == which(mz == precursor),]
  run_length <- integer(n_peaks)
  
  for (lvl in seq_len(max_check)) {
    if (nrow(live) == 0) break
    
    #each living chain plus each possible loss
    cand <- merge(live, data.frame(delta = values))
    cand$offset_new <- cand$offset + cand$delta
    cand$target <- mz[cand$anchor] - cand$offset_new
    cand <- cand[order(cand$target), ]

    idx <- closest(cand$target, mz, tolerance = tol_mz, ppm = tol_ppm, duplicates = "closest")
    
    cand <- cand[!is.na(idx), ]
    if (nrow(cand) == 0) break
    
    live <- unique(cand[, c("anchor", "offset_new")])
    names(live) <- c("anchor", "offset")
    
    run_length[unique(live$anchor)] <- lvl
  }
  
  hit <- run_length >= losses_min & run_length <= losses_max

  mz[hit]
}

findMZpeak <- function(peaks, values, tol_ppm, tol_mz){
  mz <- peaks[,"mz"]
  values <- sort(as.numeric(values))
  
  idx <- closest(values, mz, tolerance = tol_mz, ppm = tol_ppm, duplicates = "closest")
  res <- mz[idx]
  if(all(is.na(res))) NA
  else res[!is.na(res)]
}

#' @title Noise filtering of a Spectra object
#'
#' @description Using different Methods to filter the MS2 spectra of a Spectra
#' object.
#' 
#' @details
#' Uses MsCoreUtils::noise with the method MAD or SuperSmoother or MsCoreUtils::smooth
#' with MsCoreUtils::coefWMA based on the variable basicFilterValues$noiseFiltering_method
#' which can be: "MAD", "SuperSmoother" or "wma". After that is either uses a 
#' fixed threshold or a threshold based on the percent of the maximum peak.
#' 
#' @param ms2 A Spectra object.
#' @param basicFilterValues Named list containing the information which filters 
#' to use.
#'
#' @return Returns the cleaned Spectra object.
#'
#' @importFrom MsCoreUtils noise smooth coefWMA
#' @importFrom Spectra addProcessing filterIntensity
#' @noRd
cleaning_ref<- function(ms2, basicFilterValues){
  
  sps_filtered <- switch(
    basicFilterValues$noiseFiltering_method,
    
    "MAD" = addProcessing(ms2, function(x, spectrumMsLevel) {
      n <- MsCoreUtils::noise(x[, "mz"], x[, "intensity"], method = "MAD")
      x[x[, "intensity"] > 3 * n, , drop = FALSE]
    }),
    
    "SuperSmoother" = addProcessing(ms2, function(x, spectrumMsLevel) {
      n <- MsCoreUtils::noise(x[, "mz"], x[, "intensity"], method = "SuperSmoother")
      x[x[, "intensity"] > 3 * n, , drop = FALSE]
    }),
    
    "wma" = addProcessing(ms2, function(x, spectrumMsLevel) {
      if (nrow(x) < 3) return(x)
      x[, 2] <- MsCoreUtils::smooth(x[, 2], cf = MsCoreUtils::coefWMA(hws = 1))
      return(x)
    }),
    
    
    ms2
    
  )

  if(basicFilterValues$min_filtering_type == "fixed"){
    sps_filtered <- filterIntensity(
      sps_filtered, intensity = basicFilterValues$fixed_threshold 
    )
  }
  else if(basicFilterValues$min_filtering_type == "percent"){
    sps_filtered <- filterIntensity(sps_filtered, intensity = function(x, ...) x > max(x, na.rm = TRUE) * 0.01)
  }
  
  return(sps_filtered)
  

}

#' @title Calculation of chain Possiblities
#'
#' @description Calculates all possible formulas and masses from the chains.
#' 
#' @details
#' Divides the modifications in optional and required. Creates a dataframe of
#' combinations from C-atom counts, double bound counts and optional modifications.
#' Builds the formula based on this dataframe and calculated the mass for each
#' formula. 
#' 
#' @param chains List of chains from Lipid Definition.
#'
#' @return returns the calculated masses and formulas of the chains.
#'
#' @importFrom MetaboCoreUtils calculateMass
#' @import stringr
#' @importFrom tidyr extract
#' @noRd
calc_lipid_chain_masses <- function(chains){
  
  optional_mods <- unlist(
    lapply(chains$modifications, \(x) if (x$optional) x$formula),
    use.names = FALSE
  )
  optional_mods <- optional_mods[optional_mods != ""]
  optional_mods <- optional_mods[!is.na(optional_mods)]

  
  required_mods <- unlist(
    lapply(chains$modifications, \(x) if (!x$optional) x$formula),
    use.names = FALSE
  )
  
  required_mods <- required_mods[!required_mods == ""]
  required_mods <- required_mods[!is.null(required_mods)]
  required_mods <- required_mods[!is.na(required_mods)]
  
  mods_combinations <- c(
    list(""),
    unlist(
      lapply(
        seq_along(optional_mods),
        function(k) combn(optional_mods, k, simplify = FALSE)
      ),
      recursive = FALSE
    )
  )
  if(length(mods_combinations) == 0) mods_combinations <- ""
  
  combos <- expand.grid(
    C = seq(
      from = ifelse(chains$even_steps, 
                    ifelse(chains$C_min %% 2 == 0,chains$C_min, chains$C_min+1 ), 
                    chains$C_min),
      to   = chains$C_max,
      by   = ifelse(chains$even_steps, 2, 1)
    ),
    DB = chains$DB_min:chains$DB_max,
    optional_mod = mods_combinations
  )
  mod_H_correction <- length(required_mods)
  
  
  combos$H <- (2 * combos$C) + 2 - (2 * combos$DB) - mod_H_correction
  
  combos$formula <- mapply(
    function(C, H, optional_mod) {
      
      mod_DB <- stringr::str_count(
        paste0(c(optional_mod, required_mods)),
        "="
      )

      modDB <- sum(mod_DB, na.rm = TRUE)
      mod_len <- if (length(optional_mod) == 0) {
        0
      } else {
        length(optional_mod)
      }

      paste0(
        "C", C,
        "H", H - mod_len - modDB,
        paste0(optional_mod, collapse = ""),
        paste0(required_mods, collapse = "")
      )
    },
    combos$C,
    combos$H,
    combos$optional_mod
  )
  combos$formula <- gsub("=", "", combos$formula)

  combos$formula <- unlist(lapply(combos$formula, function(formula){
    matches <- regmatches(formula, gregexpr("([A-Z][a-z]?)([0-9]*)", formula))[[1]]
    element_counts <- data.frame(matches = matches) %>%
      tidyr::extract(
        matches,
        into = c("element", "count"),
        regex = "^([A-Za-z]+)([0-9]*)$"
      ) %>%
      dplyr::mutate(
        count = dplyr::if_else(count == "", "1", count)
      )

    refined_formula <- ""
    for(elem in unique(element_counts$element)){
      elem_count <- sum(as.numeric(element_counts[which(element_counts$element == elem),"count"]))
      refined_formula <- paste0(refined_formula, elem, ifelse(elem_count == 1, "", elem_count))
    }
    refined_formula
  }))

  combos$mass <- suppressWarnings(MetaboCoreUtils::calculateMass(combos$formula))
  combos <- unique(combos)

  if(is.null(optional_mods))combos[c("formula", "mass")]
  else combos[c("formula", "optional_mod", "mass")]
}

