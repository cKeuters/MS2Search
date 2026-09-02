
#' Loading MS2 data 
#'
#' @description Loading MS2 data
#'
#' @return returns a spectra object 
#'
#' @import Spectra
#' @import MsBackendMgf
#' @noRd
load_spectra <- function(ms2_file, add_log){
  
  #if(is.null(ms2_file))ms2_file <- "path/data/Test_dataset.mgf"
  if(is.null(ms2_file))ms2_file <- "path/data/example.mgf"
  print(ms2_file)
  
  sp <- Spectra::Spectra(ms2_file, source = MsBackendMgf::MsBackendMgf())
  sp$FEATURE_ID <- as.numeric(sp$FEATURE_ID)
  sp <- getRT(sp)
  
  add_log(paste("Loaded MS2-file with ", length(sp)," spectra."))
  sp
}

#' Analysis coordination 
#'
#' @description A function to run the analysis
#'
#' @return returns a matrix 
#'
#' @import Spectra
#' @import MsBackendMgf
#' @import SpectraQL
#' @noRd
run_analysis <- function(settings_list, sp, add_log){
  #settings_list <- load_project("path/data/test_settings.zip")$settings_list
  
  #-----------------------------------------------------------------------------
  # Basic filters and checks
  #-----------------------------------------------------------------------------
  results <- data.frame(FEATURE_ID = spectraData(sp)$FEATURE_ID)

  #check MS-Level 2 & empty spectra
  sp <- sp |>
    filterMsLevel(2) |>
    filterEmptySpectra() 
  add_log(paste("Non-empty spectra with MS-level 2: ", length(sp)))

  #check mode
  pol_col <- intersect(c("POLARITY", "polarity", "IONMODE", "ionmode"),spectraVariables(sp))
  if(settings_list$project_info$mode == "pos") mode <- c("POS", "POSITIVE", "positive") else mode <- c("NEG", "NEGATIVE", "negative")
  if (length(pol_col) > 0) {sp <- sp[sp[[pol_col[1]]] %in% mode, ]}
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
    
  sp_cleaned <- cleaning_ref(sp)
  
    
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
          find_ladder(
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
        peaks_data,
        sp$precursorMz
      )

      retained_vector <- lengths(ladders) > 0
      add_log(paste0("Spectra after ",rule$label,": ", sum(retained_vector)))
      results[rule$label] <- retained_vector
      if(rule$rule_priority == "required"){
        results <- results[retained_vector,]
        sp <- sp[retained_vector,]
        sp_cleaned <- sp_cleaned[retained_vector,]
        
      }
      
    }else{
      #Searching for specific mz values in the ms2 data
      retained_vector <- unlist(lapply(peaks_data, findMZpeak,
                values = masses, tol_ppm = settings_list$basicFilters$tolerance_ppm, tol_mz = settings_list$basicFilters$minimum_dalton))
      add_log(paste0("Spectra after ",rule$label,": ", sum(retained_vector)))
      results[rule$label] <- retained_vector
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
    for(ld in settings_list$lipid_definition){
      masses <- calc_lipid_chain_masses(ld)$mass
      
      retained_vector <- unlist(lapply(peaks_data, findMZpeak,
                                       values = masses, 
                                       tol_ppm = settings_list$basicFilters$tolerance_ppm, 
                                       tol_mz = settings_list$basicFilters$minimum_dalton))
      add_log(paste0("Spectra after Chain ",chain_id,": ", sum(retained_vector)))
      results[paste("chain",chain_id)] <- retained_vector
      
      chain_id <- chain_id + 1
      if(ld$priority == "required"){
        
        results <- results[retained_vector,]
        sp <- sp[retained_vector,]
        sp_cleaned <- sp_cleaned[retained_vector,]
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
      sapply(settings_list$lipid_definition, \(x) priority_weights[x$priority])
    )
    priorities <- unname(sapply(priorities, function(x) if (length(x) == 0) 0 else x))

    matches <- as.matrix(results[, 2:ncol(results)])
    matches <- apply(matches,2, function(x) as.integer(x))
    
    print(head(matches))
    print(sum(priorities))
    results$score <- rowSums((sweep(matches, MARGIN=2,priorities, `*`))) / sum(priorities)
    
    as.data.frame(results)
}




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

#' Finding sequential neutral loss pattern 
#'
#' @description Finding sequential neutral loss pattern
#'
#' @return returns the values of the ladder
#'
#' @import MsCoreUtils
#' @noRd
find_ladder <- function(peaks, precursor = NULL, values, losses_min, losses_max, tol_ppm, tol_mz, rule_min_mz, rule_max_mz, just_after_precursor = FALSE, noise_filter = TRUE) {
  
  mz <- peaks[,"mz"]
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
  !is.na(any(idx))
}

#' Noise filtering 
#'
#' @description Noise filtering from https://doi.org/10.1021/acs.analchem.5c02109
#'
#' @return returns the cleaned spectra
#'
#' @import MASS
#' @import tidyverse
#' @noRd
cleaning_ref<- function(ms2){
  sps_filtered <- addProcessing(ms2, function(x, spectrumMsLevel) {
    n <- MsCoreUtils::noise(x[, "mz"], x[, "intensity"], method = "SuperSmoother")
    x[x[, "intensity"] > 3 * n, , drop = FALSE]
  })
  
  # sps_filtered <- filterIntensity(
  #   sps_filtered,
  #   intensity = 20  
  # )
  # 
  return(sps_filtered)
  

}

#' Chains to masses 
#'
#' @description calculates the masses from the lipid definition
#'
#' @return returns the possible masses
#'
#' @import MetaboCoreUtils
#' @noRd
calc_lipid_chain_masses <- function(ld){
  
  optional_mods <- unlist(
    lapply(ld$modifications, \(x) if (x$optional) x$formula),
    use.names = FALSE
  )
  required_mods <- unlist(
    lapply(ld$modifications, \(x) if (!x$optional) x$formula),
    use.names = FALSE
  )
  mods_combinations <- unlist(
    lapply(
      0:length(optional_mods),
      function(k) {
        combn(optional_mods, k, simplify = FALSE)
      }
    ),
    recursive = FALSE
  )
  if(length(mods_combinations) == 0) mods_combinations <- ""
  combos <- expand.grid(
    C  = ld$C_min:ld$C_max,
    DB = ld$DB_min:ld$DB_max,
    optional_mod = mods_combinations)
  
  
  mod_H_correction <- length(required_mods)
  
  combos$H <- 2 * combos$C - 2 * combos$DB - mod_H_correction
  

  combos$formula <- mapply(
    function(C, H, optional_mod) {
      paste0(
        "C", C,
        "H", H - length(optional_mod) +1,
        paste0(optional_mod, collapse = ""),
        paste0(required_mods, collapse = "")
      )
    },
    combos$C,
    combos$H,
    combos$optional_mod
  )
    
  combos$mass <- calculateMass(combos$formula)
  
  if(is.null(optional_mods))combos[c("formula", "mass")]
  else combos[c("formula", "optional_mod", "mass")]
}

