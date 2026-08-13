
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
  print(ms2_file)
  if(is.null(ms2_file))ms2_file <- "path/data/Test_dataset.mgf"
  #if(is.null(ms2_file))ms2_file <- "path/data/example.mgf"
  
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
  #settings_list <- load_project_settings("path/data/Test_settings.zip")
  
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
    
  #check peak count
    all_PeakCount <- lengths(sp)
    retained_vector <- all_PeakCount >= settings_list$basicFilters$n_peaks_min & all_PeakCount <= settings_list$basicFilters$n_peaks_max
    retained_vector[is.na(retained_vector)] <- FALSE
    add_log(paste0("Spectra in peak-count-range (", settings_list$basicFilters$n_peaks_min,", ", settings_list$basicFilters$n_peaks_max, "): ", sum(retained_vector)))
    results$peakCount_filter <- retained_vector
    if(settings_list$basicFilters$n_peaks_priority == "required"){
      results <- results[retained_vector,]
      sp <- sp[retained_vector,]
    }

  #-----------------------------------------------------------------------------
  # Rule filters
  #-----------------------------------------------------------------------------

  for(rule in settings_list$specific_filters){
    
    masses <- vapply(rule$entries, function(e) {
      resolve_mass(list(input_type = rule$input_type, value = e$value, is_mzPeak = rule$is_mzPeak),settings_list$project_info$mode)
    }, numeric(1))
    masses <- abs(as.numeric(masses[!is.na(masses)]))
    
    if(!rule$is_mzPeak){
      
      #searching for multiple neutral losses
      if(rule$sequential & rule$max_repeat > 1){
        ladders <- lapply(peaksData(sp), find_ladder,
               values = masses, losses_min = rule$min_repeat, losses_max = rule$max_repeat, tol_ppm = settings_list$basicFilters$tolerance_ppm, tol_mz = settings_list$basicFilters$minimum_dalton)
      }else{
        ladders <- lapply(peaksData(sp), find_ladder,
                          values = masses, losses_min = 1, losses_max = 1, tol_ppm = settings_list$basicFilters$tolerance_ppm, tol_mz = settings_list$basicFilters$minimum_dalton)
      }
      retained_vector <- lengths(ladders) > 0
      add_log(paste0("Spectra after ",rule$label,": ", sum(retained_vector)))
      results[rule$label] <- retained_vector
      if(rule$priority == "required"){
        results <- results[retained_vector,]
        sp <- sp[retained_vector,]
      }
      
    }else{
      #Searching for specific mz values in the ms2 data
      
      retained_vector <- unlist(lapply(peaksData(sp), findMZpeak,
                values = masses, tol_ppm = settings_list$basicFilters$tolerance_ppm, tol_mz = settings_list$basicFilters$minimum_dalton))
      add_log(paste0("Spectra after ",rule$label,": ", sum(retained_vector)))
      results[rule$label] <- retained_vector
      if(rule$priority == "required"){
        results <- results[retained_vector,]
        sp <- sp[retained_vector,]
      }

    }

  }
  
  #------------------------- creating dataframe  ----------------------

    results$score <- rowSums(
      as.matrix(results[, 2:ncol(results)]),
      na.rm = TRUE
    ) / (ncol(results) - 1)
    
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
find_ladder <- function(peaksData, values, losses_min, losses_max, tol_ppm, tol_mz, just_after_precursor = FALSE) {
  mz <- peaksData[, "mz"]
  
  values <- as.numeric(values)
  n_peaks <- length(mz)
  if (n_peaks < losses_min) return(numeric(0))
  
  max_check <- losses_max + 1L
  live <- data.frame(anchor = seq_len(n_peaks), offset = 0)
  if(just_after_precursor) live <- live[live$anchor == n_peaks,]
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

findMZpeak <- function(peaksData, values, tol_ppm, tol_mz){
  mz <- peaksData[, "mz"]
  values <- as.numeric(values)
  idx <- closest(values, mz, tolerance = tol_mz, ppm = tol_ppm, duplicates = "closest")
  !is.na(any(idx))
}


