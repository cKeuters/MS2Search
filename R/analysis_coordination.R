
#' Loading MS2 data 
#'
#' @description Loading MS2 data
#'
#' @return returns a spectra object 
#'
#' @import Spectra
#' @import MsBackendMgf
#' @noRd
load_spectra <- function(ms2_file){
  print(ms2_file)
  if(is.null(ms2_file))ms2_file <- "path/data/Test_dataset.mgf"
  #if(is.null(ms2_file))ms2_file <- "path/data/example.mgf"
  
  sp <- Spectra::Spectra(ms2_file, source = MsBackendMgf::MsBackendMgf())
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
run_analysis <- function(settings_list, sp){
  #settings_list <- load_project_settings("path/data/Test_settings.zip")
  
  #-----------------------------------------------------------------------------
  # Basic filters and checks
  #-----------------------------------------------------------------------------

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
  sp <- filterPrecursorMzRange(sp,c(settings_list$basicFilters$mz_min,settings_list$basicFilters$mz_max))
  add_log(paste0("Spectra in mz-Range (", settings_list$basicFilters$mz_min,", ", settings_list$basicFilters$mz_max, "): ", length(sp)))
  
  #check RT
  sp <- getRT(sp)
  sp <- filterRt(sp, rt = c(settings_list$basicFilters$rt_min,settings_list$basicFilters$rt_max))
  add_log(paste0("Spectra in RT-Range (", settings_list$basicFilters$rt_min,", ", settings_list$basicFilters$rt_max, "): ", length(sp)))
  
  #check peak count
  sp <- sp[
    lengths(sp) > settings_list$basicFilters$n_peaks_min &
    lengths(sp) < settings_list$basicFilters$n_peaks_max,
  ]  
  add_log(paste0("Spectra in peak-count-range (", settings_list$basicFilters$n_peaks_min,", ", settings_list$basicFilters$n_peaks_max, "): ", length(sp)))
  
  
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
        
        
        has_ladder <- lengths(ladders) > 0
        
        sp <- sp[has_ladder]
        add_log(paste0("Spectra after ",rule$label,": ", length(sp)))
        

      }else{
        
        ladders <- lapply(peaksData(sp), find_ladder,
                          values = masses, losses_min = 1, losses_max = 1, tol_ppm = settings_list$basicFilters$tolerance_ppm, tol_mz = settings_list$basicFilters$minimum_dalton)
        
        has_ladder <- lengths(ladders) > 0
        
        sp <- sp[has_ladder]
        add_log(paste0("Spectra after ",rule$label,": ", length(sp)))
      }
      
    }else{
      #Searching for specific mz values in the ms2 data
      
      QUERY <- "QUERY scaninfo(MS2DATA) WHERE "
      parts <- c()
      for(mass in masses){
        ppm_da <- as.numeric(mass) * settings_list$basicFilters$tolerance_ppm * 1e-6
        
        if(ppm_da < settings_list$basicFilters$minimum_dalton){
          parts <- c(parts,paste0("MS2PROD=", mass,":TOLERANCEMZ=", settings_list$basicFilters$minimum_dalton))
        } else {
          parts <- c(parts,paste0("MS2PROD=", mass,":TOLERANCEPPM=", settings_list$basicFilters$tolerance_ppm))
        }
      }
      
      QUERY <- paste(QUERY, paste(parts, collapse = " OR "))
      add_log(QUERY)
      sp <- query(sp, QUERY)
      add_log(paste0("Spectra after ",rule$label,": ", lengths(sp$FEATURE_ID)))
    }

  }
  
  #------------------------- creating dataframe from sp ----------------------
  
  results <- spectraData(sp)
  results <- results[,c("FEATURE_ID","msLevel", "rtime", "precursorMz","POLARITY", "ION")]
  results$score <- NA
  results
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



#-----------------------------------------------------------------------------
# Logging
#-----------------------------------------------------------------------------

#' @noRd
add_log <- function(txt) {
  txt <- paste(txt , "\n")
  cat(txt)
  log_text(paste(log_text(), txt))
}



