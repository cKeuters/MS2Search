
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
  #if(is.null(ms2_file))ms2_file <- "path/data/Test_dataset.mgf"
  if(is.null(ms2_file))ms2_file <- "path/data/example.mgf"
  
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
  #settings_list <- load_project("path/data/Test-Settings.zip")$settings_list
  
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
    
    sp_cleaned <- cleaning_ref(sp)
   
    
  #-----------------------------------------------------------------------------
  # Rule filters
  #-----------------------------------------------------------------------------

  for(rule in settings_list$specific_filters){
    
    if(rule$filter_noise)  peaksData <- peaksData(sp)
    else  peaksData <- peaksData(sp_cleaned)
    
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
      ladders <- lapply(peaksData, find_ladder,
                        values = masses, 
                        losses_min = rule$min_repeat, 
                        losses_max = rule$max_repeat, 
                        rule_min_mz = rule$rule_mz_min, 
                        rule_max_mz = rule$rule_mz_max, 
                        tol_ppm = settings_list$basicFilters$tolerance_ppm, 
                        tol_mz = settings_list$basicFilters$minimum_dalton)
      
      retained_vector <- lengths(ladders) > 0
      add_log(paste0("Spectra after ",rule$label,": ", sum(retained_vector)))
      results[rule$label] <- retained_vector
      if(rule$priority == "required"){
        results <- results[retained_vector,]
        sp <- sp[retained_vector,]
      }
      
    }else{
      #Searching for specific mz values in the ms2 data
      
      retained_vector <- unlist(lapply(peaks, findMZpeak,
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
find_ladder <- function(peaks, values, losses_min, losses_max, tol_ppm, tol_mz, rule_min_mz, rule_max_mz, just_after_precursor = FALSE, noise_filter = TRUE) {
  
  #if(noise_filter)peaks <- compute_Noise(peaks)
  mz <- peaks[,"mz"]

  if(!is.null(rule_min_mz))mz <- mz[mz > rule_min_mz]
  if(!is.null(rule_max_mz))mz <- mz[mz < rule_max_mz]
  

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

findMZpeak <- function(peaks, values, tol_ppm, tol_mz){
  mz <- peaks[,"mz"]
  values <- as.numeric(values)
  idx <- closest(values, mz, tolerance = tol_mz, ppm = tol_ppm, duplicates = "closest")
  !is.na(any(idx))
}

compute_Noise <- function(peaks){
  peaks_df <- as.data.frame(peaks)
  groups <- split(
    peaks_df,
    cut(
      peaks_df$mz,
      breaks = seq(0, ceiling(max(peaks_df$mz)/100) * 100, 100),
      include.lowest = TRUE
    )
  )
  
  filtered_groups <- lapply(groups, function(set) {
    med <- median(set$intensity)
    noise <- median(abs(set$intensity - med))
    
    set$snr <- set$intensity / noise
    })
  
  peaks_df$snr <- do.call(c, filtered_groups)
  
  peaks[!peaks_df$snr > 3, , drop = FALSE]
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

  ms2_tb <- tibble::tibble(
    id = ms2$FEATURE_ID,
    dataOrigin = ms2$dataOrigin,
    scanIndex = ms2$scanIndex,
    polarity = ms2$polarity,
    msLevel = ms2$msLevel,
    precursorMz = ms2$precursorMz,
    rtime = ms2$rtime,
    #inchi = ms2$inchi,
    precursorIntensity = ms2$precursorIntensity,
    ms2_raw = as.list(peaksData(ms2, columns = c("mz", "intensity")))
  ) %>%
    dplyr::mutate(ms2_raw = purrr::map(ms2_raw, function(t) {
      t %>%
        as.data.frame() %>%
        dplyr::mutate(mz = as.numeric(mz), intensity = as.numeric(intensity))
    })) %>%
    dplyr::mutate(ms2_RLM = purrr::map(ms2_raw, function(t){
      t <- t %>%
        dplyr::filter(intensity>0) %>%
        dplyr::arrange(intensity) %>%
        dplyr::mutate(order = seq(nrow(.)))
      p75 <- t %>%
        dplyr::pull(intensity) %>%
        quantile(0.75)
      u <- t %>%
        dplyr::filter(intensity < p75)
      if(nrow(u)>3){
        mdl <- MASS::rlm(data = u, log10(intensity) ~ order)
        t_filt<- t %>%
          dplyr::mutate(pred = predict(mdl, t)) %>%
          dplyr::mutate(res = abs(pred - log10(intensity))) %>%
          dplyr::filter(res > sd(residuals(mdl))*3) %>%
          dplyr::arrange(mz)
        if (nrow(t_filt) > 0) {
          t<- t_filt
        } else {
          t<- t %>% dplyr::arrange(mz)
        }
      } else {
        return(t %>% dplyr::arrange(mz)) # Not enough points to model
      }
    }))
  ms2_cl <- ms2_tb %>%
    dplyr::mutate(precursorIntensity = precursorIntensity) %>%
    dplyr::mutate(
      mz = purrr::map(ms2_RLM, ~ .x$mz),
      intensity = purrr::map(ms2_RLM, ~ .x$intensity)
    ) %>%
    dplyr::select(
      id, dataOrigin, scanIndex, rtime, precursorMz, precursorIntensity,
      polarity, mz, intensity
    ) %>%
    dplyr::filter(purrr::map_lgl(mz, ~ length(.x) > 0)) %>%
    as.data.frame() %>%
    Spectra()
  return(ms2_cl)
}
