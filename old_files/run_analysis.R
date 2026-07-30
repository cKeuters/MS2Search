
#https://rformassspectrometry.github.io/Spectra/articles/Spectra.html
filterMS2data <- function(ms2_file, rule_sets = NULL){
  
  sp <- Spectra(ms2_file, source = MsBackendMgf())
  n_spectra <- length(sp)
  
  sp <- sp |>
    filterMsLevel(2) |>
    #filterPolarity(polarity) |>
    filterEmptySpectra() 
  
  pol_col <- intersect(
    c("POLARITY", "polarity", "IONMODE", "ionmode"),
    spectraVariables(sp)
  )
  if(state$mode == "pos") mode <- c("POS", "POSITIVE", "positive") else mode <- c("NEG", "NEGATIVE", "negative")
  if (length(pol_col) > 0) {
    sp <- sp[sp[[pol_col[1]]] %in% mode, ]
  }
    
  sp <- filterPrecursorMzRange(sp,c(basicFilterValues$mz_min,basicFilterValues$mz_max))
  
  sp <- getRT(sp)
  sp <- filterRt(sp, rt = c(basicFilterValues$rt_min,basicFilterValues$rt_max))
  
  print(paste0("Basic Filters completed: ", n_spectra, " -> ", length(sp) ))
  n_spectra <- length(sp)
  
  if(!is.null(rule_sets)){
    
    #create Neutral Loss spectra
    nl_filter <- PrecursorMzParam(
      filterPeaks = "belowPrecursor",
      msLevel = c(2L, NA_integer_),
      ppm = 10,
      tolerance = 0
    )
    sps_nl <- neutralLoss(sp,nl_filter)
    print("Neutral Loss Spectra created")
    
    overall_result <- rep(TRUE, length(sp))
    
    for(rule in rule_sets){
      if (rule$input_type == "mass") {
        rule$entries <- as.numeric(unlist(rule$entries))
      } else {
        names_vec <- unlist(rule$entries)
        idx <- match(names_vec, METABOLITES$name)
        if (any(is.na(idx))) {
          warning("Unbekannte Metabolit-Namen: ", paste(names_vec[is.na(idx)], collapse = ", "))
        }
        rule$entries <- METABOLITES$mass[idx]
      }
      
      containing_list <- lapply(
        unlist(rule$entries),
        function(mass) containsNeutralLoss(sp, neutralLoss = mass)
      )
      containing_matrix <- do.call(cbind, containing_list)
      True_rows <- rowSums(containing_matrix, na.rm = TRUE) > 0
      
      rule_result <- rep(FALSE, length(sp))
      rule_result[True_rows] <- filter_NL_count(sps_nl[True_rows, ], rule)
      
      overall_result <- overall_result & rule_result
      
    }

    if(sum(overall_result) == 0){
      result_text <- paste0("Rule Analysis completed: ", n_spectra, " -> 0")
    }else{
      sp_rules <- sp[overall_result,]
      result_text <- paste0("Rule Analysis completed: ", n_spectra, " -> ", length(sp_rules), "\n",
      "Feature IDs:", paste(sp_rules$FEATURE_ID, collapse = ", "), "\n")
    }
    cat(result_text)
    return(result_text)
  }
  return("No Result Text")
}

filter_NL_count <- function(sps_nl, rule) {
  pd_nl <- peaksData(sps_nl)
  entries_vec <- unlist(rule$entries)
    
    vapply(pd_nl, function(peaks) {
      nl_vals <- as.numeric(peaks[, "mz"])
      
      if (length(nl_vals) == 0) {
        counter <- 0
      } else {
        matches <- vapply(nl_vals, function(v) {
          any(abs(v - entries_vec) <= 0.1)
        }, logical(1))
        counter <- sum(matches)
      }
      
      counter >= rule$min_repeat & counter <= rule$max_repeat
    }, logical(1))
  
}



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
