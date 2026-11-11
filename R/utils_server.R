
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
#'
#' @return Returns the list of all parameters
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
#' @importFrom jsonlite write_json
#' @noRd
setup_and_save_project <- function(settings_list, res_settings, results, spectra, output_zip_path_file) {
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
    
    jsonlite::write_json(res_settings, file.path(tmp_dir, "res_settings.json"), auto_unbox = TRUE)
    
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
#' @importFrom jsonlite fromJSON
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
  res_settings <- NULL
  if (("res_settings.json" %in% list.files(tmp_dir))) {
    res_settings <- jsonlite::fromJSON(
      file.path(tmp_dir, "res_settings.json"),
      simplifyVector = FALSE  
    )
  }
  results <- NULL
  if (("results.csv" %in% list.files(tmp_dir))) {
    results <- read.csv(file.path(tmp_dir, "results.csv"))
    
    results[] <- lapply(results, function(col) {
      if (is.character(col)) {
        lapply(col, function(x) {
          if (grepl("/", x, fixed = TRUE)) {
            as.numeric(strsplit(x, "/", fixed = TRUE)[[1]])
          } else {
            as.numeric(x)
          }
        })
      } else {
        col
      }
    })
    colnames(results) <- gsub("\\.", " ", colnames(results))
  }
  spectra <- NULL
  if (("spectra.mgf.gz" %in% list.files(tmp_dir))) {
    spectra <- Spectra::Spectra(file.path(tmp_dir, "spectra.mgf.gz"), source = MsBackendMgf::MsBackendMgf())
  }
  message("File loaded")
  unlink(tmp_dir, recursive = TRUE)
  return(list(settings_list = settings_list, res_settings = res_settings, results = results, spectra = spectra))
}



#' Get a single spectrum from a Spectra object by feature ID
#'
#' @param spectra A `Spectra` object containing multiple spectra.
#' @param feature_id Character or numeric. The ID to match against.
#'
#' @return A `Spectra` object subset to the matching spectrum/spectra.
#'   Throws an informative error if no match is found.
#'
#' @importFrom Spectra spectraData
get_spectrum_by_id <- function(spectra, feature_id) {
  
  spectraData <- Spectra::spectraData(spectra)
  
  if (!"FEATURE_KEY" %in% colnames(spectraData)) {
    stop(sprintf(
      "Column '%s' not found in spectraData. Available columns: %s",
      "FEATURE_KEY", paste(colnames(spectraData), collapse = ", ")
    ))
  }
  
  match_idx <- which(spectraData$FEATURE_KEY == feature_id)
  
  if (length(match_idx) == 0) {
    stop(sprintf("No spectrum found for %s = '%s'", "", feature_id))
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
#' @importFrom jsonlite toJSON
#' @noRd
create_brkrrba_file <- function(settings_list, file_path){
  
  searchedSpecies <- settings_list$lipid_definition$shorthand
  searchedSpecies <- ifelse(searchedSpecies == "", "Unk" ,searchedSpecies)
  annotationName <- sprintf("%s {carbonsInAllChains}:{DBEinAllChains}", searchedSpecies )
  
  #Variables of optional modifictaions..........................................
  
  chain_opt_mods <- unname(lapply(seq_along(settings_list$lipid_definition$chains), function(i) {
    x <- settings_list$lipid_definition$chains[[i]]
    if(length(x$modifications) == 0)return(NULL)
    
    mods <- x$modifications
    formulas <- lapply(mods, function(m){if(m$optional) m$formula})

    counts <- tibble(formula = unlist(formulas)) %>% count(formula)
    lapply(seq_len(nrow(counts)), function(j) {
      list(
        name = gsub( "=", "DB_",paste0(counts$formula[j], "_", i)),
        valueType = "intRange",
        value = c(0, counts$n[j])
      )
    })
    
  }))
  
  #Formulas of the chains for the descriptors ..................................

  chain_formulas <- lapply(seq_along(settings_list$lipid_definition$chains), function(i) {
    x <- settings_list$lipid_definition$chains[[i]]

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
      "C{carbonsPerChain", i, "} + H({carbonsPerChain", i, "}*2) + H2 - H({DBEPerChain", i, "}*2)",
      if (length(req_mods) > 0){
        matches <- regmatches(req_mods, gregexpr("([A-Z][a-z]?)([0-9]*)", req_mods))[[1]]
        paste0(paste0(" + ", paste(matches, collapse = " + ")), " - H",length(req_mods))
      }
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
      formula <- paste0(formula, paste(lapply(vars[grepl("DB_", vars)], function(x){
        paste0(" - H(",x,")")
      })))
    }
    formula
    
  })
  
  #descriptors of the chains ...................................................
  
  chain_descriptors <- lapply(seq_along(chain_formulas),function(i){
    list(
      label = paste0("Chain ", i),
      query = sprintf(
        "QUERY scaninfo(MS2DATA) WHERE MS2PROD=formula(%s):TOLERANCEMZ=%s",
        chain_formulas[i], settings_list$basicFilters$minimum_dalton
      ),
      signalInfoLabel = paste0("{carbonsPerChain", i, "}:{DBEPerChain", i, "}")
    )
  })
  
  #descriptors of the rules  ..................................................
  
  #function to get the combination of masses for the sequential masses
  get_sums <- function(masses, k) {
    if (k == 1) return(masses)
    combos <- expand.grid(rep(list(seq_along(masses)), k))
    unique(apply(combos, 1, function(idx) sum(masses[idx])))
  }
  
  rule_descriptors <- lapply(settings_list$specific_filters, function(rule){

      masses <- vapply(rule$entries, function(e) {
        resolve_mass(list(input_type = rule$input_type, value = e$value, is_mzPeak = rule$is_mzPeak),settings_list$project_info$mode)
      }, numeric(1))
      NL_or_peak <- if (rule$is_mzPeak) "MS2PROD" else "MS2NL"
      
      if (isTRUE(rule$sequential)) {
        
        repeat_queries <- lapply(seq_len(rule$min_repeat), function(k) {
          
          vals <- sort(get_sums(masses, k))
          
          paste0(
            NL_or_peak,
            "=",
            vals,
            collapse = " OR "
          )
        })
        
        query <- paste0(
          "QUERY scaninfo(MS2DATA) WHERE ",
          paste(repeat_queries, collapse = " AND "),
          ":TOLERANCEMZ=",
          settings_list$basicFilters$minimum_dalton
        )
        
      } else {
        
        query <- paste0(
          "QUERY scaninfo(MS2DATA) WHERE ",
          paste0(NL_or_peak, "=", masses, collapse = " OR "),
          ":TOLERANCEMZ=",
          settings_list$basicFilters$minimum_dalton
        )
        
      }

      list(
        label = rule$label,
        query = query,
        signalInfoLabel = paste0(NL_or_peak, ": ", paste(sapply(rule$entries, `[[`, "value"), collapse = ", "))
      )
    })
  
  #Precursor validation formula  ..................................................
  
  #Precursor validation from the validation list
  validation_names <- sapply(settings_list$lipid_definition$validations, function(x) x$name)
  chain_ids <- sapply(settings_list$lipid_definition$chains, function(x) x$chain_id)
  if("precursor validation" %in% validation_names){
    pre_val <- settings_list$lipid_definition$validations[which(validation_names == "precursor validation")]
    parts <- c()
    
    for (com in pre_val[[1]]$components) {
      if (com == "Backbone") {
        parts <- c(parts, "{backbone}")
      } else {
        i <- which(chain_ids == com)
        parts <- c(parts, chain_formulas[i])
      }
    }
    neutralPrecursorFormula <- paste(parts, collapse = " + ")
    
    if(is.null(settings_list$lipid_definition$head_group) || settings_list$lipid_definition$head_group == ""){ 
      neutralPrecursorFormula <- paste(neutralPrecursorFormula, "-Na")}

    if(pre_val[[1]]$adjustments != ""){
      add <- ifelse(substr(pre_val[[1]]$adjustments,1,1) == "-", " - ", " + ")
      matches <- regmatches(pre_val[[1]]$adjustments, gregexpr("([A-Z][a-z]?)([0-9]*)", pre_val[[1]]$adjustments))[[1]]
      element_counts <- data.frame(matches = matches) %>%
        tidyr::extract(
          matches,
          into = c("element", "count"),
          regex = "^([A-Za-z]+)([0-9]*)$"
        ) %>%
        dplyr::mutate(
          count = dplyr::if_else(count == "", "1", count)
        )
      adjustments <- ""
      for(elem in unique(element_counts$element)){
        elem_count <- sum(as.numeric(element_counts[which(element_counts$element == elem),"count"]))
        adjustments <- paste0(adjustments, add, elem, ifelse(elem_count == 1, "", elem_count))
      }
      neutralPrecursorFormula <- paste(neutralPrecursorFormula, adjustments)
      }


  }
  else neutralPrecursorFormula <- paste0("{backbone} + ", paste(chain_formulas, collapse = " + "))
  
  #Molecular annotation name  ..................................................
  
  chain_opt_mod_names <- unlist(lapply(chain_opt_mods, function(x)
      sapply(x, `[[`, "name")
    ))
    molecular_annotationName <- searchedSpecies 
    chain_strings <- c()
    chain_seq <- seq_along(settings_list$lipid_definition$chains)
    for (i in chain_seq) {
      mods_i <- chain_opt_mod_names[
        grepl(paste0("_", i, "$"), chain_opt_mod_names)
      ]
      mod_elements <- sub(
        "H.*", "",
        sub(".*?_", "", sub("_[^_]+$", "", mods_i))
      )

      chain_strings <- c(
        chain_strings,
        sprintf(
          "{carbonsPerChain%d}:{DBEPerChain%d};%s",
          i, i,
          paste(
            paste0(mod_elements, "{", mods_i, "}"),
            collapse = ","
          )
        )
      )
    }
    molecular_annotationName <- paste(
      molecular_annotationName,
      paste(chain_strings, collapse = "/")
    )
    
    #Sublevel of the chains  ..................................................
    
    req_chains <- sum(
      sapply(settings_list$lipid_definition$chains, function(x) x$priority) == "required"
    )
    
    chain_subLevel <- list(
      level = "MOLECULAR_SPECIES_LEVEL",
      annotationName = molecular_annotationName,
      descriptors = chain_descriptors,
      requiredSignalsForAnnotation = req_chains
    )
    
    #Setting the first level to rules or chains, if no rules ...................
    
    has_rules <- length(rule_descriptors) > 0
    
    if(has_rules) {
      req_rules <- sum(
        sapply(settings_list$specific_filters, function(x) x$rule_priority) == "required"
      )
    }else req_rules <- 0
    
    validation_rule <- paste(
      paste0("{carbonsPerChain", chain_seq, "}", collapse = " + "),
      "== {carbonsInAllChains} &&",
      paste0("{DBEPerChain", chain_seq, "}", collapse = " + "),
      "== {DBEinAllChains}"
    )
    
    level <- list(
      level = "Species Level",
      annotationName = if (has_rules) annotationName else molecular_annotationName,
      descriptors = if (has_rules) rule_descriptors else chain_descriptors,
      requiredSignalsForAnnotation = if (has_rules) req_rules else req_chains,
      subLevel = if (has_rules) chain_subLevel,
      validationRule = validation_rule
    )

    #Combining all elements to the full structure ..............................
    
    lipid_def <- list(
    variables = c(list(
      list(
        name = "backbone",
        valueType = "formula",
        value = ifelse(is.null(settings_list$lipid_definition$head_group) || settings_list$lipid_definition$head_group == "" ,"Na", settings_list$lipid_definition$head_group)
      ),
      list(name = "numberOfChains",
           valueType = "int",
           value = length(settings_list$lipid_definition$chains)
           ),
      list(
        name = "carbonsInAllChains",
        valueType = "intRange",
        value = c(sum(unlist(lapply(settings_list$lipid_definition$chains, function(x){x$C_min}))),
                  sum(unlist(lapply(settings_list$lipid_definition$chains, function(x){x$C_max}))))
      ),
      list(
        name = "DBEinAllChains",
        valueType = "intRange",
        value = c(sum(unlist(lapply(settings_list$lipid_definition$chains, function(x){x$DB_min}))),
                  sum(unlist(lapply(settings_list$lipid_definition$chains, function(x){x$DB_max}))))
      )
    ),
    unname(Map(function(x, i) {
      list(
        name = paste0("carbonsPerChain", i),
        valueType = "intRange",
        value = if (isTRUE(x$even_steps)) {
            c(x$C_min, x$C_max, 2)
          } else {
            c(x$C_min, x$C_max)
          }
      )
    },
    settings_list$lipid_definition$chains,
    seq_along(settings_list$lipid_definition$chains)))
    ,
    unname(Map(function(x, i) {
      list(
        name = paste0("DBEPerChain", i),
        valueType = "intRange",
        value = c(x$DB_min, x$DB_max)
      )
    },
    settings_list$lipid_definition$chains,
    seq_along(settings_list$lipid_definition$chains))),
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
            annotationConfiguration = level
            
          ))
      )
      )
  )
  
    #Saving file ...............................................................
    
  json_file <- file_path
  write(toJSON(lipid_def, auto_unbox = TRUE, pretty = TRUE), json_file)
  
  txt <- readLines(json_file)
  txt <- gsub('"([^"]+)":', '\\1:', txt)
  writeLines(txt, json_file)
  
}


#-------------------------------------------------------------------------------
# Result Detail calculations
#-------------------------------------------------------------------------------

#' @title Calculate Formula from mass
#' 
#' @param mass Mass for the formula calculation
#' 
#' @importFrom Rdisop decomposeMass
#' @importFrom MetaboCoreUtils mz2mass
#' @noRd
calculateFormulaFromMass <- function(mass, adduct, ppm = 5,
                                     max_el = "C100H150N10O20P2S2") {
  
  
  neutral_mass <- MetaboCoreUtils::mz2mass(mass, adduct = adduct)[1]
  
  res <- Rdisop::decomposeMass(
    mass        = neutral_mass,
    ppm         = ppm,
    z           = 0,
    elements    = Rdisop::initializeElements(c("C", "H", "N", "O", "P", "S")),
    maxElements = max_el
  )
  
  if (is.null(res) || length(res$formula) == 0) {
    return(data.frame(formula = "-", score = NA_real_,
                      exactmass = NA_real_, ppm_error = NA_real_))
  }
  

 cand <- data.frame(
    formula   = res$formula,
    score     = res$score,
    exactmass = res$exactmass,
    valid     = res$valid,
    stringsAsFactors = FALSE
  )
  cand <- cand[cand$valid == "Valid", , drop = FALSE]
  if (nrow(cand) == 0) return(cand)
  
  # element counts for plausibility filtering
  count_el <- function(f, el) {
    m <- regmatches(f, regexpr(paste0(el, "(?![a-z])\\d*"), f, perl = TRUE))
    if (length(m) == 0) return(0)
    n <- sub(el, "", m)
    if (n == "") 1 else as.numeric(n)
  }
  cand$C <- vapply(cand$formula, count_el, numeric(1), "C")
  cand$H <- vapply(cand$formula, count_el, numeric(1), "H")
  cand$N <- vapply(cand$formula, count_el, numeric(1), "N")
  cand$O <- vapply(cand$formula, count_el, numeric(1), "O")
  cand$P <- vapply(cand$formula, count_el, numeric(1), "P")
  cand$S <- vapply(cand$formula, count_el, numeric(1), "S")
  
  keep <- cand$C > 0 &
    cand$H / cand$C >= 0.2 & cand$H / cand$C <= 3.1 &
    cand$N / cand$C <= 1.3 & cand$O / cand$C <= 1.2 &
    cand$P / cand$C <= 0.3 & cand$S / cand$C <= 0.8
  cand <- cand[keep, , drop = FALSE]
  
  
  cand$ppm_error <- (as.numeric(cand$exactmass) - as.numeric(neutral_mass)) /
    as.numeric(neutral_mass) * 1e6
  cand[order(abs(cand$ppm_error)), , drop = FALSE]
}




