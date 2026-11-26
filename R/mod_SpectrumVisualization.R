#' SpectrumVisualization UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
mod_SpectrumVisualization_ui <- function(id) {
  ns <- NS(id)
  tagList(
    checkboxGroupInput(ns("overlays"), NULL, inline = TRUE,
                       choices = c("Neutral losses" = "nl", "Chains" = "chain", "Validations" = "val")),
    plotly::plotlyOutput(ns("spectrum_plot"))
  )
}
#' SpectrumVisualization Server Functions
#'
#' @noRd 
mod_SpectrumVisualization_server <- function(id, spectrum_r, settings_list_r, results_r,
                                set_isotope_sum_r, ctx_r) {
  moduleServer(id, function(input, output, session) {
    ctx <- ctx_r      
    
    built <- reactive({
      req(spectrum_r(), nrow(results_r()) > 0)
      plottingSpectra(spectrum_r(), set_isotope_sum_r(), settings_list_r(), results_r(),
                      ctx(), show = isolate(input$overlays))
    })
    
    output$spectrum_plot <- plotly::renderPlotly(built()$fig)
    
    proxy <- plotly::plotlyProxy("spectrum_plot", session)
    
    observeEvent(input$overlays, ignoreNULL = FALSE, ignoreInit = TRUE, {
      idx <- built()$idx
      for (grp in names(idx)) {
        if (length(idx[[grp]]) == 0) next
        plotly::plotlyProxyInvoke(proxy, "restyle",
                                  list(visible = grp %in% input$overlays),
                                  as.list(idx[[grp]]))
      }
    })
  })
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
#' @importFrom data.table data.table rbindlist
#' 
#' @noRd
neutral_losses_from_results <- function(results, feature_id, prepped, spectrum,
                                        tol = 0.01, include_precursor = TRUE) {
  row <- results[results$FEATURE_ID == feature_id, , drop = FALSE]
  if (nrow(row) == 0) return(NULL)
  clean <- function(x) gsub("_", " ", gsub("\\.", " ", x))
  colnames(row) <- clean(colnames(row))
  nm    <- clean(names(results))
  first <- which(nm == "peakCount filter") + 1
  last  <- which(nm == "score") - 1
  if (length(first) == 0 || length(last) == 0 || first > last) return(NULL)
  filter_cols <- nm[first:last]
  
  mz <- unlist(spectrum$mz); int <- unlist(spectrum$intensity)
  cand_mz <- mz; cand_int <- int
  if (include_precursor) {
    cand_mz  <- c(cand_mz, spectrum$precursorMz)
    cand_int <- c(cand_int, max(int))
  }
  mz_sorted <- sort(cand_mz)               
  keys <- vapply(prepped, `[[`, "", "key")
  
  out <- list()
  for (col in filter_cols) {
    r_idx <- match(col, keys)
    if (is.na(r_idx)) r_idx <- match(sub(" [0-9]+$", "", col), keys)
    if (is.na(r_idx)) next
    p <- prepped[[r_idx]]
    
    for (f in unlist(row[[col]][1])) {
      if (is.na(f)) next
      for (g in p$grids) {                         # largest repeat count first
        target <- f - g$sums
        ord <- order(target)
        idx <- integer(length(target))
        idx[ord] <- closest(target[ord], mz_sorted, tolerance = tol)
        hit_rows <- which(!is.na(idx))
        if (length(hit_rows) == 0) next
        partner <- mz_sorted[idx[hit_rows]]
        out[[length(out) + 1]] <- data.table::data.table(
          mz_found    = f,
          int_found   = cand_int[match(f, cand_mz)],
          mz_partner  = partner,
          int_partner = cand_int[match(partner, cand_mz)],
          delta       = g$sums[hit_rows],
          label       = p$label,
          masses      = lapply(hit_rows, function(r) p$searched[g$idx[, r]])
        )
        break                                      # first count with hits wins
      }
    }
  }
  if (length(out) == 0) return(NULL)
  as.data.frame(data.table::rbindlist(out))
}

#once per settings change......................................................
combos_with_repl <- function(n, k) {
  if (k == 1) return(matrix(seq_len(n), nrow = 1))
  utils::combn(n + k - 1, k) - (0:(k - 1))
}

prep_rule <- function(rule, mode) {
  masses <- vapply(rule$entries, function(e)
    resolve_mass(list(input_type = rule$input_type, value = e$value,
                      is_mzPeak = rule$is_mzPeak), mode), numeric(1))
  searched <- abs(masses[!is.na(masses)])
  if (!isTRUE(rule$sequential)) rule$min_repeat <- rule$max_repeat <- 1
  counts <- if (length(searched)) rev(seq.int(rule$min_repeat, rule$max_repeat)) else integer(0)
  list(
    label = as.character(rule$label)[1],
    key   = gsub("\\.", " ", make.names(as.character(rule$label)[1])),
    searched = searched,
    grids = lapply(counts, function(k) {
      idx <- combos_with_repl(length(searched), k)
      list(idx = idx, sums = colSums(matrix(searched[idx], nrow = k)))
    })
  )
}

#' @importFrom grDevices rainbow
#' @noRd
prepare_plot_context <- function(settings_list) {
  ld    <- settings_list$lipid_definition
  bf    <- settings_list$basicFilters
  mode  <- settings_list$project_info$mode
  chains <- ld$chains
  validations <- ld$validations
  if(length(chains) == 0)validations <- NULL
  
  chain_calcs <- lapply(chains, calc_lipid_chain_masses)
  names(chain_calcs) <- vapply(chains, function(x) as.character(x$chain_id)[1], character(1))
  
  headgroup_mass <- if (is.null(ld$head_group) || ld$head_group == "") 0
  else calculate_adjustment_mass(ld$head_group)
  
  validation_values <- lapply(validations, function(val) {
    ids    <- setdiff(val$components, "Backbone")
    combos <- do.call(expand.grid, lapply(ids, function(id) chain_calcs[[as.character(id)]]$mass))
    masses <- rowSums(combos) + calculate_adjustment_mass(val$adjustments) +
      if ("Backbone" %in% val$components) headgroup_mass else 0
    list(combos = combos, masses = masses)
  })
  
  levs <- sort(unique(c(vapply(settings_list$specific_filters, function(x) as.character(x$label), ""),
                        paste("chain", seq_along(chains)))))
  colors <- setNames(rainbow(length(levs)), levs)
  colors["none"] <- "#404040"
  colors["Intensity Filter"] <- "#D3D3D3"
  
  list(
    rules_prepped = lapply(settings_list$specific_filters, prep_rule, mode = mode),
    chain_calcs   = chain_calcs,
    validations   = validation_values,
    levs = levs, colors = colors,
    priorities = c(unlist(lapply(settings_list$specific_filters, function(x) x$rule_priority)),
                   unlist(lapply(chains, function(x) x$priority)),
                   unlist(lapply(validations, function(x) x$val_priority)),
                   "low"),                      # for the intensity filter
    tol = function(x) pmax(x * bf$tolerance_ppm / 1e6, bf$minimum_dalton)
  )
}

#overlay trace builders ........................................................
nl_label <- function(m) {
  tb <- table(m)
  if (sum(tb) > 1) paste(paste0(tb, "x ", names(tb)), collapse = ", ") else names(tb)
}

nl_specs <- function(nl, colors) {
  if (is.null(nl)) return(list())
  nl$y   <- pmax(nl$int_found, nl$int_partner)
  nl$col <- unname(colors[nl$label]); nl$col[is.na(nl$col)] <- "black"
  nl$txt <- vapply(nl$masses, nl_label, character(1))
  
  lines <- lapply(split(nl, nl$label), function(d) list(
    type = "scatter", mode = "lines",
    x = as.vector(rbind(d$mz_found, d$mz_partner, NA)),
    y = as.vector(rbind(d$y, d$y, NA)),
    line = list(color = d$col[1], width = 1.5, dash = "dot"),
    hoverinfo = "skip", showlegend = FALSE))
  txt <- list(list(
    type = "scatter", mode = "text",
    x = (nl$mz_found + nl$mz_partner) / 2, y = nl$y, text = nl$txt,
    textposition = "top center", textfont = list(size = 10, color = nl$col),
    hoverinfo = "none", showlegend = FALSE))
  c(unname(lines), txt)
}

points_spec <- function(p) {
  if (is.null(p) || nrow(p) == 0)
    return(list(type = "scatter", mode = "markers",
                x = numeric(0), y = numeric(0),
                hoverinfo = "skip", showlegend = FALSE))
  list(type = "scatter", mode = "markers+text",
       x = p$x, y = p$y * 1.05, text = p$text, textposition = "top center",
       marker = list(symbol = "triangle-down", size = 7, color = p$col),
       textfont = list(size = 10, color = p$col),
       hoverinfo = "text", showlegend = FALSE)
}

#' @title Chain Points
#'
#' @description Compares the chain values found in the results with the
#' corresponding chain masses and identifies the formula
#' 
#' @param results The results of the analysis.
#' @param mz mz values of one spectra
#' @param intensities intensities of the mz values
#' @param ctx containing relevant data for plotting like calculated chains, colors, priorities, etc.
#' 
#' @return A dataframe containing the x and y coordinates, the chains formula and color
#' 
#' @noRd
chain_points <- function(results, mz, intensities, ctx) {
  if (!length(mz)) return(NULL)
  
  cols <- grep("chain", names(results))
  pts <- list()
  
  for (k in seq_len(min(length(cols), length(ctx$chain_calcs)))) {
    found <- as.numeric(unlist(results[[cols[k]]])) #found values from results
    found <- unique(found[!is.na(found)])
    if (!length(found)) next
    
    ref  <- ctx$chain_calcs[[k]]
    near <- abs(outer(ref$mass, found, "-")) <=
      matrix(ctx$tol(found), nrow(ref), length(found), byrow = TRUE)
    keep <- which(colSums(near) > 0)
    if (!length(keep)) next
    
    fk  <- found[keep]
    pos <- vapply(fk, function(f) which.min(abs(mz - f)), integer(1))
    ok  <- abs(mz[pos] - fk) <= ctx$tol(fk)
    if (!any(ok)) next
    keep <- keep[ok]; pos <- pos[ok]
    
    col <- unname(ctx$colors[gsub("\\.", " ", names(results)[cols[k]])])
    pts[[length(pts) + 1]] <- data.frame(
      x = mz[pos], y = intensities[pos],
      text = vapply(keep, function(j) paste(ref$formula[near[, j]], collapse = " / "), ""),
      col = if (is.na(col)) "black" else col, stringsAsFactors = FALSE)
  }
  if (length(pts)) do.call(rbind, pts) else NULL
}

#' @title Validation Points
#'
#' @description Compares the validation values found in the results with the
#' corresponding validation masses and identifies the chains from which they
#' were derived.
#' 
#' @param results The results of the analysis.
#' @param mz mz values of one spectra
#' @param intensities intensities of the mz values
#' @param ctx containing relevant data for plotting like calculated chains, colors, priorities, etc.
#' 
#' @return A dataframe containing the x and y coordinates, the chain combination text and color
#' 
#' @noRd
validation_points <- function(results, mz, intensities, ctx) {
  n_val <- length(ctx$validations)
  if (n_val == 0) return(NULL)
  
  rn <- names(results)
  val_cols <- rn[(length(rn) - n_val):(length(rn) - 1)]
  pts <- list()
  
  #find the found chains for each validation
  for (k in seq_len(n_val)) {
    found <- as.numeric(unlist(results[[val_cols[k]]])) #found validation values from results
    found <- unique(found[!is.na(found)])
    if (!length(found)) next
    
    validation <- ctx$validations[[k]]
    
    #comparison of validation masses and found result values within tolerance
    near <- abs(outer(validation$masses, found, "-")) <=
      matrix(ctx$tol(found), length(validation$masses), length(found), byrow = TRUE)
    keep <- which(colSums(near) > 0)
    if (!length(keep)) next
    
    #Getting combination of chains for the validation
    pts[[length(pts) + 1]] <- data.frame(
      x = found[keep], y = intensities[match(found[keep], mz)],
      text = vapply(keep, function(j)
        paste0(round(as.numeric(unlist(validation$combos[which(near[, j])[1], ])),2), collapse = "+"), ""),
      col = "#404040", stringsAsFactors = FALSE)
  }
  if (length(pts)) do.call(rbind, pts) else NULL
}


#' @title Plotting one spectrum
#'
#' @description Creating a plotly plot with colors for the rule filters and 
#' neutral loss visualization.
#' 
#' @param sp A Spectra object with all spectra.
#' @param settings_list list of all setting parameters.
#' @param results The results of the analysis.
#' 
#' @return The plotly figure.
#' 
#' @importFrom plotly plot_ly add_segments add_trace layout
#' @importFrom grDevices rainbow
#' @noRd
# main plot ...................................................................
plottingSpectra <- function(spectrum, set_isotope_sum, settings_list, results,
                            ctx, show = character(0)) {
  mz <- unlist(spectrum$mz)
  intensities <- unlist(spectrum$intensity)
  filtered_spectrum <- cleaning_ref(spectrum, settings_list$basicFilters)
  if (set_isotope_sum) results <- summarize_isotope_peaks(results, spectrum)
  
  start_filtercols <- which(colnames(results) == "peakCount_filter") + 1
  filter_results <- results[, start_filtercols:(ncol(results) - 1), drop = FALSE]
  filter_results[["Intensity Filter"]] <- list(setdiff(mz, unlist(filtered_spectrum$mz)))
  
  category <- resolve_filter_category(mz, filter_results, ctx$priorities)
  category <- factor(category, levels = c("none", "Intensity Filter", ctx$levs))
  nl <- neutral_losses_from_results(results, results$FEATURE_ID[1], ctx$rules_prepped,
                                    spectrum, tol = settings_list$basicFilters$minimum_dalton)
  
  overlays <- list(
    nl    = nl_specs(nl, ctx$colors),
    chain = list(points_spec(chain_points(results, mz, intensities, ctx))),
    val   = list(points_spec(validation_points(results, mz, intensities, ctx)))
  )

  # overlay traces come first -> their indices are known (0-based, as plotly.js expects)
  fig <- plot_ly()
  idx <- list(); n <- 0L
  for (grp in names(overlays)) {
    specs <- overlays[[grp]]
    idx[[grp]] <- n + seq_along(specs) - 1L
    for (s in specs)
      fig <- do.call(plotly::add_trace, c(list(p = fig), s, list(visible = grp %in% show)))
    n <- n + length(specs)
  }

  fig <- fig %>%
    add_segments(x = spectrum$precursorMz, xend = spectrum$precursorMz, 
                 y = 0, yend = max(intensities),
                 line = list(color = "red", dash = "dash", width = 1),
                 name = "Precursor", inherit = FALSE)
  
  hover <- paste0("m/z: ", mz, "<br>Intensity: ", intensities)
  for (lv in levels(category)) {
    i <- which(category == lv)
    if(length(mz[i])==0) next
    
    if (!length(i)) next
    col <- unname(ctx$colors[lv])

    fig <- fig %>% add_segments(
      x = mz[i], xend = mz[i], y = 0, yend = intensities[i],
      name = lv, legendgroup = lv,
      line = list(color = col, width = 1),
      text = hover[i], hoverinfo = "text", showlegend = TRUE, inherit = FALSE
    )
  }

  fig <- fig %>% layout(hovermode = "x", xaxis = list(title = "m/z"),
                        yaxis = list(title = "Intensity"))
  
  list(fig = fig, idx = idx)
}




#' @title Ordering the filter results by priority
#'
#' @description For each mz peak in the spectrum the results are compared and
#' from the fitting ones the filter name with the highest priority is taken.
#' 
#' @details
#' Creates a Tuple of the name of the filters and their priority. Then goes through
#' each mz value and compares the masses of the filter results to the mz value.
#' If a filter has a match to the mz value, it is counted as a hit. At the end 
#' the filter with the highest priority is selected for the mz value. After all
#' mz values have a filter, their names are returned as category.
#' 
#' 
#' @param mz The mz values of the spectrum
#' @param filter_results results of the analysis
#' @param filter_priority priority of the result filters
#' 
#' @return The names of the prioritized rules for each mz value.
#' 
#' @noRd
resolve_filter_category <- function(mz, filter_results, filter_priority) {
  priority_levels <- c("low", "basic", "important", "required")
  
  keep_cols <- !sapply(filter_results, function(x) all(is.na(x)))
  filter_results <- filter_results[, keep_cols, drop = FALSE]
  filter_priority <- filter_priority[keep_cols]

  rank <- match(filter_priority, priority_levels)
  rank[is.na(rank)] <- 1L
  
  hit <- vapply(filter_results,
                function(vals) mz %in% unlist(vals, use.names = FALSE),
                logical(length(mz)))
  if (!is.matrix(hit)) hit <- matrix(hit, nrow = length(mz))
  
  best <- max.col(sweep(hit, 2, rank, `*`), ties.method = "first")
  out  <- gsub("\\.", " ", names(filter_results)[best])
  out[rowSums(hit) == 0] <- "none"
  out
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
  
  if("MaxIntensity" %in% colnames(ms1_data)){
    samples <- as.data.frame(ms1_data[, (which(colnames(ms1_data) == "MaxIntensity")+1):ncol(ms1_data)])
  } else  samples <- as.data.frame(ms1_data[, 1:(ncol(ms1_data)-1)])
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

collapse_peaks <- function(x, mz, int, tol = 1.1) {
  if (all(is.na(x))) return(NA_real_)   
  x <- sort(unique(x[!is.na(x)]))
  if (length(x) < 2) return(x)
  
  grp <- cumsum(c(TRUE, diff(x) >= tol))
  
  unname(vapply(split(x, grp), function(g) {
    idx <- vapply(g, function(v) which.min(abs(mz - v)), integer(1))
    g[which.max(int[idx])]
  }, numeric(1)))
}

summarize_isotope_peaks <- function(results, spectrum) {
  mz  <- unlist(spectrum$mz)
  int <- unlist(spectrum$intensity)
  
  nm    <- gsub("\\.", " ", names(results))
  first <- which(nm == "peakCount_filter") + 1
  last  <- which(nm == "score") - 1
  val_cols <- nm[first:last]
  
  results_collapsed <- results
  for (j in seq_along(results)) {
    if (nm[j] %in% val_cols) {
      results_collapsed[[j]] <- lapply(results[[j]], collapse_peaks,
                                       mz = mz, int = int, tol = 1.1)
    }
  }
  results_collapsed
}
