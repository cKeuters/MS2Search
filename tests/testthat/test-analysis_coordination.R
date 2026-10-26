test_that("Loading MS2 data from mgf works", {
  sp <- load_spectra("Test_dataset.mgf")
  expect_equal("rtime" %in% spectraVariables(sp), TRUE)
})

test_that("Finding Neutral Losses works", {
  nl <- find_neutralLosses(peaks = c(10,20,30,33,33.1, 50, 60, 60.3,100,110,300,350,400,800,20000),
                     precursor = 1000,
                     values = c(10,100),
                     losses_min = 1,
                     losses_max = 2,
                     tol_ppm = 10,
                     tol_mz = 0.1,
                     rule_min_mz = 0,
                     rule_max_mz = 1000,
                     just_after_precursor = FALSE)
  

  expect_equal(nl, c(60,110,400))
})

test_that("Running complete Analysis", {
  sp <- load_spectra("Test_dataset.mgf")
  settings_list <- jsonlite::fromJSON("settings.json", simplifyVector = FALSE)
  write_log <- function(msg) {
    cat(format(Sys.time(), "[%H:%M:%S]"), msg, "\n")
  }

  res <- run_analysis(settings_list, sp, add_log = write_log)
  
  res_test <- data.frame(
    FEATURE_ID = 5715,
    mz_filter = TRUE,
    rt_filter = TRUE,
    peakCount_filter = TRUE,
    score = 0.8,
    check.names = FALSE
  )

  res_test[["Water Losses"]] <- list(c(
    300.287959, 318.300211, 638.608324,
    639.612216, 656.618718, 657.621928
  ))
  res_test[["Glucose Loss"]] <- list(c(
    800.658054, 801.665675, 818.671180
  ))
  res_test[["chain 1"]] <- list(c(318.300211))
  res_test[["chain 2"]] <- list(c(283.2843429, 300.2879595, 356.3508931))
  res_test[["precursor validation"]] <- list(NA)

  res_test <- res_test[, names(res)]

  expect_equal(res, res_test)
})