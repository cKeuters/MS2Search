
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `{MS2Search}`

<!-- badges: start -->

<!-- badges: end -->

## Installation

You can install the development version of `{MS2Search}` like so:

``` r
# FILL THIS IN! HOW CAN PEOPLE INSTALL YOUR DEV PACKAGE?
```

## Run

You can launch the application by running:

``` r
MS2Search::run_app()
```

## About

You are reading the doc about version : 0.0.0.9000

This README has been compiled on the

``` r
Sys.time()
#> [1] "2026-07-30 10:01:01 UTC"
```

Here are the tests results and package coverage:

``` r
devtools::check(quiet = TRUE)
#> ══ Documenting ═════════════════════════════════════════════════════════════════
#> ℹ Installed roxygen2 version (8.0.0) doesn't match declared (7.1.1)
#> ✖ `check()` will not re-document this package.
#> ℹ Do you need to re-run `document()`?
#> ── R CMD check results ─────────────────────────────── MS2Search 0.0.0.9000 ────
#> Duration: 1m 15.3s
#> 
#> ❯ checking top-level files ... NOTE
#>   Non-standard file/directory found at top level:
#>     ‘old_files’
#> 
#> 0 errors ✔ | 0 warnings ✔ | 1 note ✖
```

``` r
covr::package_coverage()
#> Error in `loadNamespace()`:
#> ! there is no package called 'covr'
```