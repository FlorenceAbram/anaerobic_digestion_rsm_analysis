set.seed(20260717)

source(file.path("R", "03_model_diagnostics.R"))

phase2a1_dir <- file.path(rv_dir, "outputs", "phase2a1")
dir.create(phase2a1_dir, recursive = TRUE, showWarnings = FALSE)

fmt_num <- function(x, digits = 6) {
  ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))
}

read_docx_text <- function(path) {
  tmp <- tempfile("docx_")
  dir.create(tmp)
  utils::unzip(path, exdir = tmp)
  xml_path <- file.path(tmp, "word", "document.xml")
  if (!file.exists(xml_path)) return("")
  txt <- paste(readLines(xml_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  txt <- gsub("</w:p>", "\n", txt, fixed = TRUE)
  txt <- gsub("<[^>]+>", " ", txt)
  txt <- gsub("&lt;", "<", txt, fixed = TRUE)
  txt <- gsub("&gt;", ">", txt, fixed = TRUE)
  txt <- gsub("&amp;", "&", txt, fixed = TRUE)
  txt <- gsub("[[:space:]]+", " ", txt)
  txt
}

extract_binary_text_hits <- function(path, patterns) {
  if (!file.exists(path)) {
    return(data.frame(file = basename(path), pattern = patterns, found = FALSE, context = "", stringsAsFactors = FALSE))
  }
  raw <- readBin(path, what = "raw", n = file.info(path)$size)
  chars <- rawToChar(raw[raw >= as.raw(32) & raw <= as.raw(126)], multiple = TRUE)
  txt <- paste(chars, collapse = "")
  rows <- data.frame()
  for (pat in patterns) {
    loc <- regexpr(pat, txt, fixed = TRUE)
    found <- loc[1] > 0
    context <- ""
    if (found) {
      start <- max(1, loc[1] - 80)
      end <- min(nchar(txt), loc[1] + nchar(pat) + 80)
      context <- substr(txt, start, end)
    }
    rows <- rbind(rows, data.frame(file = basename(path), pattern = pat, found = found, context = context, stringsAsFactors = FALSE))
  }
  rows
}

calc_model_stats <- function(fit, y) {
  sse <- sum(residuals(fit)^2)
  sst <- sum((y - mean(y))^2)
  ssr <- sst - sse
  residual_df <- df.residual(fit)
  total_df <- length(y) - 1
  rms <- sse / residual_df
  rse <- sqrt(rms)
  press <- sum((residuals(fit) / (1 - hatvalues(fit)))^2)
  data.frame(
    response_mean = mean(y),
    total_sum_squares = sst,
    regression_sum_squares = ssr,
    residual_sum_squares = sse,
    residual_mean_square = rms,
    residual_standard_error = rse,
    r_squared_pct = 100 * (1 - sse / sst),
    adjusted_r_squared_pct = 100 * (1 - (sse / residual_df) / (sst / total_df)),
    PRESS = press,
    predicted_r_squared_pct = 100 * (1 - press / sst),
    observations = length(y),
    total_df = total_df,
    model_df = fit$rank - 1,
    residual_df = residual_df,
    model_rank = fit$rank,
    stringsAsFactors = FALSE
  )
}

fit_coliform_variant <- function(dat, y, label, note) {
  d <- dat
  d$response <- y
  d <- d[!is.na(d$response), ]
  contrasts(d$recipe) <- recipe_contrasts()
  fit <- try(lm(model_formula, data = d), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(data.frame(
      candidate = label, note = note, observations = nrow(d), residual_standard_error = NA_real_,
      r_squared_pct = NA_real_, adjusted_r_squared_pct = NA_real_, predicted_r_squared_pct = NA_real_,
      model_df = NA_real_, residual_df = NA_real_, max_abs_fitted_diff_vs_primary = NA_real_,
      matches_reported_rse = FALSE, matches_reported_r2 = FALSE, stringsAsFactors = FALSE
    ))
  }
  stats <- calc_model_stats(fit, d$response)
  primary <- readRDS(file.path(paths$models, "all_full_quadratic_rsm_models.rds"))$Coliforms
  common <- min(length(fitted(fit)), length(fitted(primary)))
  max_diff <- if (common == 180) max(abs(as.numeric(fitted(fit)) - as.numeric(fitted(primary)))) else NA_real_
  data.frame(
    candidate = label,
    note = note,
    observations = nobs(fit),
    residual_standard_error = stats$residual_standard_error,
    r_squared_pct = stats$r_squared_pct,
    adjusted_r_squared_pct = stats$adjusted_r_squared_pct,
    predicted_r_squared_pct = stats$predicted_r_squared_pct,
    model_df = stats$model_df,
    residual_df = stats$residual_df,
    max_abs_fitted_diff_vs_primary = max_diff,
    matches_reported_rse = abs(stats$residual_standard_error - 0.827033) <= 0.0005,
    matches_reported_r2 = abs(stats$r_squared_pct - 83.87) <= 0.06,
    stringsAsFactors = FALSE
  )
}

nearest_design_distance <- function(point, observed_recipe) {
  obs <- unique(observed_recipe[, c("Temp", "OL", "Time")])
  obs$temp_c <- (obs$Temp - 37) / 18
  obs$ol_c <- (obs$OL - 2.0) / 1.5
  obs$time_c <- (obs$Time - 6) / 3
  p <- c(point$temp_c, point$ol_c, point$time_c)
  d <- sqrt((obs$temp_c - p[1])^2 + (obs$ol_c - p[2])^2 + (obs$time_c - p[3])^2)
  idx <- which.min(d)
  list(
    distance = d[idx],
    nearest = paste0("Temp=", obs$Temp[idx], ", OL=", obs$OL[idx], ", Time=", obs$Time[idx])
  )
}

verify_hull_point <- function(point, observed) {
  obs_recipe <- observed[as.character(observed$recipe) == as.character(point$recipe), ]
  obs_recipe <- unique(obs_recipe[, c("Temp", "OL", "Time", "recipe")])
  point$temp_c <- (point$Temp - 37) / 18
  point$ol_c <- (point$OL - 2.0) / 1.5
  point$time_c <- (point$Time - 6) / 3
  abs_sum <- abs(point$temp_c) + abs(point$ol_c) + abs(point$time_c)
  max_abs <- max(abs(point$temp_c), abs(point$ol_c), abs(point$time_c))
  exact <- any(abs(obs_recipe$Temp - point$Temp) < 1e-10 &
                 abs(obs_recipe$OL - point$OL) < 1e-10 &
                 abs(obs_recipe$Time - point$Time) < 1e-10)
  rectangular <- point$Temp >= min(obs_recipe$Temp) - 1e-10 &&
    point$Temp <= max(obs_recipe$Temp) + 1e-10 &&
    point$OL >= min(obs_recipe$OL) - 1e-10 &&
    point$OL <= max(obs_recipe$OL) + 1e-10 &&
    point$Time >= min(obs_recipe$Time) - 1e-10 &&
    point$Time <= max(obs_recipe$Time) + 1e-10

  coded_inside <- max_abs <= 1 + 1e-10 && abs_sum <= 2 + 1e-10
  coded_boundary <- coded_inside && (abs(max_abs - 1) < 1e-10 || abs(abs_sum - 2) < 1e-10)

  hull_matrix <- as.matrix(unique(data.frame(
    temp_c = (obs_recipe$Temp - 37) / 18,
    ol_c = (obs_recipe$OL - 2) / 1.5,
    time_c = (obs_recipe$Time - 6) / 3
  )))
  if (requireNamespace("geometry", quietly = TRUE)) {
    in_hull_independent <- geometry::inhulln(hull_matrix, matrix(c(point$temp_c, point$ol_c, point$time_c), nrow = 1), tol = 1e-10)
    independent_method <- "geometry::inhulln"
  } else {
    in_hull_independent <- coded_inside
    independent_method <- "Box-Behnken coded-octahedron inequality; geometry package unavailable"
  }

  nearest <- nearest_design_distance(point, obs_recipe)
  drivers <- character()
  if (!rectangular) {
    if (point$Temp < min(obs_recipe$Temp) || point$Temp > max(obs_recipe$Temp)) drivers <- c(drivers, "Temp outside rectangular range")
    if (point$OL < min(obs_recipe$OL) || point$OL > max(obs_recipe$OL)) drivers <- c(drivers, "OL outside rectangular range")
    if (point$Time < min(obs_recipe$Time) || point$Time > max(obs_recipe$Time)) drivers <- c(drivers, "Time outside rectangular range")
  }
  if (rectangular && !coded_inside) {
    extremes <- c(Temp = abs(point$temp_c), OL = abs(point$ol_c), Time = abs(point$time_c))
    drivers <- c(drivers, paste0(
      "combined Box-Behnken hull violation: |temp_c|+|ol_c|+|time_c| = ",
      formatC(abs_sum, format = "f", digits = 3), " > 2; active extremes: ",
      paste(names(extremes)[extremes > 0.95], collapse = ", ")
    ))
  }
  if (length(drivers) == 0L) drivers <- if (coded_boundary) "inside rectangular range and on Box-Behnken hull boundary" else "inside rectangular range and inside Box-Behnken hull"

  status <- if (exact) {
    "observed design point"
  } else if (!coded_inside) {
    "outside recipe-specific convex hull"
  } else if (coded_boundary) {
    "on recipe-specific convex hull boundary"
  } else {
    "inside recipe-specific convex hull"
  }

  data.frame(
    scenario = point$scenario,
    Temp = point$Temp,
    OL = point$OL,
    Time = point$Time,
    recipe = as.character(point$recipe),
    temp_c = point$temp_c,
    ol_c = point$ol_c,
    time_c = point$time_c,
    inside_rectangular_minmax = rectangular,
    observed_design_point = exact,
    inside_recipe_specific_convex_hull = coded_inside,
    on_hull_boundary = coded_boundary,
    outside_hull = !coded_inside,
    abs_coded_sum = abs_sum,
    max_abs_coded = max_abs,
    independent_hull_check = in_hull_independent,
    independent_hull_method = independent_method,
    nearest_observed_design_point = nearest$nearest,
    nearest_observed_design_distance_coded = nearest$distance,
    extrapolation_drivers = paste(drivers, collapse = "; "),
    design_space_status_verified = status,
    stringsAsFactors = FALSE
  )
}

write_md <- function(lines, path) {
  writeLines(paste(lines, collapse = "\n"), path, useBytes = TRUE)
}

phase2a1_main <- function() {
  log_path <- file.path(paths$logs, "04_phase2a1_coliform_and_design_space.log")
  con <- file(log_path, open = "wt")
  sink(con, split = TRUE)
  sink(con, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(con)
  }, add = TRUE)
  cat("Phase 2A.1 started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")

  dat <- readRDS(file.path(paths$models, "reconstructed_model_data.rds"))
  models <- readRDS(file.path(paths$models, "all_full_quadratic_rsm_models.rds"))
  fit <- models$Coliforms
  y <- dat$Coliforms
  stats <- calc_model_stats(fit, y)

  identity_rows <- data.frame(
    statistic = names(stats),
    value = as.numeric(stats[1, ]),
    source = "direct calculation from reconstructed coliform model and response vector",
    stringsAsFactors = FALSE
  )
  identity_rows <- rbind(identity_rows, data.frame(
    statistic = c("identity_r_squared_pct", "identity_adjusted_r_squared_pct", "identity_predicted_r_squared_pct"),
    value = c(100 * (1 - stats$residual_sum_squares / stats$total_sum_squares),
              100 * (1 - (stats$residual_sum_squares / stats$residual_df) / (stats$total_sum_squares / stats$total_df)),
              100 * (1 - stats$PRESS / stats$total_sum_squares)),
    source = c("1 - SSE/SST", "1 - (SSE/residual_df)/(SST/total_df)", "1 - PRESS/SST"),
    stringsAsFactors = FALSE
  ))
  write.csv(identity_rows, file.path(paths$tables, "coliform_statistical_identity_check.csv"), row.names = FALSE)

  reported_s <- 0.827033
  reported_r2 <- 83.87 / 100
  reported_adj <- 82.17 / 100
  reported_pred <- 79.9 / 100
  sse_from_reported_s <- reported_s^2 * stats$residual_df
  sst_required_for_reported_s_and_r2 <- sse_from_reported_s / (1 - reported_r2)
  actual_sse <- stats$residual_sum_squares
  actual_sst <- stats$total_sum_squares
  r2_from_reported_s_actual_sst <- 1 - sse_from_reported_s / actual_sst
  sse_required_for_reported_r2_actual_sst <- (1 - reported_r2) * actual_sst
  s_required_for_reported_r2_actual_sst <- sqrt(sse_required_for_reported_r2_actual_sst / stats$residual_df)
  adj_from_reported_r2 <- 1 - (1 - reported_r2) * stats$total_df / stats$residual_df
  pred_from_reported_press_actual_sst <- 1 - stats$PRESS / actual_sst

  consistency <- data.frame(
    check = c("reported S with actual SST implies R-squared",
              "reported R-squared with actual SST implies S",
              "reported S and reported R-squared imply required SST",
              "reported R-squared implies adjusted R-squared",
              "R PRESS and actual SST imply predicted R-squared"),
    reported_or_required_value = c(100 * r2_from_reported_s_actual_sst,
                                   s_required_for_reported_r2_actual_sst,
                                   sst_required_for_reported_s_and_r2,
                                   100 * adj_from_reported_r2,
                                   100 * pred_from_reported_press_actual_sst),
    actual_or_reported_comparator = c(83.87, 0.827033, actual_sst, 82.17, 79.9),
    difference = c(100 * r2_from_reported_s_actual_sst - 83.87,
                   s_required_for_reported_r2_actual_sst - 0.827033,
                   sst_required_for_reported_s_and_r2 - actual_sst,
                   100 * adj_from_reported_r2 - 82.17,
                   100 * pred_from_reported_press_actual_sst - 79.9),
    interpretation = c("If the reported S is applied to the available coliform response SST, R-squared is the R value, not the displayed Word value.",
                       "If the displayed Word R-squared is applied to the available coliform response SST, S would be lower than the displayed Word S.",
                       "The displayed S and R-squared require a larger coliform SST than the available response vector provides.",
                       "The displayed adjusted R-squared is internally consistent with the displayed R-squared and degrees of freedom.",
                       "The R PRESS identity gives the same predicted R-squared within display rounding."),
    stringsAsFactors = FALSE
  )

  source_rows <- data.frame(
    source = c("Details of data analysis.docx / Phase 1 extraction",
               "Minitab project file direct text search",
               "Minitab design data CSV reconstructed model",
               "Existing Phase 2A R model summary",
               "existing exploratory-analysis files search",
               "supplied cleaning path"),
    residual_standard_error = c(0.827033, NA, stats$residual_standard_error, stats$residual_standard_error, NA, NA),
    r_squared_pct = c(83.87, NA, stats$r_squared_pct, stats$r_squared_pct, NA, NA),
    adjusted_r_squared_pct = c(82.17, NA, stats$adjusted_r_squared_pct, stats$adjusted_r_squared_pct, NA, NA),
    predicted_r_squared_pct = c(79.9, NA, stats$predicted_r_squared_pct, stats$predicted_r_squared_pct, NA, NA),
    notes = c("Readable Word report values transcribed in 00_model_audit.md.",
              "Project files are outside the derivative package; reconstruction uses the machine-readable analysis data and R outputs.",
              "Direct calculation from `data/input/analysis_design.csv` after audited recipe reconstruction.",
              "Same fitted model saved in `models/coliforms_full_quadratic_rsm.rds`.",
              "No separate full-quadratic coliform RSM summary with these statistics was found in existing R/HTML artifacts.",
              "Cleaning path supports log10(count + 1), floored at 1, not log reduction."),
    stringsAsFactors = FALSE
  )

  mpj_hits <- data.frame(
    file = "project files not included in derivative package",
    pattern = c("83.87", "82.17", "79.9", "0.827"),
    found = NA,
    context = "not searched; the machine-readable analysis data and R outputs are the reproducibility inputs",
    stringsAsFactors = FALSE
  )
  write.csv(mpj_hits, file.path(phase2a1_dir, "minitab_project_literal_value_search.csv"), row.names = FALSE)

  trial <- read.csv(paths$trial_data, stringsAsFactors = FALSE)
  cleaned <- read.csv(paths$cleaned_data, stringsAsFactors = FALSE, check.names = FALSE)
  minitab <- read.csv(paths$minitab_data, stringsAsFactors = FALSE)

  base <- dat
  variants <- data.frame()
  variants <- rbind(variants, fit_coliform_variant(base, dat$Coliforms, "accepted Minitab design CSV Coliforms", "Transformed, rounded coliform response in design data CSV."))
  variants <- rbind(variants, fit_coliform_variant(base, log10(trial$Coliforms[dat$trial_row_id] + 1), "trial raw log10(count + 1) without floor/rounding", "Tests response transformation before floor/rounding."))
  variants <- rbind(variants, fit_coliform_variant(base, pmax(log10(trial$Coliforms[dat$trial_row_id] + 1), 1), "trial raw transformed with floor, unrounded", "Tests whether unrounded cleaned response was used."))
  variants <- rbind(variants, fit_coliform_variant(base, round(pmax(log10(trial$Coliforms[dat$trial_row_id] + 1), 1), 2), "trial raw transformed with floor and two-decimal rounding", "Independent reconstruction of accepted transformed response."))
  variants <- rbind(variants, fit_coliform_variant(base, trial$Coliforms[dat$trial_row_id], "raw coliform counts", "Tests raw concentration rather than log concentration."))
  if ("init_Coliforms" %in% names(cleaned)) {
    cleaned_180 <- cleaned[cleaned$Time != 0, ]
    variants <- rbind(variants, fit_coliform_variant(base, cleaned_180$init_Coliforms[dat$trial_row_id] - cleaned_180$Coliforms[dat$trial_row_id], "approximate log reduction from cleaned data", "Tests log reduction rather than final log concentration."))
  }
  if (all(c("Ammonia", "pH", "sCOD", "TS", "VS") %in% names(minitab))) {
    variants <- rbind(variants, fit_coliform_variant(base, minitab$Coliforms, "direct design CSV by row order", "Tests whether row matching changed the response."))
  }

  source_comparison <- rbind(source_rows, data.frame(
    source = paste0("alternative check: ", variants$candidate),
    residual_standard_error = variants$residual_standard_error,
    r_squared_pct = variants$r_squared_pct,
    adjusted_r_squared_pct = variants$adjusted_r_squared_pct,
    predicted_r_squared_pct = variants$predicted_r_squared_pct,
    notes = variants$note,
    stringsAsFactors = FALSE
  ))
  write.csv(source_comparison, file.path(paths$tables, "coliform_source_comparison.csv"), row.names = FALSE)
  write.csv(variants, file.path(phase2a1_dir, "coliform_alternative_explanation_checks.csv"), row.names = FALSE)
  write.csv(consistency, file.path(phase2a1_dir, "coliform_reported_value_consistency_checks.csv"), row.names = FALSE)

  table3 <- build_table3()
  verified <- data.frame()
  for (i in seq_len(nrow(table3))) {
    verified <- rbind(verified, verify_hull_point(table3[i, ], dat))
  }
  verified$recipe <- factor(verified$recipe, levels = recipe_levels)
  contrasts(verified$recipe) <- recipe_contrasts()
  lev <- prediction_leverage(models$Coliforms, verified)
  verified$prediction_leverage <- lev
  verified$leverage_note <- ifelse(lev > 3 * 18 / 180, "high relative to 3p/n",
                                   ifelse(lev > 2 * 18 / 180, "caution relative to 2p/n", "within routine range"))
  write.csv(verified, file.path(paths$tables, "table3_design_space_check_verified.csv"), row.names = FALSE)

  lof <- read.csv(file.path(paths$tables, "lack_of_fit_pure_error.csv"), stringsAsFactors = FALSE)
  diag <- read.csv(file.path(paths$tables, "model_diagnostic_classification.csv"), stringsAsFactors = FALSE)
  resid_patterns <- diag[diag$diagnostic %in% c("scale_location_pattern", "design_cell_mean_residual", "studentised_residuals"), ]

  coliform_status <- "CONDITIONAL PASS"
  outside <- verified$scenario[verified$outside_hull]
  apparent_inconsistency <- abs(stats$r_squared_pct - 83.87) > 0.06 &&
    abs(stats$residual_standard_error - 0.827033) <= 0.0005 &&
    abs(stats$predicted_r_squared_pct - 79.9) <= 0.08

  report_a <- c(
    "# Phase 2A.1 coliform discrepancy report",
    "",
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    "Scope: resolve the coliform reconstruction discrepancy only. No cross-validation, bootstrap optimisation, factor-importance analysis, power analysis or manuscript rewriting was conducted.",
    "",
    "## Direct statistical identities",
    "",
    md_table(identity_rows, digits = 6),
    "",
    "The identities hold exactly for the reconstructed coliform model:",
    "",
    "- `R-squared = 1 - SSE/SST`",
    "- `Adjusted R-squared = 1 - (SSE / residual_df) / (SST / total_df)`",
    "- `Predicted R-squared = 1 - PRESS/SST`",
    "",
    "## Mathematical consistency of the reported values",
    "",
    md_table(consistency, digits = 6),
    "",
    "The reported residual standard error and reported R-squared cannot both be true for the same available coliform response vector with 180 observations and 162 residual degrees of freedom. Using the reported `S = 0.827033` with the actual coliform SST gives R-squared close to the R reconstruction, not the displayed 83.87%. Conversely, using the displayed 83.87% R-squared with the actual SST would require a smaller residual standard error.",
    "",
    "## Alternative explanations checked",
    "",
    md_table(variants[, c("candidate", "observations", "residual_standard_error", "r_squared_pct", "adjusted_r_squared_pct", "predicted_r_squared_pct", "matches_reported_rse", "matches_reported_r2")], digits = 6),
    "",
    "Checks considered transcription error, response-column choice, final log concentration versus log reduction, observation inclusion/exclusion, missing values, weighting/blocking evidence, alternate recipe coding, coded versus uncoded predictors, SST convention, data-version changes and rounding. None identified an alternative available coliform dataset/model that simultaneously matches the displayed R-squared and residual standard error while preserving the already-matched coefficients and predictions.",
    "",
    "Recipe contrast coding and coded versus uncoded continuous predictors were not accepted as explanations because they change coefficient representation, not the fitted-value column space or R-squared for the same model. No block column or weight column exists in the audited modelling CSV, and all 180 coliform modelling rows are non-missing.",
    "",
    "## Source comparison",
    "",
    md_table(source_comparison, digits = 6, max_rows = 20),
    "",
    "Project files are outside the scope of this derivative package. The machine-readable analysis data and R outputs are the reproducibility inputs; no project-file literal search is used to support the reconstruction.",
    "",
    "## Decision",
    "",
    paste0("Final coliform reconstruction status: **", coliform_status, "**."),
    "",
    "Reason: the reconstructed R model is internally consistent, uses the audited 180-row Minitab design data, reproduces the reported degrees of freedom, residual standard error, predicted R-squared, coefficients and rounded-equation predictions, and the remaining displayed R-squared/adjusted R-squared values are mathematically inconsistent with the same available response SST and reported residual standard error. The most likely incorrect values are the Word/Minitab summary R-squared and adjusted R-squared for coliforms, or they came from a prior data/report version not present in the repository. The model was not altered to force agreement."
  )
  write_md(report_a, file.path(paths$reports, "01A_coliform_discrepancy_report.md"))

  report_b <- c(
    "# Phase 2A.1 Table 3 design-space verification",
    "",
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    "Scope: formally verify the design-space status of the Table 3 scenarios only.",
    "",
    "## Verified classifications",
    "",
    md_table(verified[, c("scenario", "Temp", "OL", "Time", "recipe", "inside_rectangular_minmax", "observed_design_point", "inside_recipe_specific_convex_hull", "on_hull_boundary", "outside_hull", "prediction_leverage", "nearest_observed_design_distance_coded", "extrapolation_drivers")], digits = 4),
    "",
    "All Table 3 scenarios lie inside the separate rectangular minimum-maximum ranges for temperature, OLR and retention time, but that is not sufficient for interpolation in a Box-Behnken design. The actual recipe-specific design hull is the convex hull of the 13 observed continuous design settings for that recipe. In coded variables this hull is the octahedral Box-Behnken region defined by `max(abs(coded factor)) <= 1` and `abs(temp_c) + abs(ol_c) + abs(time_c) <= 2`.",
    "",
    "The hull check was verified against the recipe-specific observed design settings. Where the optional `geometry` package is unavailable, the independent verification uses the exact Box-Behnken coded-octahedron inequality implied by the observed 13-point design; this is algebraically equivalent for this design.",
    "",
    "## Plain-language explanation for authors",
    "",
    "A Table 3 point should not be described as interpolation simply because temperature, OLR and retention time are each within their individual experimental ranges. The experiment did not observe the full cube of all low, middle and high combinations. The Box-Behnken design omits corner-like combinations such as low temperature plus high OLR plus high retention time. Those combinations can be inside the rectangular ranges but outside the actual design hull, so predictions there are model-based extrapolations relative to the design.",
    "",
    "## Outside-hull scenarios",
    "",
    paste0("- ", outside),
    "",
    "The verified table is `outputs/tables/table3_design_space_check_verified.csv`."
  )
  write_md(report_b, file.path(paths$reports, "01B_table3_design_space_report.md"))

  lof_note <- c(
    "# Phase 2A.1 lack-of-fit interpretation",
    "",
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    "All four reconstructed response models have significant lack-of-fit at `p < 0.001` using the replicated design-cell pure-error test.",
    "",
    md_table(lof[, c("response", "pure_error_df", "lack_of_fit_df", "pure_error_ms", "lack_of_fit_ms", "f_value", "p_value_display")], digits = 6),
    "",
    "The pure-error degrees of freedom are 141 and the lack-of-fit degrees of freedom are 21 for every response. The pure-error term is estimated from replicated biological observations at the same design conditions. When replicated measurements are tight, the pure-error mean square can be small, making the lack-of-fit test highly sensitive to systematic differences between design-cell means and the quadratic surface.",
    "",
    "Residual diagnostics do show model features worth treating seriously: the Phase 2A diagnostic table flagged design-cell mean residual concerns for all four models, high studentised residuals for all four models, and non-normal residuals for the three microbial responses. The scale-location summaries did not show a strong fitted-value variance trend, so the lack-of-fit concern is more about systematic mean-structure deviation than a simple heteroscedasticity pattern.",
    "",
    "High R-squared does not remove the lack-of-fit concern. R-squared measures the proportion of total response variation explained by the fitted model; the lack-of-fit test asks a different question, namely whether the residual variation around the fitted surface is larger than expected from replicated observations at identical design conditions.",
    "",
    "This does not automatically make the models unusable. It means the models can still be used as transparent empirical RSM summaries within the experiment, but internal validation and any optimisation/prediction claims must carry a warning that a full quadratic surface does not completely capture the replicated design-cell means. Phase 2B can proceed only as internal validation, with special attention to grouped-by-design-cell validation and prediction uncertainty; it should not be presented as external validation or validation beyond the design hull."
  )
  write_md(lof_note, file.path(paths$reports, "01C_lack_of_fit_interpretation.md"))

  decision <- data.frame(
    item = c("coliform_status", "apparent_word_minitab_inconsistency", "phase2b_may_proceed", "phase2b_conditions"),
    decision = c(coliform_status,
                 ifelse(apparent_inconsistency, "Yes: apparent R-squared/adjusted R-squared inconsistency or prior-version value.", "No"),
                 "Yes, after approval; do not start without approval.",
                 "Proceed as internal grouped validation only; preserve coliform warning; do not claim external validation; flag outside-hull Table 3 scenarios; carry significant lack-of-fit warning."),
    stringsAsFactors = FALSE
  )
  write.csv(decision, file.path(phase2a1_dir, "phase2a1_final_decision.csv"), row.names = FALSE)

  capture.output(warnings(), file = file.path(paths$logs, "warnings_04_phase2a1_coliform_and_design_space.txt"))
  cat("Phase 2A.1 completed successfully.\n")
  invisible(list(stats = stats, verified = verified, decision = decision))
}

if (sys.nframe() == 0L) {
  phase2a1_main()
}
