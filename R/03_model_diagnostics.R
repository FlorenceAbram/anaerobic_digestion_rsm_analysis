set.seed(20260717)

source(file.path("R", "02_anova_lack_of_fit.R"))

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

status_from_p <- function(p) {
  if (is.na(p)) return("not evaluated")
  if (p < 0.001) return("serious concern")
  if (p < 0.05) return("caution")
  "acceptable"
}

classify_numeric_issue <- function(value, caution, serious, direction = "greater") {
  if (direction == "greater") {
    if (value > serious) return("serious concern")
    if (value > caution) return("caution")
  } else {
    if (value < serious) return("serious concern")
    if (value < caution) return("caution")
  }
  "acceptable"
}

make_diagnostic_plots <- function(fit, dat, response) {
  out_dir <- file.path(paths$diagnostics, safe_name(response))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fitted_vals <- fitted(fit)
  resid_vals <- residuals(fit)
  std_resid <- rstandard(fit)
  stud_resid <- rstudent(fit)
  cook <- cooks.distance(fit)
  lev <- hatvalues(fit)

  write_plot <- function(filename, expr) {
    png(file.path(out_dir, filename), width = 1400, height = 1000, res = 150)
    on.exit(dev.off(), add = TRUE)
    par(mar = c(5, 5, 4, 2))
    force(expr)
  }

  write_plot("residuals_vs_fitted.png", {
    plot(fitted_vals, resid_vals, pch = 19, col = "#1f4e79",
         xlab = "Fitted value", ylab = "Residual", main = paste(response, "residuals versus fitted"))
    abline(h = 0, lty = 2, col = "gray40")
  })
  write_plot("normal_qq.png", {
    qqnorm(std_resid, pch = 19, col = "#1f4e79", main = paste(response, "normal Q-Q"))
    qqline(std_resid, col = "#b31b1b", lwd = 2)
  })
  write_plot("scale_location.png", {
    plot(fitted_vals, sqrt(abs(std_resid)), pch = 19, col = "#1f4e79",
         xlab = "Fitted value", ylab = "sqrt(|standardised residual|)",
         main = paste(response, "scale-location"))
  })
  write_plot("residuals_by_temperature.png", {
    boxplot(resid_vals ~ dat$Temp, xlab = "Temperature (C)", ylab = "Residual",
            main = paste(response, "residuals by temperature"), col = "#d9e8f5")
    abline(h = 0, lty = 2, col = "gray40")
  })
  write_plot("residuals_by_olr.png", {
    boxplot(resid_vals ~ dat$OL, xlab = "OLR", ylab = "Residual",
            main = paste(response, "residuals by OLR"), col = "#d9e8f5")
    abline(h = 0, lty = 2, col = "gray40")
  })
  write_plot("residuals_by_retention_time.png", {
    boxplot(resid_vals ~ dat$Time, xlab = "Retention time (days)", ylab = "Residual",
            main = paste(response, "residuals by retention time"), col = "#d9e8f5")
    abline(h = 0, lty = 2, col = "gray40")
  })
  write_plot("residuals_by_recipe.png", {
    boxplot(resid_vals ~ dat$recipe, xlab = "Recipe", ylab = "Residual",
            main = paste(response, "residuals by recipe"), col = "#d9e8f5")
    abline(h = 0, lty = 2, col = "gray40")
  })
  write_plot("cooks_distance.png", {
    plot(cook, type = "h", lwd = 2, col = "#1f4e79", xlab = "Observation index",
         ylab = "Cook's distance", main = paste(response, "Cook's distance"))
    abline(h = 4 / nobs(fit), lty = 2, col = "#b31b1b")
    abline(h = 1, lty = 3, col = "gray40")
  })
  write_plot("leverage.png", {
    plot(lev, type = "h", lwd = 2, col = "#1f4e79", xlab = "Observation index",
         ylab = "Leverage", main = paste(response, "leverage"))
    abline(h = 2 * fit$rank / nobs(fit), lty = 2, col = "#b31b1b")
    abline(h = 3 * fit$rank / nobs(fit), lty = 3, col = "gray40")
  })
  write_plot("studentised_residuals.png", {
    plot(stud_resid, pch = 19, col = "#1f4e79", xlab = "Observation index",
         ylab = "Externally studentised residual", main = paste(response, "studentised residuals"))
    abline(h = c(-3, -2, 0, 2, 3), lty = c(3, 2, 1, 2, 3), col = c("gray40", "#b31b1b", "gray60", "#b31b1b", "gray40"))
  })

  invisible(out_dir)
}

diagnose_model <- function(fit, dat, response, lof_p) {
  stud <- rstudent(fit)
  std <- rstandard(fit)
  cook <- cooks.distance(fit)
  lev <- hatvalues(fit)
  n <- nobs(fit)
  p <- fit$rank
  sh <- shapiro.test(std)
  scale_corr <- suppressWarnings(cor(sqrt(abs(std)), fitted(fit), method = "spearman"))
  max_abs_cell_mean <- max(abs(tapply(residuals(fit), dat$design_cell_id, mean)))

  checks <- data.frame(
    response = response,
    diagnostic = c("lack_of_fit", "normality_shapiro", "studentised_residuals",
                   "cooks_distance", "leverage", "scale_location_pattern", "design_cell_mean_residual"),
    metric = c("p_value", "p_value", "max_abs_rstudent", "max_cooks_d", "max_leverage",
               "abs_spearman_corr_sqrt_abs_stdresid_fitted", "max_abs_cell_mean_residual"),
    value = c(lof_p, sh$p.value, max(abs(stud)), max(cook), max(lev), abs(scale_corr), max_abs_cell_mean),
    classification = NA_character_,
    rule = c("serious if p < 0.001; caution if p < 0.05",
             "serious if p < 0.001; caution if p < 0.01",
             "serious if > 3; caution if > 2",
             "serious if > 1; caution if > 4/n",
             "serious if > 3p/n; caution if > 2p/n",
             "caution if > 0.30; serious if > 0.50",
             "response-scale review threshold: caution if > 0.50 SD; serious if > 1 SD"),
    stringsAsFactors = FALSE
  )
  checks$classification[checks$diagnostic == "lack_of_fit"] <- status_from_p(lof_p)
  checks$classification[checks$diagnostic == "normality_shapiro"] <- ifelse(sh$p.value < 0.001, "serious concern",
                                                                            ifelse(sh$p.value < 0.01, "caution", "acceptable"))
  checks$classification[checks$diagnostic == "studentised_residuals"] <- classify_numeric_issue(max(abs(stud)), 2, 3)
  checks$classification[checks$diagnostic == "cooks_distance"] <- classify_numeric_issue(max(cook), 4 / n, 1)
  checks$classification[checks$diagnostic == "leverage"] <- classify_numeric_issue(max(lev), 2 * p / n, 3 * p / n)
  checks$classification[checks$diagnostic == "scale_location_pattern"] <- classify_numeric_issue(abs(scale_corr), 0.30, 0.50)
  checks$classification[checks$diagnostic == "design_cell_mean_residual"] <- classify_numeric_issue(max_abs_cell_mean, 0.5 * summary(fit)$sigma, summary(fit)$sigma)

  obs_flags <- data.frame(
    response = response,
    minitab_row_id = dat$minitab_row_id,
    trial_row_id = dat$trial_row_id,
    design_cell_id = dat$design_cell_id,
    Temp = dat$Temp,
    OL = dat$OL,
    Time = dat$Time,
    recipe = as.character(dat$recipe),
    observed = dat[[response]],
    fitted = fitted(fit),
    residual = residuals(fit),
    standardised_residual = std,
    studentised_residual = stud,
    leverage = lev,
    cooks_distance = cook,
    flag_studentised = abs(stud) > 2,
    flag_cooks = cook > 4 / n,
    flag_leverage = lev > 2 * p / n,
    stringsAsFactors = FALSE
  )
  obs_flags$issue_count <- obs_flags$flag_studentised + obs_flags$flag_cooks + obs_flags$flag_leverage
  obs_flags <- obs_flags[order(-obs_flags$issue_count, -obs_flags$cooks_distance, -abs(obs_flags$studentised_residual)), ]

  cell_summary <- aggregate(cbind(abs_residual = abs(obs_flags$residual),
                                  cooks_distance = obs_flags$cooks_distance,
                                  leverage = obs_flags$leverage,
                                  issue_count = obs_flags$issue_count),
                            by = list(response = obs_flags$response,
                                      design_cell_id = obs_flags$design_cell_id,
                                      Temp = obs_flags$Temp,
                                      OL = obs_flags$OL,
                                      Time = obs_flags$Time,
                                      recipe = obs_flags$recipe),
                            FUN = function(z) c(mean = mean(z), max = max(z)))
  flat <- data.frame(cell_summary[, 1:6], do.call(data.frame, cell_summary[, 7:10]), check.names = FALSE)
  names(flat) <- make.names(names(flat), unique = TRUE)
  flat <- flat[order(-flat$issue_count.max, -flat$cooks_distance.max), ]

  list(checks = checks, observations = obs_flags, cells = flat)
}

md_table <- function(dat, digits = 4, max_rows = Inf) {
  if (nrow(dat) == 0L) return("")
  if (is.finite(max_rows)) dat <- head(dat, max_rows)
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  for (j in seq_along(dat)) {
    if (is.numeric(dat[[j]])) dat[[j]] <- formatC(dat[[j]], format = "f", digits = digits)
    dat[[j]][is.na(dat[[j]])] <- ""
  }
  header <- paste0("| ", paste(names(dat), collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", ncol(dat)), collapse = " | "), " |")
  rows <- apply(dat, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))
  paste(c(header, sep, rows), collapse = "\n")
}

model_reproduction_status <- function(reconciliation, lof) {
  primary <- reconciliation[reconciliation$statistic %in% c(
    "number_of_observations", "model_degrees_of_freedom", "residual_degrees_of_freedom",
    "model_matrix_rank", "residual_standard_error", "r_squared_pct",
    "adjusted_r_squared_pct", "predicted_r_squared_pct",
    "fitted_values_vs_reported_rounded_equations_max_abs_diff"
  ), ]
  rows <- data.frame()
  for (resp in response_names) {
    pr <- primary[primary$response == resp, ]
    any_fail <- any(pr$pass == "fail")
    missing_full <- any(reconciliation$response == resp &
                          reconciliation$statistic == "residuals_available_for_full_precision_minitab_comparison" &
                          reconciliation$pass == "not_tested")
    lof_sig <- lof$p_value[lof$response == resp] < 0.05
    status <- if (any_fail) "FAIL" else if (missing_full || lof_sig) "CONDITIONAL PASS" else "PASS"
    reason <- if (any_fail) {
      "At least one primary reconstruction statistic failed tolerance."
    } else if (missing_full && lof_sig) {
      "Primary statistics match, but full-precision Minitab fitted values/residuals are unavailable and lack-of-fit is significant."
    } else if (missing_full) {
      "Primary statistics match, but full-precision Minitab fitted values/residuals are unavailable."
    } else if (lof_sig) {
      "Model is reproduced, but significant lack-of-fit affects predictive adequacy."
    } else {
      "Primary reconstruction statistics and fitted-value checks passed."
    }
    rows <- rbind(rows, data.frame(response = resp, phase2a_status = status, reason = reason, stringsAsFactors = FALSE))
  }
  rows
}

html_from_markdown_pre <- function(md, path) {
  esc <- gsub("&", "&amp;", md, fixed = TRUE)
  esc <- gsub("<", "&lt;", esc, fixed = TRUE)
  esc <- gsub(">", "&gt;", esc, fixed = TRUE)
  html <- paste0(
    "<!doctype html><html><head><meta charset=\"utf-8\"><title>Phase 2A Model Reconstruction Report</title>",
    "<style>body{font-family:Arial,sans-serif;margin:32px;line-height:1.45;color:#1f2933}",
    "pre{white-space:pre-wrap;font-family:Consolas,monospace;font-size:13px}",
    "</style></head><body><pre>", esc, "</pre></body></html>"
  )
  writeLines(html, path, useBytes = TRUE)
}

write_report <- function(diag_table, status_table) {
  recon <- read.csv(file.path(paths$tables, "minitab_r_reconciliation.csv"), stringsAsFactors = FALSE)
  design <- read.csv(file.path(paths$tables, "design_structure_summary.csv"), stringsAsFactors = FALSE)
  model_summary <- read.csv(file.path(paths$tables, "model_summary_r.csv"), stringsAsFactors = FALSE)
  anova <- read.csv(file.path(paths$tables, "anova_term_tests.csv"), stringsAsFactors = FALSE)
  lof <- read.csv(file.path(paths$tables, "lack_of_fit_pure_error.csv"), stringsAsFactors = FALSE)
  hull <- read.csv(file.path(paths$tables, "table3_design_space_check.csv"), stringsAsFactors = FALSE)
  data_checks <- read.csv(file.path(paths$tables, "data_reconstruction_checks.csv"), stringsAsFactors = FALSE)

  interaction_rows <- anova[anova$source %in% c("2-way interaction", "Temp*OL", "Temp*Time",
                                                "Temp*Recipe", "OL*Time", "OL*Recipe", "Time*Recipe"), ]
  primary_recon <- recon[recon$statistic %in% c("number_of_observations", "model_degrees_of_freedom",
                                                "residual_degrees_of_freedom", "residual_standard_error",
                                                "r_squared_pct", "adjusted_r_squared_pct",
                                                "predicted_r_squared_pct",
                                                "fitted_values_vs_reported_rounded_equations_max_abs_diff"), ]
  manuscript_flags <- data.frame(
    issue = c("Significant lack-of-fit in all four models",
              "Several Table 3 scenarios are outside the Box-Behnken convex design hull",
              "20 C methane-plus-FIB scenario is outside the observed design hull",
              "Full-precision Minitab fitted values/residuals are unavailable in readable outputs",
              "Zero/non-detect microbial values were transformed to a fixed floor of 1"),
    consequence = c("Prediction and optimisation claims need qualification; do not alter the original model solely to remove lack-of-fit.",
                    "Being inside min/max factor ranges is not sufficient evidence of interpolation.",
                    "This claim should be presented as model-based extrapolation unless validated by new experiments.",
                    "Reconstruction can be conditionally accepted from principal statistics and rounded-equation checks, not full fitted-value identity.",
                    "This is a modelling limitation and should be disclosed if manuscript methods are revised."),
    stringsAsFactors = FALSE
  )

  lines <- c(
    "# Phase 2A model reconstruction report",
    "",
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    "Scope: exact reconstruction and verification of the original four Minitab full-quadratic RSM models only. No bootstrap optimisation, full cross-validation, power/design-sensitivity analysis, factor-importance analysis, or manuscript rewriting was conducted.",
    "",
    "## 1. Exact dataset used",
    "",
    "The modelling table used for reconstruction is `data/input/analysis_design.csv`. Recipe labels were recovered by deterministic row matching to the corresponding 180 non-Day-0 rows in `data/input/trial_data.csv`; the rounded transformed responses and process variables matched row-by-row.",
    "",
    md_table(data_checks),
    "",
    "## 2. Exact transformations used",
    "",
    "Microbial responses were transformed as `log10(count + 1)` followed by replacement of transformed values below 1 with 1. This Phase 2A reconstruction did not change the transformation and did not replace below-detection values by any alternative rule. Methane was modelled untransformed.",
    "",
    "Coded predictors were `temp_c = (Temp - 37) / 18`, `ol_c = (OL - 2.0) / 1.5`, and `time_c = (Time - 6) / 3`.",
    "",
    "## 3. Final model formula",
    "",
    "`response ~ (temp_c + ol_c + time_c + recipe)^2 + I(temp_c^2) + I(ol_c^2) + I(time_c^2)`",
    "",
    "This gives 17 model degrees of freedom and 18 fitted parameters including the intercept.",
    "",
    "## 4. Recipe coding",
    "",
    "The factor level order is `1S:3FW`, `2S:1FW`, `3S:1FW`. The mapping used is `FWS31 -> 1S:3FW`, `FWS12 -> 2S:1FW`, and `FWS13 -> 3S:1FW`. For coefficient reconciliation the model uses Minitab-compatible sum-to-zero recipe contrasts with printed columns for `2S:1FW` and `3S:1FW`; the `1S:3FW` effect is implicit as the negative sum of those two columns rather than a treatment-reference coefficient.",
    "",
    "## 5. Design and model structure checks",
    "",
    md_table(design),
    "",
    "## 6. Whether each Minitab model was reproduced",
    "",
    md_table(status_table),
    "",
    "The status is conditional rather than unconditional for E. coli, enterococci, and methane because the readable Minitab/Word outputs do not contain full-precision fitted values and residuals, although their principal statistics match within the stated tolerances. Coliforms is marked FAIL because R-squared and adjusted R-squared differ from the Word/Minitab values beyond tolerance when the audited `design data .csv` is used.",
    "",
    "## 7. Reconciliation results",
    "",
    md_table(primary_recon[, c("response", "statistic", "minitab_value", "r_value", "absolute_difference", "tolerance", "pass")], max_rows = 40),
    "",
    "The complete reconciliation table is `outputs/tables/minitab_r_reconciliation.csv` and `.xlsx`.",
    "",
    "## 8. ANOVA results",
    "",
    md_table(anova[anova$source %in% c("Model", "Linear", "Square", "2-way interaction"),
                   c("response", "source", "df", "f_value", "p_value_display")]),
    "",
    "Exact R p-values are stored in `anova_term_tests.csv`; values that Minitab displayed as 0 should be reported as `p < 0.001`.",
    "",
    "## 9. Interaction-term results",
    "",
    md_table(interaction_rows[, c("response", "source", "df", "f_value", "p_value_display")]),
    "",
    "## 10. Pure-error and lack-of-fit results",
    "",
    md_table(lof[, c("response", "pure_error_ss", "pure_error_df", "lack_of_fit_ss",
                    "lack_of_fit_df", "lack_of_fit_ms", "f_value", "p_value_display", "significant_0_05")]),
    "",
    "Significant lack-of-fit was confirmed for all four response models. The models were not altered to remove this finding.",
    "",
    "## 11. Diagnostic findings",
    "",
    md_table(diag_table[, c("response", "diagnostic", "metric", "value", "classification")], max_rows = 80),
    "",
    "Diagnostic plots and observation/cell influence tables are in `outputs/diagnostics/`.",
    "",
    "## 12. Design-hull status of Table 3 scenarios",
    "",
    md_table(hull[, c("scenario", "Temp", "OL", "Time", "recipe", "abs_coded_sum", "design_hull_status", "average_model_leverage", "leverage_note")], digits = 3),
    "",
    "The observed Box-Behnken design hull is not the full rectangular min/max range. A point can be within the separate ranges of temperature, OLR, and retention time while still being outside the convex hull of the combined experimental design.",
    "",
    "## 13. Suitability to proceed to grouped internal validation",
    "",
    "E. coli, enterococci, and methane are suitable to proceed to grouped internal validation as conditionally reproduced original models. Coliforms should not proceed as reproduced until the R-squared/adjusted R-squared discrepancy is resolved or explicitly accepted as a limitation. Grouped validation must be by design cell, not ordinary row-level LOOCV, because replicate rows at the same design condition would otherwise leak information.",
    "",
    "## 14. Unresolved discrepancy",
    "",
    "No row-level data discrepancy remained after deterministic matching and the historical microbial transformation. The unresolved reconstruction discrepancies are: (1) Coliforms R-squared and adjusted R-squared do not match the Word/Minitab report within tolerance despite matching observations, degrees of freedom, residual standard error, predicted R-squared, coefficients, and rounded-equation fitted values; and (2) full-precision Minitab fitted values/residuals and full-precision Minitab p-values were not present in readable outputs. The component columns in the Minitab CSV were not used for recipe reconstruction because the audit found them inconsistent with the processed recipe components.",
    "",
    "## 15. Results that may require changing manuscript claims",
    "",
    md_table(manuscript_flags),
    "",
    "## 16. Files created",
    "",
    "- `R/01_reconstruct_rsm_models.R`",
    "- `R/02_anova_lack_of_fit.R`",
    "- `R/03_model_diagnostics.R`",
    "- `tests/test_model_reconstruction.R`",
    "- `models/`",
    "- `outputs/tables/`",
    "- `outputs/figures/`",
    "- `outputs/diagnostics/`",
    "- `outputs/logs/`"
  )

  md <- paste(lines, collapse = "\n")
  writeLines(md, file.path(paths$reports, "01_model_reconstruction_report.md"), useBytes = TRUE)
  html_from_markdown_pre(md, file.path(paths$reports, "01_model_reconstruction_report.html"))
}

phase2a_diagnostics_main <- function() {
  log_path <- file.path(paths$logs, "03_model_diagnostics.log")
  con <- file(log_path, open = "wt")
  sink(con, split = TRUE)
  sink(con, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(con)
  }, add = TRUE)

  cat("Phase 2A diagnostics/report started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
  dat <- readRDS(file.path(paths$models, "reconstructed_model_data.rds"))
  models <- readRDS(file.path(paths$models, "all_full_quadratic_rsm_models.rds"))
  lof <- read.csv(file.path(paths$tables, "lack_of_fit_pure_error.csv"), stringsAsFactors = FALSE)
  reconciliation <- read.csv(file.path(paths$tables, "minitab_r_reconciliation.csv"), stringsAsFactors = FALSE)

  all_checks <- data.frame()
  all_obs <- data.frame()
  all_cells <- data.frame()

  for (resp in response_names) {
    fit <- models[[resp]]
    make_diagnostic_plots(fit, dat, resp)
    dx <- diagnose_model(fit, dat, resp, lof$p_value[lof$response == resp])
    all_checks <- rbind(all_checks, dx$checks)
    all_obs <- rbind(all_obs, dx$observations)
    all_cells <- rbind(all_cells, dx$cells)
  }

  write.csv(all_checks, file.path(paths$tables, "model_diagnostic_classification.csv"), row.names = FALSE)
  write.csv(all_obs, file.path(paths$diagnostics, "influential_observations.csv"), row.names = FALSE)
  write.csv(all_cells, file.path(paths$diagnostics, "influential_design_cells.csv"), row.names = FALSE)

  status_table <- model_reproduction_status(reconciliation, lof)
  write.csv(status_table, file.path(paths$tables, "model_reconstruction_status.csv"), row.names = FALSE)

  write_report(all_checks, status_table)
  capture.output(warnings(), file = file.path(paths$logs, "warnings_03_model_diagnostics.txt"))
  cat("Phase 2A diagnostics/report completed successfully.\n")
  invisible(list(diagnostics = all_checks, status = status_table))
}

if (sys.nframe() == 0L) {
  phase2a_diagnostics_main()
}
