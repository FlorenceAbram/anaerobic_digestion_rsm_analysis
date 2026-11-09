set.seed(20260717)

source(file.path("R", "01_reconstruct_rsm_models.R"))

minitab_anova_p <- list(
  Coliforms = c(Model = 0, Linear = 0, temp_c = 0, ol_c = 0, time_c = 0, recipe = 0,
                Square = 0, `I(temp_c^2)` = 0, `I(ol_c^2)` = 0.904, `I(time_c^2)` = 0.352,
                `2-way interaction` = 0.121, `temp_c:ol_c` = 0.07, `temp_c:time_c` = 0.558,
                `temp_c:recipe` = 0.038, `ol_c:time_c` = 0.065, `ol_c:recipe` = 0.872,
                `time_c:recipe` = 0.874, `Lack-of-fit` = 0),
  Ecoli = c(Model = 0, Linear = 0, temp_c = 0, ol_c = 0, time_c = 0, recipe = 0,
            Square = 0, `I(temp_c^2)` = 0, `I(ol_c^2)` = 0.861, `I(time_c^2)` = 0.333,
            `2-way interaction` = 0.037, `temp_c:ol_c` = 0.002, `temp_c:time_c` = 0.918,
            `temp_c:recipe` = 0.087, `ol_c:time_c` = 0.131, `ol_c:recipe` = 0.752,
            `time_c:recipe` = 0.63, `Lack-of-fit` = 0),
  Enterococci = c(Model = 0, Linear = 0, temp_c = 0, ol_c = 0, time_c = 0.695, recipe = 0,
                  Square = 0, `I(temp_c^2)` = 0, `I(ol_c^2)` = 0, `I(time_c^2)` = 0.994,
                  `2-way interaction` = 0.011, `temp_c:ol_c` = 0.716, `temp_c:time_c` = 0.004,
                  `temp_c:recipe` = 0.019, `ol_c:time_c` = 0.41, `ol_c:recipe` = 0.203,
                  `time_c:recipe` = 0.455, `Lack-of-fit` = 0),
  Methane = c(Model = 0, Linear = 0, temp_c = 0.1, ol_c = 0, time_c = 0, recipe = 0.099,
              Square = 0, `I(temp_c^2)` = 0, `I(ol_c^2)` = 0.026, `I(time_c^2)` = 0.038,
              `2-way interaction` = 0, `temp_c:ol_c` = 0.894, `temp_c:time_c` = 0.639,
              `temp_c:recipe` = 0.771, `ol_c:time_c` = 0, `ol_c:recipe` = 0.07,
              `time_c:recipe` = 0.212, `Lack-of-fit` = 0)
)

fmt_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return("p < 0.001")
  paste0("p = ", formatC(p, format = "f", digits = 3))
}

term_columns <- function(x_names) {
  list(
    Model = setdiff(x_names, "(Intercept)"),
    Linear = c("temp_c", "ol_c", "time_c", "recipe2S:1FW", "recipe3S:1FW"),
    temp_c = "temp_c",
    ol_c = "ol_c",
    time_c = "time_c",
    recipe = c("recipe2S:1FW", "recipe3S:1FW"),
    Square = c("I(temp_c^2)", "I(ol_c^2)", "I(time_c^2)"),
    `I(temp_c^2)` = "I(temp_c^2)",
    `I(ol_c^2)` = "I(ol_c^2)",
    `I(time_c^2)` = "I(time_c^2)",
    `2-way interaction` = c("temp_c:ol_c", "temp_c:time_c", "ol_c:time_c",
                            "temp_c:recipe2S:1FW", "temp_c:recipe3S:1FW",
                            "ol_c:recipe2S:1FW", "ol_c:recipe3S:1FW",
                            "time_c:recipe2S:1FW", "time_c:recipe3S:1FW"),
    `temp_c:ol_c` = "temp_c:ol_c",
    `temp_c:time_c` = "temp_c:time_c",
    `temp_c:recipe` = c("temp_c:recipe2S:1FW", "temp_c:recipe3S:1FW"),
    `ol_c:time_c` = "ol_c:time_c",
    `ol_c:recipe` = c("ol_c:recipe2S:1FW", "ol_c:recipe3S:1FW"),
    `time_c:recipe` = c("time_c:recipe2S:1FW", "time_c:recipe3S:1FW")
  )
}

partial_f <- function(fit, cols) {
  x <- model.matrix(fit)
  y <- model.response(model.frame(fit))
  col_idx <- match(cols, colnames(x))
  if (any(is.na(col_idx))) stop("Requested model column not found: ", paste(cols[is.na(col_idx)], collapse = ", "))
  full_sse <- sum(residuals(fit)^2)
  reduced_x <- x[, -col_idx, drop = FALSE]
  reduced_fit <- lm.fit(reduced_x, y)
  reduced_sse <- sum(reduced_fit$residuals^2)
  df_num <- fit$rank - reduced_fit$rank
  df_den <- df.residual(fit)
  f_value <- ((reduced_sse - full_sse) / df_num) / (full_sse / df_den)
  p_value <- stats::pf(f_value, df_num, df_den, lower.tail = FALSE)
  c(df = df_num, ss = reduced_sse - full_sse, ms = (reduced_sse - full_sse) / df_num,
    f_value = f_value, p_value = p_value)
}

lack_of_fit <- function(fit, dat, response) {
  y <- dat[[response]]
  cell <- dat$design_cell_id
  cell_mean <- ave(y, cell, FUN = mean)
  pure_error_ss <- sum((y - cell_mean)^2)
  pure_error_df <- length(y) - length(unique(cell))
  model_error_ss <- sum(residuals(fit)^2)
  lack_ss <- model_error_ss - pure_error_ss
  lack_df <- length(unique(cell)) - fit$rank
  lack_ms <- lack_ss / lack_df
  pure_error_ms <- pure_error_ss / pure_error_df
  f_value <- lack_ms / pure_error_ms
  p_value <- stats::pf(f_value, lack_df, pure_error_df, lower.tail = FALSE)
  data.frame(
    response = response,
    residual_ss = model_error_ss,
    residual_df = df.residual(fit),
    pure_error_ss = pure_error_ss,
    pure_error_df = pure_error_df,
    pure_error_ms = pure_error_ms,
    lack_of_fit_ss = lack_ss,
    lack_of_fit_df = lack_df,
    lack_of_fit_ms = lack_ms,
    f_value = f_value,
    p_value = p_value,
    p_value_display = fmt_p(p_value),
    significant_0_05 = p_value < 0.05,
    stringsAsFactors = FALSE
  )
}

build_table3 <- function() {
  data.frame(
    scenario = c("EU FIB only, unconstrained temperature",
                 "EU FIB + methane, unconstrained temperature",
                 "EU FIB + methane + fixed 20 C",
                 "EU FIB only + fixed 20 C",
                 "EU FIB only + fixed 25 C",
                 "EU FIB only + fixed 30 C",
                 "EU FIB only + fixed 35 C",
                 "EU FIB only + fixed 45 C",
                 "EU FIB only + fixed 48 C",
                 "EU FIB only + fixed 50 C"),
    Temp = c(49.54, 45.89, 20, 20, 25, 30, 35, 45, 48, 50),
    OL = c(0.5, 3.5, 3.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
    Time = c(5.45, 9.00, 9.00, 4.18, 5.55, 7.00, 9.00, 8.96, 7.50, 4.27),
    recipe = factor(c("1S:3FW", "2S:1FW", "2S:1FW", "2S:1FW", "2S:1FW",
                      "2S:1FW", "2S:1FW", "2S:1FW", "1S:3FW", "1S:3FW"),
                    levels = recipe_levels),
    source = "CSCEE-D-26-00383 Table 3 / audit extraction",
    stringsAsFactors = FALSE
  )
}

classify_hull <- function(newdat, observed) {
  newdat$temp_c <- (newdat$Temp - 37) / 18
  newdat$ol_c <- (newdat$OL - 2.0) / 1.5
  newdat$time_c <- (newdat$Time - 6) / 3
  coded <- as.matrix(newdat[, c("temp_c", "ol_c", "time_c")])
  abs_sum <- rowSums(abs(coded))
  max_abs <- apply(abs(coded), 1, max)
  observed_key <- paste(observed$Temp, observed$OL, observed$Time, observed$recipe, sep = "__")
  new_key <- paste(newdat$Temp, newdat$OL, newdat$Time, newdat$recipe, sep = "__")
  exact_observed <- new_key %in% unique(observed_key)
  inside <- max_abs <= 1 + 1e-10 & abs_sum <= 2 + 1e-10
  boundary <- inside & (abs(max_abs - 1) <= 1e-10 | abs(abs_sum - 2) <= 1e-10)
  status <- ifelse(exact_observed, "observed design point",
                   ifelse(!inside, "extrapolation outside the convex design hull",
                          ifelse(boundary, "boundary prediction", "interpolation inside the convex design hull")))
  cbind(newdat, data.frame(
    abs_coded_sum = abs_sum,
    max_abs_coded = max_abs,
    design_hull_status = status,
    stringsAsFactors = FALSE
  ))
}

prediction_leverage <- function(fit, newdat) {
  tt <- delete.response(terms(fit))
  x0 <- model.matrix(tt, newdat)
  xtxi <- chol2inv(qr.R(qr(model.matrix(fit))))
  as.numeric(rowSums((x0 %*% xtxi) * x0))
}

phase2a_anova_main <- function() {
  log_path <- file.path(paths$logs, "02_anova_lack_of_fit.log")
  con <- file(log_path, open = "wt")
  sink(con, split = TRUE)
  sink(con, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(con)
  }, add = TRUE)

  cat("Phase 2A ANOVA/lack-of-fit started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
  dat <- readRDS(file.path(paths$models, "reconstructed_model_data.rds"))
  models <- readRDS(file.path(paths$models, "all_full_quadratic_rsm_models.rds"))

  anova_rows <- data.frame()
  lof_rows <- data.frame()
  reconciliation <- read.csv(file.path(paths$tables, "minitab_r_reconciliation.csv"), stringsAsFactors = FALSE)

  for (resp in response_names) {
    fit <- models[[resp]]
    x_names <- colnames(model.matrix(fit))
    groups <- term_columns(x_names)
    for (term in names(groups)) {
      pf <- partial_f(fit, groups[[term]])
      term_label <- switch(term,
                           temp_c = "Temp",
                           ol_c = "OL",
                           time_c = "Time",
                           recipe = "Recipe",
                           `I(temp_c^2)` = "Temp*Temp",
                           `I(ol_c^2)` = "OL*OL",
                           `I(time_c^2)` = "Time*Time",
                           `2-way interaction` = "2-way interaction",
                           `temp_c:ol_c` = "Temp*OL",
                           `temp_c:time_c` = "Temp*Time",
                           `temp_c:recipe` = "Temp*Recipe",
                           `ol_c:time_c` = "OL*Time",
                           `ol_c:recipe` = "OL*Recipe",
                           `time_c:recipe` = "Time*Recipe",
                           term)
      anova_rows <- rbind(anova_rows, data.frame(
        response = resp,
        source = term_label,
        df = unname(pf["df"]),
        ss = unname(pf["ss"]),
        ms = unname(pf["ms"]),
        f_value = unname(pf["f_value"]),
        p_value = unname(pf["p_value"]),
        p_value_display = fmt_p(unname(pf["p_value"])),
        test_type = "partial F test against full model subspace",
        stringsAsFactors = FALSE
      ))
      mt_p <- minitab_anova_p[[resp]][[term]]
      if (!is.null(mt_p)) {
        recon_p <- if (mt_p == 0) NA_real_ else mt_p
        explanation <- if (mt_p == 0) {
          "Minitab Word output displays 0; exact R p-value recorded and should be reported as p < 0.001."
        } else {
          "Compared with displayed Minitab term p-value; tolerance allows three-decimal reporting."
        }
        tolerance <- if (mt_p == 0) NA_real_ else 0.0015
        reconciliation <- append_recon(reconciliation, resp, paste0("anova_p_", term_label),
                                       recon_p, unname(pf["p_value"]), tolerance, explanation)
      }
    }
    lof <- lack_of_fit(fit, dat, resp)
    lof_rows <- rbind(lof_rows, lof)
    reconciliation <- append_recon(reconciliation, resp, "pure_error_degrees_of_freedom", 141, lof$pure_error_df, 0,
                                   "Calculated from replicated design cells.")
    reconciliation <- append_recon(reconciliation, resp, "lack_of_fit_degrees_of_freedom", 21, lof$lack_of_fit_df, 0,
                                   "Calculated as design-cell means DF minus model rank.")
    reconciliation <- append_recon(reconciliation, resp, "lack_of_fit_p_value", NA_real_, lof$p_value, NA_real_,
                                   "Minitab Word output displays 0; exact R p-value recorded and should be reported as p < 0.001.")
  }

  write.csv(anova_rows, file.path(paths$tables, "anova_term_tests.csv"), row.names = FALSE)
  write.csv(lof_rows, file.path(paths$tables, "lack_of_fit_pure_error.csv"), row.names = FALSE)

  table3 <- build_table3()
  hull <- classify_hull(table3, dat)
  levs <- sapply(models, prediction_leverage, newdat = hull)
  colnames(levs) <- paste0("leverage_", names(models))
  hull <- cbind(hull, as.data.frame(levs))
  hull$average_model_leverage <- rowMeans(levs)
  hull$leverage_note <- ifelse(hull$average_model_leverage > 3 * 18 / 180, "high relative to 3p/n",
                               ifelse(hull$average_model_leverage > 2 * 18 / 180, "caution relative to 2p/n", "within routine range"))
  write.csv(hull, file.path(paths$tables, "table3_design_space_check.csv"), row.names = FALSE)

  write.csv(reconciliation, file.path(paths$tables, "minitab_r_reconciliation.csv"), row.names = FALSE)
  write_simple_xlsx(reconciliation, file.path(paths$tables, "minitab_r_reconciliation.xlsx"))

  capture.output(warnings(), file = file.path(paths$logs, "warnings_02_anova_lack_of_fit.txt"))
  cat("Phase 2A ANOVA/lack-of-fit completed successfully.\n")
  invisible(list(anova = anova_rows, lack_of_fit = lof_rows, hull = hull))
}

if (sys.nframe() == 0L) {
  phase2a_anova_main()
}
