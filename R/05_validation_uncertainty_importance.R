set.seed(20260717)

source(file.path("R", "04_coliform_and_design_space.R"))

phase2b_dir <- file.path(rv_dir, "outputs", "phase2b")
supp_fig_dir <- file.path(paths$figures, "supplementary")
dir.create(phase2b_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_fig_dir, recursive = TRUE, showWarnings = FALSE)

fib_responses <- c("Coliforms", "Ecoli", "Enterococci")
regulatory_log_threshold <- log10(1000 + 1)
strict_table3_log_threshold <- 2.5

prep_model_data <- function(dat) {
  dat$recipe <- factor(dat$recipe, levels = recipe_levels)
  contrasts(dat$recipe) <- recipe_contrasts()
  dat
}

prep_newdata <- function(newdat) {
  newdat$temp_c <- (newdat$Temp - 37) / 18
  newdat$ol_c <- (newdat$OL - 2.0) / 1.5
  newdat$time_c <- (newdat$Time - 6) / 3
  newdat$recipe <- factor(as.character(newdat$recipe), levels = recipe_levels)
  contrasts(newdat$recipe) <- recipe_contrasts()
  newdat
}

fit_response <- function(dat, response) {
  d <- prep_model_data(dat)
  d$response <- d[[response]]
  lm(model_formula, data = d)
}

metric_summary <- function(observed, predicted) {
  err <- predicted - observed
  data.frame(
    n = length(observed),
    rmse = sqrt(mean(err^2)),
    mae = mean(abs(err)),
    bias = mean(err),
    q2 = 1 - sum(err^2) / sum((observed - mean(observed))^2),
    stringsAsFactors = FALSE
  )
}

aggregate_cell_mean <- function(dat, observed, predicted) {
  tmp <- data.frame(
    design_cell_id = dat$design_cell_id,
    Temp = dat$Temp,
    OL = dat$OL,
    Time = dat$Time,
    recipe = as.character(dat$recipe),
    observed = observed,
    predicted = predicted,
    stringsAsFactors = FALSE
  )
  aggregate(cbind(observed, predicted) ~ design_cell_id + Temp + OL + Time + recipe, tmp, mean)
}

grouped_condition_cv <- function(dat, response) {
  cells <- sort(unique(dat$design_cell_id))
  pred <- rep(NA_real_, nrow(dat))
  fold_rank <- data.frame()
  for (cell in cells) {
    train <- dat[dat$design_cell_id != cell, ]
    test <- dat[dat$design_cell_id == cell, ]
    fit <- fit_response(train, response)
    nd <- prep_newdata(test)
    pred[dat$design_cell_id == cell] <- as.numeric(predict(fit, newdata = nd))
    fold_rank <- rbind(fold_rank, data.frame(
      response = response,
      held_out_design_cell_id = cell,
      training_n = nrow(train),
      test_n = nrow(test),
      model_rank = fit$rank,
      residual_df = df.residual(fit),
      stringsAsFactors = FALSE
    ))
  }
  obs_metrics <- metric_summary(dat[[response]], pred)
  obs_metrics$response <- response
  obs_metrics$validation_unit <- "observation, grouped by complete design condition"
  cell_dat <- aggregate_cell_mean(dat, dat[[response]], pred)
  cell_metrics <- metric_summary(cell_dat$observed, cell_dat$predicted)
  cell_metrics$response <- response
  cell_metrics$validation_unit <- "design-condition mean"
  list(predictions = data.frame(
    response = response,
    minitab_row_id = dat$minitab_row_id,
    design_cell_id = dat$design_cell_id,
    Temp = dat$Temp,
    OL = dat$OL,
    Time = dat$Time,
    recipe = as.character(dat$recipe),
    observed = dat[[response]],
    grouped_cv_predicted = pred,
    grouped_cv_residual = dat[[response]] - pred,
    stringsAsFactors = FALSE
  ), cell_predictions = cbind(response = response, cell_dat),
  metrics = rbind(obs_metrics, cell_metrics), fold_rank = fold_rank)
}

row_loo_metrics <- function(models, dat) {
  out <- data.frame()
  for (resp in response_names) {
    fit <- models[[resp]]
    loo_resid <- residuals(fit) / (1 - hatvalues(fit))
    loo_pred <- dat[[resp]] - loo_resid
    m <- metric_summary(dat[[resp]], loo_pred)
    m$response <- resp
    m$validation_unit <- "ordinary row-level leave-one-out; replicate leakage risk"
    out <- rbind(out, m)
  }
  out
}

make_interval_table <- function(models, table3) {
  nd <- prep_newdata(table3)
  rows <- data.frame()
  for (resp in response_names) {
    ci <- as.data.frame(predict(models[[resp]], newdata = nd, interval = "confidence", level = 0.95))
    pi <- as.data.frame(predict(models[[resp]], newdata = nd, interval = "prediction", level = 0.95))
    is_fib <- resp %in% fib_responses
    regulatory_threshold <- if (is_fib) rep(regulatory_log_threshold, nrow(table3)) else rep(NA_real_, nrow(table3))
    strict_threshold <- if (is_fib) rep(strict_table3_log_threshold, nrow(table3)) else rep(NA_real_, nrow(table3))
    rows <- rbind(rows, data.frame(
      scenario = table3$scenario,
      Temp = table3$Temp,
      OL = table3$OL,
      Time = table3$Time,
      recipe = as.character(table3$recipe),
      response = resp,
      fit = ci$fit,
      confidence_lwr = ci$lwr,
      confidence_upr = ci$upr,
      prediction_lwr = pi$lwr,
      prediction_upr = pi$upr,
      regulatory_threshold_log10_1000_cfu_g = regulatory_threshold,
      strict_table3_threshold_log10 = strict_threshold,
      point_meets_regulatory_threshold = if (is_fib) ci$fit <= regulatory_log_threshold else rep(NA, nrow(table3)),
      confidence_upr_meets_regulatory_threshold = if (is_fib) ci$upr <= regulatory_log_threshold else rep(NA, nrow(table3)),
      prediction_upr_meets_regulatory_threshold = if (is_fib) pi$upr <= regulatory_log_threshold else rep(NA, nrow(table3)),
      point_meets_strict_2_5_threshold = if (is_fib) ci$fit <= strict_table3_log_threshold else rep(NA, nrow(table3)),
      confidence_upr_meets_strict_2_5_threshold = if (is_fib) ci$upr <= strict_table3_log_threshold else rep(NA, nrow(table3)),
      prediction_upr_meets_strict_2_5_threshold = if (is_fib) pi$upr <= strict_table3_log_threshold else rep(NA, nrow(table3)),
      stringsAsFactors = FALSE
    ))
  }
  rows
}

make_prediction_grid <- function() {
  g <- expand.grid(
    Temp = seq(19, 55, by = 2),
    OL = seq(0.5, 3.5, by = 0.25),
    Time = seq(3, 9, by = 0.5),
    recipe = recipe_levels,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  g <- prep_newdata(g)
  g$abs_coded_sum <- abs(g$temp_c) + abs(g$ol_c) + abs(g$time_c)
  g$max_abs_coded <- pmax(abs(g$temp_c), abs(g$ol_c), abs(g$time_c))
  g[g$max_abs_coded <= 1 + 1e-10 & g$abs_coded_sum <= 2 + 1e-10, ]
}

predict_all_responses <- function(fits, newdat) {
  nd <- prep_newdata(newdat)
  out <- data.frame(row_id = seq_len(nrow(nd)))
  for (resp in response_names) out[[resp]] <- as.numeric(predict(fits[[resp]], newdata = nd))
  cbind(newdat[, c("Temp", "OL", "Time", "recipe")], out[, response_names])
}

choose_optima <- function(preds, regulatory_threshold = regulatory_log_threshold, strict_threshold = strict_table3_log_threshold) {
  p <- preds
  p$max_fib <- apply(p[, fib_responses], 1, max)
  p$boundary <- abs((p$Temp - 37) / 18) + abs((p$OL - 2) / 1.5) + abs((p$Time - 6) / 3) >= 1.95 |
    pmax(abs((p$Temp - 37) / 18), abs((p$OL - 2) / 1.5), abs((p$Time - 6) / 3)) >= 0.95
  ord_fib <- order(p$max_fib, -p$Methane)
  fib_min <- p[ord_fib[1], ]
  reg_ok <- p[p$max_fib <= regulatory_threshold, ]
  reg_best <- if (nrow(reg_ok) > 0) reg_ok[order(-reg_ok$Methane, reg_ok$max_fib)[1], ] else fib_min
  strict_ok <- p[p$max_fib <= strict_threshold, ]
  strict_best <- if (nrow(strict_ok) > 0) strict_ok[order(-strict_ok$Methane, strict_ok$max_fib)[1], ] else fib_min
  rbind(
    cbind(criterion = "minimise worst FIB", fib_min),
    cbind(criterion = "maximise methane subject to regulatory FIB", reg_best),
    cbind(criterion = "maximise methane subject to stricter 2.5 log FIB", strict_best)
  )
}

bootstrap_models <- function(dat, B = 300) {
  split_idx <- split(seq_len(nrow(dat)), dat$design_cell_id)
  out <- vector("list", B)
  for (b in seq_len(B)) {
    idx <- unlist(lapply(split_idx, function(ii) sample(ii, length(ii), replace = TRUE)), use.names = FALSE)
    boot_dat <- dat[idx, ]
    out[[b]] <- lapply(response_names, function(resp) fit_response(boot_dat, resp))
    names(out[[b]]) <- response_names
  }
  out
}

summarise_bootstrap_table3 <- function(boot_fits, table3) {
  nd <- prep_newdata(table3)
  rows <- data.frame()
  for (b in seq_along(boot_fits)) {
    for (resp in response_names) {
      rows <- rbind(rows, data.frame(
        bootstrap_id = b,
        scenario = table3$scenario,
        response = resp,
        predicted = as.numeric(predict(boot_fits[[b]][[resp]], newdata = nd)),
        stringsAsFactors = FALSE
      ))
    }
  }
  rows
}

summarise_quantiles <- function(x) {
  c(mean = mean(x), sd = sd(x), q025 = unname(quantile(x, 0.025)),
    median = unname(quantile(x, 0.5)), q975 = unname(quantile(x, 0.975)))
}

bootstrap_optima <- function(boot_fits, grid) {
  rows <- data.frame()
  for (b in seq_along(boot_fits)) {
    preds <- predict_all_responses(boot_fits[[b]], grid)
    opt <- choose_optima(preds)
    opt$bootstrap_id <- b
    rows <- rbind(rows, opt)
  }
  rows
}

summarise_optima <- function(opt) {
  rows <- data.frame()
  for (crit in unique(opt$criterion)) {
    d <- opt[opt$criterion == crit, ]
    mode_recipe <- names(sort(table(as.character(d$recipe)), decreasing = TRUE))[1]
    rows <- rbind(rows, data.frame(
      criterion = crit,
      n_bootstrap = nrow(d),
      temp_mean = mean(d$Temp),
      temp_q025 = unname(quantile(d$Temp, 0.025)),
      temp_median = unname(quantile(d$Temp, 0.5)),
      temp_q975 = unname(quantile(d$Temp, 0.975)),
      ol_mean = mean(d$OL),
      ol_q025 = unname(quantile(d$OL, 0.025)),
      ol_median = unname(quantile(d$OL, 0.5)),
      ol_q975 = unname(quantile(d$OL, 0.975)),
      time_mean = mean(d$Time),
      time_q025 = unname(quantile(d$Time, 0.025)),
      time_median = unname(quantile(d$Time, 0.5)),
      time_q975 = unname(quantile(d$Time, 0.975)),
      modal_recipe = mode_recipe,
      modal_recipe_pct = 100 * max(table(as.character(d$recipe))) / nrow(d),
      boundary_pct = 100 * mean(as.logical(d$boundary)),
      max_fib_mean = mean(d$max_fib),
      methane_mean = mean(d$Methane),
      stringsAsFactors = FALSE
    ))
  }
  rows
}

factor_importance <- function(dat) {
  formulas <- list(
    Temperature = response ~ (ol_c + time_c + recipe)^2 + I(ol_c^2) + I(time_c^2),
    Organic_loading = response ~ (temp_c + time_c + recipe)^2 + I(temp_c^2) + I(time_c^2),
    Retention_time = response ~ (temp_c + ol_c + recipe)^2 + I(temp_c^2) + I(ol_c^2),
    Recipe = response ~ (temp_c + ol_c + time_c)^2 + I(temp_c^2) + I(ol_c^2) + I(time_c^2)
  )
  rows <- data.frame()
  for (resp in response_names) {
    d <- prep_model_data(dat)
    d$response <- d[[resp]]
    full <- lm(model_formula, data = d)
    full_sse <- sum(residuals(full)^2)
    total_delta <- 0
    tmp <- data.frame()
    for (fac in names(formulas)) {
      red <- lm(formulas[[fac]], data = d)
      a <- anova(red, full)
      delta <- sum(residuals(red)^2) - full_sse
      total_delta <- total_delta + delta
      tmp <- rbind(tmp, data.frame(
        response = resp,
        factor_group = fac,
        reduced_model_df = red$rank - 1,
        removed_df = full$rank - red$rank,
        delta_sse = delta,
        partial_f = a$F[2],
        p_value = a$`Pr(>F)`[2],
        partial_eta_squared = delta / (delta + full_sse),
        stringsAsFactors = FALSE
      ))
    }
    tmp$hierarchical_importance_pct <- 100 * tmp$delta_sse / sum(tmp$delta_sse)
    rows <- rbind(rows, tmp)
  }
  rows
}

replicate_sensitivity <- function(dat, models, table3, iterations = 80) {
  split_idx <- split(seq_len(nrow(dat)), dat$design_cell_id)
  sizes <- data.frame(noncentre_n = c(1, 2, 3, 4), centre_n = c(3, 6, 9, 12))
  nd <- prep_newdata(table3)
  full_pred <- data.frame()
  for (resp in response_names) {
    full_pred <- rbind(full_pred, data.frame(
      response = resp,
      scenario = table3$scenario,
      full_prediction = as.numeric(predict(models[[resp]], newdata = nd)),
      stringsAsFactors = FALSE
    ))
  }
  rows <- data.frame()
  for (s in seq_len(nrow(sizes))) {
    for (i in seq_len(iterations)) {
      idx <- unlist(lapply(split_idx, function(ii) {
        n_take <- if (length(ii) > 4) sizes$centre_n[s] else sizes$noncentre_n[s]
        sample(ii, min(n_take, length(ii)), replace = FALSE)
      }), use.names = FALSE)
      sub_dat <- dat[idx, ]
      for (resp in response_names) {
        fit <- fit_response(sub_dat, resp)
        pr <- as.numeric(predict(fit, newdata = nd))
        fp <- full_pred[full_pred$response == resp, ]
        rows <- rbind(rows, data.frame(
          response = resp,
          iteration = i,
          noncentre_n = sizes$noncentre_n[s],
          centre_n = sizes$centre_n[s],
          table3_rmse_vs_full_model = sqrt(mean((pr - fp$full_prediction)^2)),
          table3_max_abs_drift_vs_full_model = max(abs(pr - fp$full_prediction)),
          model_rank = fit$rank,
          residual_df = df.residual(fit),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  rows
}

write_figures <- function(cv_cell, intervals, importance, optima, sensitivity) {
  png(file.path(supp_fig_dir, "S1_grouped_condition_cv_observed_vs_predicted.png"), width = 1800, height = 1400, res = 180)
  par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))
  for (resp in response_names) {
    d <- cv_cell[cv_cell$response == resp, ]
    lim <- range(c(d$observed, d$predicted), na.rm = TRUE)
    plot(d$observed, d$predicted, pch = 19, col = "#1F4E79", xlab = "Observed design-condition mean",
         ylab = "Grouped-CV predicted mean", main = resp, xlim = lim, ylim = lim)
    abline(0, 1, col = "#9B1C1C", lwd = 2)
  }
  dev.off()

  fib <- intervals[intervals$response %in% fib_responses, ]
  scenario_order <- unique(fib$scenario)
  worst <- data.frame()
  for (sc in scenario_order) {
    d <- fib[fib$scenario == sc, ]
    worst <- rbind(worst, data.frame(
      scenario = sc,
      fit = max(d$fit),
      confidence_upr = max(d$confidence_upr),
      prediction_upr = max(d$prediction_upr),
      stringsAsFactors = FALSE
    ))
  }
  png(file.path(supp_fig_dir, "S2_table3_worst_fib_prediction_intervals.png"), width = 1900, height = 1100, res = 180)
  par(mar = c(8, 5, 3, 1))
  x <- seq_len(nrow(worst))
  ylim <- range(c(worst$fit, worst$confidence_upr, worst$prediction_upr, regulatory_log_threshold, strict_table3_log_threshold), na.rm = TRUE)
  plot(x, worst$fit, pch = 19, col = "#1F4E79", ylim = ylim, xaxt = "n", xlab = "", ylab = "Worst predicted FIB log10(CFU/g)",
       main = "Table 3 worst-FIB uncertainty")
  segments(x, worst$fit, x, worst$confidence_upr, col = "#1F4E79", lwd = 3)
  segments(x + 0.12, worst$fit, x + 0.12, worst$prediction_upr, col = "#9B1C1C", lwd = 2)
  abline(h = regulatory_log_threshold, col = "#7A5A00", lwd = 2, lty = 2)
  abline(h = strict_table3_log_threshold, col = "#555555", lwd = 2, lty = 3)
  axis(1, at = x, labels = paste0("T3-", x), las = 2)
  legend("topright", legend = c("Mean to upper 95% confidence", "Mean to upper 95% prediction", "1000 CFU/g", "2.5 log"),
         col = c("#1F4E79", "#9B1C1C", "#7A5A00", "#555555"), lty = c(1, 1, 2, 3), lwd = c(3, 2, 2, 2), bty = "n")
  dev.off()

  png(file.path(supp_fig_dir, "S3_hierarchical_factor_importance.png"), width = 1800, height = 1200, res = 180)
  par(mfrow = c(2, 2), mar = c(5, 5, 3, 1))
  for (resp in response_names) {
    d <- importance[importance$response == resp, ]
    cols <- c("#1F4E79", "#6A8A3A", "#7A5A00", "#9B1C1C")
    barplot(d$hierarchical_importance_pct, names.arg = gsub("_", "\n", d$factor_group), col = cols,
            ylim = c(0, max(100, d$hierarchical_importance_pct) * 1.05), ylab = "Importance (%)", main = resp)
  }
  dev.off()

  png(file.path(supp_fig_dir, "S4_bootstrap_optimum_uncertainty.png"), width = 1800, height = 1300, res = 180)
  par(mfrow = c(2, 2), mar = c(5, 5, 3, 1))
  crits <- unique(optima$criterion)
  boxplot(Temp ~ criterion, data = optima, las = 2, col = "#D9E8F5", ylab = "Temperature (C)", main = "Bootstrap optima: temperature")
  boxplot(OL ~ criterion, data = optima, las = 2, col = "#E8EEF5", ylab = "OLR", main = "Bootstrap optima: OLR")
  boxplot(Time ~ criterion, data = optima, las = 2, col = "#F2F4F7", ylab = "Retention time (days)", main = "Bootstrap optima: retention time")
  barplot(t(table(optima$recipe, optima$criterion)), beside = TRUE, las = 2, col = c("#1F4E79", "#6A8A3A", "#9B1C1C"),
          ylab = "Bootstrap selections", main = "Bootstrap optima: recipe")
  dev.off()

  sens_summary <- aggregate(cbind(table3_rmse_vs_full_model, table3_max_abs_drift_vs_full_model) ~ response + noncentre_n + centre_n,
                            sensitivity, mean)
  png(file.path(supp_fig_dir, "S5_retrospective_design_sensitivity.png"), width = 1700, height = 1100, res = 180)
  par(mar = c(5, 5, 3, 1))
  plot(NA, xlim = range(sens_summary$noncentre_n), ylim = range(sens_summary$table3_rmse_vs_full_model),
       xlab = "Non-centre replicates retained per design condition", ylab = "Mean Table 3 prediction RMSE vs full model",
       main = "Retrospective replicate-count sensitivity")
  cols <- c(Coliforms = "#1F4E79", Ecoli = "#6A8A3A", Enterococci = "#7A5A00", Methane = "#9B1C1C")
  for (resp in response_names) {
    d <- sens_summary[sens_summary$response == resp, ]
    lines(d$noncentre_n, d$table3_rmse_vs_full_model, type = "b", pch = 19, lwd = 2, col = cols[[resp]])
  }
  legend("topright", legend = response_names, col = cols, pch = 19, lwd = 2, bty = "n")
  dev.off()
}

phase2b_main <- function() {
  log_path <- file.path(paths$logs, "05_validation_uncertainty_importance.log")
  con <- file(log_path, open = "wt")
  sink(con, split = TRUE)
  sink(con, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(con)
  }, add = TRUE)

  cat("Phase 2B validation/uncertainty/importance started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
  dat <- readRDS(file.path(paths$models, "reconstructed_model_data.rds"))
  dat <- prep_model_data(dat)
  models <- readRDS(file.path(paths$models, "all_full_quadratic_rsm_models.rds"))
  table3 <- build_table3()
  table3_verified <- read.csv(file.path(paths$tables, "table3_design_space_check_verified.csv"), stringsAsFactors = FALSE)

  cv_preds <- data.frame()
  cv_cells <- data.frame()
  cv_metrics <- data.frame()
  fold_rank <- data.frame()
  for (resp in response_names) {
    cv <- grouped_condition_cv(dat, resp)
    cv_preds <- rbind(cv_preds, cv$predictions)
    cv_cells <- rbind(cv_cells, cv$cell_predictions)
    cv_metrics <- rbind(cv_metrics, cv$metrics)
    fold_rank <- rbind(fold_rank, cv$fold_rank)
  }
  row_metrics <- row_loo_metrics(models, dat)
  validation_metrics <- rbind(cv_metrics, row_metrics)

  intervals <- make_interval_table(models, table3)
  grid <- make_prediction_grid()
  full_grid_predictions <- predict_all_responses(models, grid)
  full_optima <- choose_optima(full_grid_predictions)
  boot_fits <- bootstrap_models(dat, B = 300)
  boot_table3 <- summarise_bootstrap_table3(boot_fits, table3)
  boot_summary <- aggregate(predicted ~ scenario + response, boot_table3, function(x) paste(summarise_quantiles(x), collapse = "|"))
  boot_quant <- do.call(rbind, strsplit(as.character(boot_summary$predicted), "\\|", fixed = FALSE))
  boot_quant <- as.data.frame(apply(boot_quant, 2, as.numeric))
  names(boot_quant) <- c("bootstrap_mean", "bootstrap_sd", "bootstrap_q025", "bootstrap_median", "bootstrap_q975")
  boot_summary <- cbind(boot_summary[, c("scenario", "response")], boot_quant)

  boot_wide <- reshape(boot_table3, idvar = c("bootstrap_id", "scenario"), timevar = "response", direction = "wide")
  names(boot_wide) <- sub("^predicted\\.", "", names(boot_wide))
  boot_wide$max_fib <- apply(boot_wide[, fib_responses], 1, max)
  threshold_probs <- aggregate(cbind(regulatory_pass = max_fib <= regulatory_log_threshold,
                                     strict_2_5_pass = max_fib <= strict_table3_log_threshold) ~ scenario,
                               boot_wide, mean)
  threshold_probs$regulatory_pass <- 100 * threshold_probs$regulatory_pass
  threshold_probs$strict_2_5_pass <- 100 * threshold_probs$strict_2_5_pass

  boot_opt <- bootstrap_optima(boot_fits, grid)
  opt_summary <- summarise_optima(boot_opt)
  importance <- factor_importance(dat)
  sensitivity <- replicate_sensitivity(dat, models, table3, iterations = 80)
  sensitivity_summary <- aggregate(cbind(table3_rmse_vs_full_model, table3_max_abs_drift_vs_full_model, model_rank, residual_df) ~
                                     response + noncentre_n + centre_n, sensitivity, mean)

  write.csv(cv_preds, file.path(paths$tables, "grouped_condition_cv_predictions.csv"), row.names = FALSE)
  write.csv(cv_cells, file.path(paths$tables, "grouped_condition_cv_cell_predictions.csv"), row.names = FALSE)
  write.csv(validation_metrics, file.path(paths$tables, "validation_metrics_grouped_vs_row_loo.csv"), row.names = FALSE)
  write.csv(fold_rank, file.path(phase2b_dir, "grouped_cv_fold_model_rank.csv"), row.names = FALSE)
  write.csv(intervals, file.path(paths$tables, "table3_prediction_confidence_intervals.csv"), row.names = FALSE)
  write.csv(full_optima, file.path(paths$tables, "grid_optima_full_model.csv"), row.names = FALSE)
  write.csv(boot_table3, file.path(phase2b_dir, "table3_bootstrap_predictions_long.csv"), row.names = FALSE)
  write.csv(boot_summary, file.path(paths$tables, "table3_bootstrap_prediction_uncertainty.csv"), row.names = FALSE)
  write.csv(threshold_probs, file.path(paths$tables, "table3_bootstrap_threshold_probabilities.csv"), row.names = FALSE)
  write.csv(boot_opt, file.path(phase2b_dir, "bootstrap_optimum_conditions.csv"), row.names = FALSE)
  write.csv(opt_summary, file.path(paths$tables, "bootstrap_optimum_uncertainty_summary.csv"), row.names = FALSE)
  write.csv(importance, file.path(paths$tables, "hierarchical_factor_importance.csv"), row.names = FALSE)
  write.csv(sensitivity, file.path(phase2b_dir, "replicate_count_sensitivity_long.csv"), row.names = FALSE)
  write.csv(sensitivity_summary, file.path(paths$tables, "retrospective_design_sensitivity_summary.csv"), row.names = FALSE)

  write_figures(cv_cells, intervals, importance, boot_opt, sensitivity)

  grouped_only <- validation_metrics[validation_metrics$validation_unit != "ordinary row-level leave-one-out; replicate leakage risk", ]
  lof <- read.csv(file.path(paths$tables, "lack_of_fit_pure_error.csv"), stringsAsFactors = FALSE)
  status <- read.csv(file.path(paths$tables, "model_reconstruction_status.csv"), stringsAsFactors = FALSE)

  report <- c(
    "# Phase 2B validation, uncertainty and factor-importance report",
    "",
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    "Scope: Reviewer 1 Comment 4 and Reviewer 2 Comments 3, 6, 7 and 13 only. This is internal validation of the existing full-quadratic RSM models; no external experiment is claimed.",
    "",
    "## Model reconstruction status carried forward",
    "",
    md_table(status, digits = 4),
    "",
    "Coliforms is treated as a conditional-use model after Phase 2A.1: the fitted R model is internally consistent and reproduces the available data/model outputs except for the displayed R-squared and adjusted R-squared values, which are mathematically inconsistent with the available response SST and reported residual standard error.",
    "",
    "## Grouped internal validation by complete design condition",
    "",
    md_table(grouped_only[, c("response", "validation_unit", "n", "rmse", "mae", "bias", "q2")], digits = 4),
    "",
    "Grouped leave-one-design-condition-out validation is the appropriate internal validation because all biological replicates from the same temperature, OLR, retention-time and recipe condition must be held out together. Ordinary row-level leave-one-out is also tabulated for comparison but should not be used as the primary validation statistic because replicate rows from the same design condition remain in training.",
    "",
    "## Lack-of-fit",
    "",
    md_table(lof[, c("response", "pure_error_df", "lack_of_fit_df", "f_value", "p_value_display", "significant_0_05")], digits = 4),
    "",
    "The significant pure-error lack-of-fit tests mean that high R-squared values should not be presented alone as proof that the quadratic surface fully captures the replicated design-cell means.",
    "",
    "## Table 3 prediction and uncertainty",
    "",
    md_table(intervals[intervals$response %in% fib_responses,
                       c("scenario", "response", "fit", "confidence_lwr", "confidence_upr", "prediction_lwr", "prediction_upr",
                         "point_meets_regulatory_threshold", "confidence_upr_meets_regulatory_threshold",
                         "point_meets_strict_2_5_threshold", "confidence_upr_meets_strict_2_5_threshold")],
             digits = 4, max_rows = 40),
    "",
    paste0("The regulatory FIB threshold is assessed as log10(1000 + 1) = ", formatC(regulatory_log_threshold, digits = 4, format = "f"), ". The stricter Table 3 optimisation criterion is assessed separately as 2.5 log10(CFU/g). A point estimate below either threshold is not proof of compliance; the confidence and prediction intervals quantify uncertainty around that point estimate."),
    "",
    "## Table 3 design-hull status",
    "",
    md_table(table3_verified[, c("scenario", "Temp", "OL", "Time", "recipe", "inside_recipe_specific_convex_hull", "on_hull_boundary", "outside_hull", "extrapolation_drivers")], digits = 4),
    "",
    "## Bootstrap uncertainty around optima",
    "",
    md_table(opt_summary, digits = 3),
    "",
    "The bootstrap resampled observations within each complete design condition, preserving the original design settings and replicate counts. The optimisation grid was restricted to the observed Box-Behnken design hull, so these summaries quantify internal optimum stability rather than validated operating optima for plants.",
    "",
    "## Retrospective design-sensitivity analysis",
    "",
    md_table(sensitivity_summary, digits = 4),
    "",
    "This replaces a misleading post-hoc power calculation. It asks how much the Table 3 predictions would drift if fewer biological replicates had been available per design condition, while preserving the same design geometry.",
    "",
    "## Hierarchical operational-factor importance",
    "",
    md_table(importance[, c("response", "factor_group", "removed_df", "delta_sse", "partial_f", "p_value", "partial_eta_squared", "hierarchical_importance_pct")], digits = 4),
    "",
    "Factor importance was calculated by removing each operational factor and all hierarchy-related terms involving that factor, then comparing the reduced model to the full quadratic model.",
    "",
    "## Files created",
    "",
    "- `outputs/tables/grouped_condition_cv_predictions.csv`",
    "- `outputs/tables/grouped_condition_cv_cell_predictions.csv`",
    "- `outputs/tables/validation_metrics_grouped_vs_row_loo.csv`",
    "- `outputs/tables/table3_prediction_confidence_intervals.csv`",
    "- `outputs/tables/table3_bootstrap_prediction_uncertainty.csv`",
    "- `outputs/tables/table3_bootstrap_threshold_probabilities.csv`",
    "- `outputs/tables/bootstrap_optimum_uncertainty_summary.csv`",
    "- `outputs/tables/retrospective_design_sensitivity_summary.csv`",
    "- `outputs/tables/hierarchical_factor_importance.csv`",
    "- `outputs/figures/supplementary/S1_grouped_condition_cv_observed_vs_predicted.png`",
    "- `outputs/figures/supplementary/S2_table3_worst_fib_prediction_intervals.png`",
    "- `outputs/figures/supplementary/S3_hierarchical_factor_importance.png`",
    "- `outputs/figures/supplementary/S4_bootstrap_optimum_uncertainty.png`",
    "- `outputs/figures/supplementary/S5_retrospective_design_sensitivity.png`"
  )
  write_md(report, file.path(paths$reports, "02_validation_uncertainty_importance_report.md"))
  capture.output(warnings(), file = file.path(paths$logs, "warnings_05_validation_uncertainty_importance.txt"))
  cat("Phase 2B validation/uncertainty/importance completed successfully.\n")
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  phase2b_main()
}
