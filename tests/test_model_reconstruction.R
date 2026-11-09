set.seed(20260717)

source(file.path("R", "03_model_diagnostics.R"))

failures <- character()
check <- function(condition, message) {
  if (!isTRUE(condition)) {
    failures <<- c(failures, message)
    message("FAIL: ", message)
  } else {
    message("PASS: ", message)
  }
}

input_files <- unlist(paths[c("minitab_data", "trial_data", "cleaned_data")])
pre_hash <- tools::md5sum(input_files)

phase2a_reconstruct_main()
phase2a_anova_main()
phase2a_diagnostics_main()

post_hash <- tools::md5sum(input_files)

saved_hash <- read.csv(file.path(paths$tables, "input_file_checksums.csv"), stringsAsFactors = FALSE)
reconstructed <- readRDS(file.path(paths$models, "reconstructed_model_data.rds"))
models <- readRDS(file.path(paths$models, "all_full_quadratic_rsm_models.rds"))
design_summary <- read.csv(file.path(paths$tables, "design_structure_summary.csv"), stringsAsFactors = FALSE)
cell_counts <- read.csv(file.path(paths$tables, "design_cell_replicate_counts.csv"), stringsAsFactors = FALSE)
fitted_table <- read.csv(file.path(paths$tables, "model_fitted_residuals.csv"), stringsAsFactors = FALSE)

check(identical(unname(pre_hash), unname(post_hash)), "original input files were not modified during Phase 2A")
check(all(toupper(saved_hash$md5) == toupper(unname(post_hash))),
      "input-file checksums match saved checksum table")

check(nrow(reconstructed) == 180L, "reconstructed modelling table has 180 observations")
check(all(response_names %in% names(reconstructed)), "all four response columns are available")
check(all(levels(reconstructed$recipe) == recipe_levels), "recipe factor levels are in the documented order")
check(all(unique(reconstructed$RecipeCode) %in% names(recipe_map)), "recipe code mapping uses only audited recipe codes")
check(all(!is.na(reconstructed$recipe)), "all rows have a mapped recipe")

for (resp in response_names) {
  fit <- models[[resp]]
  check(fit$rank == 18L, paste(resp, "model matrix rank is 18"))
  check(fit$rank - 1L == 17L, paste(resp, "model degrees of freedom is 17"))
  check(df.residual(fit) == 162L, paste(resp, "residual degrees of freedom is 162"))
  fit_rows <- fitted_table[fitted_table$response == resp, ]
  fit_rows <- fit_rows[order(fit_rows$minitab_row_id), ]
  pred <- as.numeric(fitted(fit))
  check(max(abs(pred - fit_rows$fitted)) < 1e-10, paste(resp, "fitted values are reproducible from saved model"))
}

named_design <- setNames(design_summary$value, design_summary$quantity)
check(named_design[["unique_design_conditions"]] == 39, "39 unique design conditions")
check(named_design[["pure_error_degrees_of_freedom"]] == 141, "141 pure-error degrees of freedom")
check(named_design[["lack_of_fit_degrees_of_freedom"]] == 21, "21 lack-of-fit degrees of freedom")
check(all(cell_counts$replicate_count[cell_counts$is_centre] == 12), "centre design cells have 12 replicates")
check(all(cell_counts$replicate_count[!cell_counts$is_centre] == 4), "non-centre design cells have 4 replicates")

model_files <- file.path(paths$models, paste0(tolower(response_names), "_full_quadratic_rsm.rds"))
check(all(file.exists(model_files)), "four response-specific model RDS files were created")
check(file.exists(file.path(paths$tables, "minitab_r_reconciliation.csv")), "CSV reconciliation table was created")
check(file.exists(file.path(paths$tables, "minitab_r_reconciliation.xlsx")), "XLSX reconciliation table was created")
check(file.exists(file.path(paths$tables, "table3_design_space_check.csv")), "Table 3 design-space check was created")
check(file.exists(file.path(paths$reports, "01_model_reconstruction_report.md")), "Markdown report was created")
check(file.exists(file.path(paths$reports, "01_model_reconstruction_report.html")), "HTML report was created")

test_log <- file.path(paths$logs, "test_model_reconstruction.log")
writeLines(c(
  paste("Test run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("R version:", R.version.string),
  if (length(failures) == 0L) "All tests passed." else c("Failures:", failures)
), test_log)

if (length(failures) > 0L) {
  stop(paste(failures, collapse = "\n"))
}

message("All Phase 2A reconstruction tests passed.")
