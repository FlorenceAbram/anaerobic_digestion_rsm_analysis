set.seed(20260717)

root_dir <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
rv_dir <- root_dir

paths <- list(
  minitab_data = file.path(root_dir, "data", "input", "analysis_design.csv"),
  trial_data = file.path(root_dir, "data", "input", "trial_data.csv"),
  cleaned_data = file.path(root_dir, "data", "input", "combined_data.csv"),
  models = file.path(rv_dir, "models"),
  tables = file.path(rv_dir, "outputs", "tables"),
  figures = file.path(rv_dir, "outputs", "figures"),
  diagnostics = file.path(rv_dir, "outputs", "diagnostics"),
  logs = file.path(rv_dir, "outputs", "logs"),
  reports = file.path(rv_dir, "reports")
)

required_dirs <- unlist(paths[c("models", "tables", "figures", "diagnostics", "logs", "reports")])
invisible(lapply(required_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

response_names <- c("Coliforms", "Ecoli", "Enterococci", "Methane")
microbial_names <- c("Coliforms", "Ecoli", "Enterococci")
factor_names <- c("Temp", "OL", "Time")
model_formula <- as.formula(
  "response ~ (temp_c + ol_c + time_c + recipe)^2 + I(temp_c^2) + I(ol_c^2) + I(time_c^2)"
)

recipe_map <- c(FWS31 = "1S:3FW", FWS12 = "2S:1FW", FWS13 = "3S:1FW")
recipe_levels <- c("1S:3FW", "2S:1FW", "3S:1FW")

recipe_contrasts <- function() {
  mat <- matrix(c(-1, 1, 0,
                  -1, 0, 1),
                nrow = 3, ncol = 2,
                dimnames = list(recipe_levels, c("2S:1FW", "3S:1FW")))
  mat
}

minitab_summary <- data.frame(
  response = response_names,
  n_obs = c(180, 180, 180, 180),
  model_df = c(17, 17, 17, 17),
  residual_df = c(162, 162, 162, 162),
  residual_standard_error = c(0.827033, 0.798127, 0.635807, 16.8947),
  r_squared_pct = c(83.87, 84.0, 86.8, 87.8),
  adjusted_r_squared_pct = c(82.17, 82.3, 85.4, 86.6),
  predicted_r_squared_pct = c(79.9, 80.2, 83.5, 84.7),
  stringsAsFactors = FALSE
)

reported_coefficients <- list(
  Coliforms = c("(Intercept)" = 4.65, temp_c = -2.0607, ol_c = 0.4601, time_c = -0.4842,
                "recipe2S:1FW" = -0.4217, "recipe3S:1FW" = 0.0082,
                "I(temp_c^2)" = -1.424, "I(ol_c^2)" = 0.015, "I(time_c^2)" = -0.116,
                "temp_c:ol_c" = -0.217, "temp_c:time_c" = -0.07,
                "temp_c:recipe2S:1FW" = 0.291, "temp_c:recipe3S:1FW" = -0.056,
                "ol_c:time_c" = -0.222, "ol_c:recipe2S:1FW" = -0.061,
                "ol_c:recipe3S:1FW" = 0.021, "time_c:recipe2S:1FW" = 0.061,
                "time_c:recipe3S:1FW" = -0.023),
  Ecoli = c("(Intercept)" = 4.574, temp_c = -1.918, ol_c = 0.4816, time_c = -0.5467,
            "recipe2S:1FW" = -0.3619, "recipe3S:1FW" = -0.0432,
            "I(temp_c^2)" = -1.556, "I(ol_c^2)" = -0.021, "I(time_c^2)" = -0.117,
            "temp_c:ol_c" = -0.359, "temp_c:time_c" = -0.012,
            "temp_c:recipe2S:1FW" = 0.248, "temp_c:recipe3S:1FW" = -0.068,
            "ol_c:time_c" = -0.175, "ol_c:recipe2S:1FW" = -0.058,
            "ol_c:recipe3S:1FW" = -0.028, "time_c:recipe2S:1FW" = 0.111,
            "time_c:recipe3S:1FW" = -0.058),
  Enterococci = c("(Intercept)" = 5.676, temp_c = -1.8134, ol_c = 0.296, time_c = -0.0255,
                  "recipe2S:1FW" = -0.3993, "recipe3S:1FW" = 0.2194,
                  "I(temp_c^2)" = -1.3168, "I(ol_c^2)" = -0.4455, "I(time_c^2)" = 0.0007,
                  "temp_c:ol_c" = 0.0335, "temp_c:time_c" = -0.2702,
                  "temp_c:recipe2S:1FW" = 0.2363, "temp_c:recipe3S:1FW" = -0.0225,
                  "ol_c:time_c" = 0.0758, "ol_c:recipe2S:1FW" = -0.1439,
                  "ol_c:recipe3S:1FW" = 0.1414, "time_c:recipe2S:1FW" = -0.1039,
                  "time_c:recipe3S:1FW" = 0.0953),
  Methane = c("(Intercept)" = 92.38, temp_c = -2.85, ol_c = 16.19, time_c = 23.06,
              "recipe2S:1FW" = -2.46, "recipe3S:1FW" = -1.35,
              "I(temp_c^2)" = -74.04, "I(ol_c^2)" = -5.72, "I(time_c^2)" = -5.3,
              "temp_c:ol_c" = -0.32, "temp_c:time_c" = -1.15,
              "temp_c:recipe2S:1FW" = 1.75, "temp_c:recipe3S:1FW" = -0.99,
              "ol_c:time_c" = 14.58, "ol_c:recipe2S:1FW" = -5.67,
              "ol_c:recipe3S:1FW" = 2.98, "time_c:recipe2S:1FW" = 4.23,
              "time_c:recipe3S:1FW" = -2.86)
)

reported_equation_predict <- function(dat, response) {
  T <- dat$Temp
  O <- dat$OL
  H <- dat$Time
  R <- as.character(dat$recipe)
  out <- rep(NA_real_, nrow(dat))
  if (response == "Coliforms") {
    out[R == "2S:1FW"] <- 0.25 + 0.2509*T[R == "2S:1FW"] + 0.833*O[R == "2S:1FW"] + 0.160*H[R == "2S:1FW"] - 0.004396*T[R == "2S:1FW"]^2 + 0.0067*O[R == "2S:1FW"]^2 - 0.0129*H[R == "2S:1FW"]^2 - 0.00805*T[R == "2S:1FW"]*O[R == "2S:1FW"] - 0.00130*T[R == "2S:1FW"]*H[R == "2S:1FW"] - 0.0494*O[R == "2S:1FW"]*H[R == "2S:1FW"]
    out[R == "3S:1FW"] <- 1.45 + 0.2316*T[R == "3S:1FW"] + 0.888*O[R == "3S:1FW"] + 0.132*H[R == "3S:1FW"] - 0.004396*T[R == "3S:1FW"]^2 + 0.0067*O[R == "3S:1FW"]^2 - 0.0129*H[R == "3S:1FW"]^2 - 0.00805*T[R == "3S:1FW"]*O[R == "3S:1FW"] - 0.00130*T[R == "3S:1FW"]*H[R == "3S:1FW"] - 0.0494*O[R == "3S:1FW"]*H[R == "3S:1FW"]
    out[R == "1S:3FW"] <- 2.23 + 0.2217*T[R == "1S:3FW"] + 0.901*O[R == "1S:3FW"] + 0.127*H[R == "1S:3FW"] - 0.004396*T[R == "1S:3FW"]^2 + 0.0067*O[R == "1S:3FW"]^2 - 0.0129*H[R == "1S:3FW"]^2 - 0.00805*T[R == "1S:3FW"]*O[R == "1S:3FW"] - 0.00130*T[R == "1S:3FW"]*H[R == "1S:3FW"] - 0.0494*O[R == "1S:3FW"]*H[R == "1S:3FW"]
  }
  if (response == "Ecoli") {
    out[R == "2S:1FW"] <- -0.63 + 0.2904*T[R == "2S:1FW"] + 1.045*O[R == "2S:1FW"] + 0.096*H[R == "2S:1FW"] - 0.004801*T[R == "2S:1FW"]^2 - 0.0094*O[R == "2S:1FW"]^2 - 0.0129*H[R == "2S:1FW"]^2 - 0.01330*T[R == "2S:1FW"]*O[R == "2S:1FW"] - 0.00022*T[R == "2S:1FW"]*H[R == "2S:1FW"] - 0.0388*O[R == "2S:1FW"]*H[R == "2S:1FW"]
    out[R == "3S:1FW"] <- 0.64 + 0.2729*T[R == "3S:1FW"] + 1.065*O[R == "3S:1FW"] + 0.040*H[R == "3S:1FW"] - 0.004801*T[R == "3S:1FW"]^2 - 0.0094*O[R == "3S:1FW"]^2 - 0.0129*H[R == "3S:1FW"]^2 - 0.01330*T[R == "3S:1FW"]*O[R == "3S:1FW"] - 0.00022*T[R == "3S:1FW"]*H[R == "3S:1FW"] - 0.0388*O[R == "3S:1FW"]*H[R == "3S:1FW"]
    out[R == "1S:3FW"] <- 1.16 + 0.2666*T[R == "1S:3FW"] + 1.140*O[R == "1S:3FW"] + 0.041*H[R == "1S:3FW"] - 0.004801*T[R == "1S:3FW"]^2 - 0.0094*O[R == "1S:3FW"]^2 - 0.0129*H[R == "1S:3FW"]^2 - 0.01330*T[R == "1S:3FW"]*O[R == "1S:3FW"] - 0.00022*T[R == "1S:3FW"]*H[R == "1S:3FW"] - 0.0388*O[R == "1S:3FW"]*H[R == "1S:3FW"]
  }
  if (response == "Enterococci") {
    out[R == "2S:1FW"] <- 1.405 + 0.2407*T[R == "2S:1FW"] + 0.747*O[R == "2S:1FW"] + 0.107*H[R == "2S:1FW"] - 0.004064*T[R == "2S:1FW"]^2 - 0.1980*O[R == "2S:1FW"]^2 + 0.0001*H[R == "2S:1FW"]^2 + 0.00124*T[R == "2S:1FW"]*O[R == "2S:1FW"] - 0.00500*T[R == "2S:1FW"]*H[R == "2S:1FW"] + 0.0168*O[R == "2S:1FW"]*H[R == "2S:1FW"]
    out[R == "3S:1FW"] <- 1.776 + 0.2263*T[R == "3S:1FW"] + 0.937*O[R == "3S:1FW"] + 0.174*H[R == "3S:1FW"] - 0.004064*T[R == "3S:1FW"]^2 - 0.1980*O[R == "3S:1FW"]^2 + 0.0001*H[R == "3S:1FW"]^2 + 0.00124*T[R == "3S:1FW"]*O[R == "3S:1FW"] - 0.00500*T[R == "3S:1FW"]*H[R == "3S:1FW"] + 0.0168*O[R == "3S:1FW"]*H[R == "3S:1FW"]
    out[R == "1S:3FW"] <- 2.489 + 0.2157*T[R == "1S:3FW"] + 0.844*O[R == "1S:3FW"] + 0.145*H[R == "1S:3FW"] - 0.004064*T[R == "1S:3FW"]^2 - 0.1980*O[R == "1S:3FW"]^2 + 0.0001*H[R == "1S:3FW"]^2 + 0.00124*T[R == "1S:3FW"]*O[R == "1S:3FW"] - 0.00500*T[R == "1S:3FW"]*H[R == "1S:3FW"] + 0.0168*O[R == "1S:3FW"]*H[R == "1S:3FW"]
  }
  if (response == "Methane") {
    out[R == "2S:1FW"] <- -287.3 + 17.000*T[R == "2S:1FW"] - 1.82*O[R == "2S:1FW"] + 10.47*H[R == "2S:1FW"] - 0.22851*T[R == "2S:1FW"]^2 - 2.54*O[R == "2S:1FW"]^2 - 0.589*H[R == "2S:1FW"]^2 - 0.0120*T[R == "2S:1FW"]*O[R == "2S:1FW"] - 0.0212*T[R == "2S:1FW"]*H[R == "2S:1FW"] + 3.240*O[R == "2S:1FW"]*H[R == "2S:1FW"]
    out[R == "3S:1FW"] <- -277.9 + 16.847*T[R == "3S:1FW"] + 3.95*O[R == "3S:1FW"] + 8.11*H[R == "3S:1FW"] - 0.22851*T[R == "3S:1FW"]^2 - 2.54*O[R == "3S:1FW"]^2 - 0.589*H[R == "3S:1FW"]^2 - 0.0120*T[R == "3S:1FW"]*O[R == "3S:1FW"] - 0.0212*T[R == "3S:1FW"]*H[R == "3S:1FW"] + 3.240*O[R == "3S:1FW"]*H[R == "3S:1FW"]
    out[R == "1S:3FW"] <- -275.8 + 16.860*T[R == "1S:3FW"] + 3.76*O[R == "1S:3FW"] + 8.60*H[R == "1S:3FW"] - 0.22851*T[R == "1S:3FW"]^2 - 2.54*O[R == "1S:3FW"]^2 - 0.589*H[R == "1S:3FW"]^2 - 0.0120*T[R == "1S:3FW"]*O[R == "1S:3FW"] - 0.0212*T[R == "1S:3FW"]*H[R == "1S:3FW"] + 3.240*O[R == "1S:3FW"]*H[R == "1S:3FW"]
  }
  out
}

make_key <- function(dat, cols) {
  do.call(paste, c(dat[cols], sep = "|"))
}

transform_microbes <- function(dat) {
  for (nm in microbial_names) {
    dat[[nm]] <- pmax(log10(dat[[nm]] + 1), 1)
  }
  dat
}

round_for_minitab <- function(dat) {
  for (nm in microbial_names) {
    dat[[nm]] <- round(dat[[nm]], 2)
  }
  dat$Methane <- round(dat$Methane, 1)
  dat
}

append_recon <- function(rows, response, statistic, minitab, r_value, tolerance, explanation) {
  abs_diff <- if (is.na(minitab) || is.na(r_value)) NA_real_ else abs(r_value - minitab)
  rel_diff <- if (is.na(abs_diff) || isTRUE(all.equal(minitab, 0))) NA_real_ else abs_diff / abs(minitab)
  pass <- if (is.na(abs_diff)) "not_tested" else if (abs_diff <= tolerance) "pass" else "fail"
  rbind(rows, data.frame(
    response = response,
    statistic = statistic,
    minitab_value = minitab,
    r_value = r_value,
    absolute_difference = abs_diff,
    relative_difference = rel_diff,
    tolerance = tolerance,
    pass = pass,
    explanation = explanation,
    stringsAsFactors = FALSE
  ))
}

write_simple_xlsx <- function(dat, path, sheet_name = "reconciliation") {
  tmp <- tempfile("xlsx_")
  dir.create(tmp)
  dir.create(file.path(tmp, "_rels"))
  dir.create(file.path(tmp, "xl"))
  dir.create(file.path(tmp, "xl", "_rels"))
  dir.create(file.path(tmp, "xl", "worksheets"))
  esc <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x
  }
  cell_ref <- function(i, j) {
    letters <- ""
    while (j > 0) {
      r <- (j - 1) %% 26
      letters <- paste0(LETTERS[r + 1], letters)
      j <- (j - r - 1) %/% 26
    }
    paste0(letters, i)
  }
  all_dat <- rbind(names(dat), as.matrix(dat))
  rows <- character(nrow(all_dat))
  for (i in seq_len(nrow(all_dat))) {
    cells <- character(ncol(all_dat))
    for (j in seq_len(ncol(all_dat))) {
      value <- esc(all_dat[i, j])
      cells[j] <- paste0("<c r=\"", cell_ref(i, j), "\" t=\"inlineStr\"><is><t>", value, "</t></is></c>")
    }
    rows[i] <- paste0("<row r=\"", i, "\">", paste(cells, collapse = ""), "</row>")
  }
  sheet <- paste0("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>",
                  "<worksheet xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\">",
                  "<sheetData>", paste(rows, collapse = ""), "</sheetData></worksheet>")
  workbook <- paste0("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>",
                     "<workbook xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\" ",
                     "xmlns:r=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships\">",
                     "<sheets><sheet name=\"", esc(sheet_name), "\" sheetId=\"1\" r:id=\"rId1\"/></sheets></workbook>")
  rels <- "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?><Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\"><Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument\" Target=\"xl/workbook.xml\"/></Relationships>"
  wb_rels <- "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?><Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\"><Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet\" Target=\"worksheets/sheet1.xml\"/></Relationships>"
  content_types <- "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?><Types xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\"><Default Extension=\"rels\" ContentType=\"application/vnd.openxmlformats-package.relationships+xml\"/><Default Extension=\"xml\" ContentType=\"application/xml\"/><Override PartName=\"/xl/workbook.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml\"/><Override PartName=\"/xl/worksheets/sheet1.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml\"/></Types>"
  writeLines(content_types, file.path(tmp, "[Content_Types].xml"), useBytes = TRUE)
  writeLines(rels, file.path(tmp, "_rels", ".rels"), useBytes = TRUE)
  writeLines(workbook, file.path(tmp, "xl", "workbook.xml"), useBytes = TRUE)
  writeLines(wb_rels, file.path(tmp, "xl", "_rels", "workbook.xml.rels"), useBytes = TRUE)
  writeLines(sheet, file.path(tmp, "xl", "worksheets", "sheet1.xml"), useBytes = TRUE)
  if (file.exists(path)) unlink(path)
  old <- setwd(tmp)
  zip_result <- try(utils::zip(zipfile = path, files = list.files(".", recursive = TRUE, all.files = TRUE, no.. = TRUE)),
                    silent = TRUE)
  setwd(old)
  if (inherits(zip_result, "try-error") || !file.exists(path)) {
    ps_quote <- function(x) paste0("'", gsub("'", "''", normalizePath(x, winslash = "\\", mustWork = FALSE)), "'")
    if (file.exists(path)) unlink(path)
    ps <- paste(
      "Add-Type -AssemblyName System.IO.Compression.FileSystem;",
      paste0("if (Test-Path -LiteralPath ", ps_quote(path), ") { Remove-Item -LiteralPath ", ps_quote(path), " -Force };"),
      paste0("[System.IO.Compression.ZipFile]::CreateFromDirectory(", ps_quote(tmp), ", ", ps_quote(path), ");")
    )
    status <- system2("powershell.exe", c("-NoProfile", "-Command", ps), stdout = TRUE, stderr = TRUE)
    if (!file.exists(path)) stop("Failed to create XLSX file: ", paste(status, collapse = "\n"))
  }
  invisible(path)
}

build_reconstructed_data <- function() {
  minitab <- read.csv(paths$minitab_data, check.names = FALSE, stringsAsFactors = FALSE)
  trial_raw <- read.csv(paths$trial_data, check.names = FALSE, stringsAsFactors = FALSE)
  cleaned <- read.csv(paths$cleaned_data, check.names = FALSE, stringsAsFactors = FALSE)

  if (nrow(minitab) != 180L) stop("Minitab modelling file does not contain 180 observations.")
  if (nrow(trial_raw) != 180L) stop("trial_data.csv does not contain the expected 180 non-Day-0 observations.")
  missing_responses <- setdiff(response_names, names(minitab))
  if (length(missing_responses) > 0L) stop("Missing response columns in Minitab file: ", paste(missing_responses, collapse = ", "))
  if (!("Recipe" %in% names(trial_raw))) stop("trial_data.csv has no Recipe column.")

  trial_transformed_full <- transform_microbes(trial_raw)
  trial_transformed <- round_for_minitab(trial_transformed_full)
  minitab_rounded <- round_for_minitab(minitab)
  key_cols <- c("Temp", "OL", "Time", "Ammonia", "pH", "sCOD", "TS", "VS", response_names)
  minitab_key <- make_key(minitab_rounded, key_cols)
  trial_key <- make_key(trial_transformed, key_cols)
  if (anyDuplicated(minitab_key) > 0L || anyDuplicated(trial_key) > 0L) {
    stop("Row matching key is not unique; manual reconciliation is required.")
  }
  matched_idx <- match(minitab_key, trial_key)
  if (any(is.na(matched_idx))) {
    write.csv(data.frame(minitab_row_id = which(is.na(matched_idx)), key = minitab_key[is.na(matched_idx)]),
              file.path(paths$tables, "unmatched_minitab_rows.csv"), row.names = FALSE)
    stop("Unexplained row-level discrepancies prevent model reconstruction.")
  }

  reconstructed <- minitab
  reconstructed$minitab_row_id <- seq_len(nrow(minitab))
  reconstructed$trial_row_id <- matched_idx
  reconstructed$RecipeCode <- trial_raw$Recipe[matched_idx]
  if (any(!(reconstructed$RecipeCode %in% names(recipe_map)))) {
    stop("Unexpected recipe code found: ", paste(unique(reconstructed$RecipeCode), collapse = ", "))
  }
  reconstructed$recipe <- factor(unname(recipe_map[reconstructed$RecipeCode]), levels = recipe_levels)
  reconstructed$temp_c <- (reconstructed$Temp - 37) / 18
  reconstructed$ol_c <- (reconstructed$OL - 2.0) / 1.5
  reconstructed$time_c <- (reconstructed$Time - 6) / 3
  reconstructed$design_cell_id <- paste(reconstructed$Temp, reconstructed$OL, reconstructed$Time, reconstructed$recipe, sep = "__")

  row_comp <- data.frame(
    minitab_row_id = reconstructed$minitab_row_id,
    trial_row_id = reconstructed$trial_row_id,
    temperature_diff = reconstructed$Temp - trial_raw$Temp[matched_idx],
    olr_diff = reconstructed$OL - trial_raw$OL[matched_idx],
    retention_time_diff = reconstructed$Time - trial_raw$Time[matched_idx],
    recipe_code = reconstructed$RecipeCode,
    recipe_label = as.character(reconstructed$recipe),
    coliform_round_diff = reconstructed$Coliforms - round(trial_transformed_full$Coliforms[matched_idx], 2),
    ecoli_round_diff = reconstructed$Ecoli - round(trial_transformed_full$Ecoli[matched_idx], 2),
    enterococci_round_diff = reconstructed$Enterococci - round(trial_transformed_full$Enterococci[matched_idx], 2),
    methane_diff = reconstructed$Methane - round(trial_raw$Methane[matched_idx], 1),
    stringsAsFactors = FALSE
  )
  write.csv(row_comp, file.path(paths$tables, "row_level_data_comparison.csv"), row.names = FALSE)

  comp_summary <- data.frame(
    field = c("temperature", "organic_loading_rate", "retention_time", "recipe", "coliforms_transformed_rounded",
              "ecoli_transformed_rounded", "enterococci_transformed_rounded", "methane"),
    max_abs_difference = c(
      max(abs(row_comp$temperature_diff)),
      max(abs(row_comp$olr_diff)),
      max(abs(row_comp$retention_time_diff)),
      NA_real_,
      max(abs(row_comp$coliform_round_diff)),
      max(abs(row_comp$ecoli_round_diff)),
      max(abs(row_comp$enterococci_round_diff)),
      max(abs(row_comp$methane_diff))
    ),
    discrepancies = c(
      sum(abs(row_comp$temperature_diff) > 0),
      sum(abs(row_comp$olr_diff) > 0),
      sum(abs(row_comp$retention_time_diff) > 0),
      sum(is.na(row_comp$recipe_label)),
      sum(abs(row_comp$coliform_round_diff) > 1e-12),
      sum(abs(row_comp$ecoli_round_diff) > 1e-12),
      sum(abs(row_comp$enterococci_round_diff) > 1e-12),
      sum(abs(row_comp$methane_diff) > 1e-12)
    ),
    stringsAsFactors = FALSE
  )
  write.csv(comp_summary, file.path(paths$tables, "data_reconstruction_checks.csv"), row.names = FALSE)

  if (any(comp_summary$discrepancies > 0, na.rm = TRUE)) {
    stop("Unexplained row-level discrepancies found; see data_reconstruction_checks.csv.")
  }

  coded_levels <- data.frame(
    predictor = rep(c("temp_c", "ol_c", "time_c"), each = 3),
    raw_level = c(sort(unique(reconstructed$Temp)), sort(unique(reconstructed$OL)), sort(unique(reconstructed$Time))),
    coded_level = c(sort(unique(reconstructed$temp_c)), sort(unique(reconstructed$ol_c)), sort(unique(reconstructed$time_c)))
  )
  write.csv(coded_levels, file.path(paths$tables, "coded_predictor_levels.csv"), row.names = FALSE)
  if (any(abs(coded_levels$coded_level - rep(c(-1, 0, 1), 3)) > 1e-12)) {
    stop("Coded predictor levels are not the expected -1, 0, +1.")
  }

  cell_counts <- as.data.frame(table(reconstructed$design_cell_id), stringsAsFactors = FALSE)
  names(cell_counts) <- c("design_cell_id", "replicate_count")
  cell_parts <- do.call(rbind, strsplit(cell_counts$design_cell_id, "__", fixed = TRUE))
  cell_counts$Temp <- as.numeric(cell_parts[, 1])
  cell_counts$OL <- as.numeric(cell_parts[, 2])
  cell_counts$Time <- as.numeric(cell_parts[, 3])
  cell_counts$recipe <- cell_parts[, 4]
  cell_counts$is_centre <- cell_counts$Temp == 37 & cell_counts$OL == 2 & cell_counts$Time == 6
  write.csv(cell_counts, file.path(paths$tables, "design_cell_replicate_counts.csv"), row.names = FALSE)

  design_summary <- data.frame(
    quantity = c("observations", "unique_design_conditions", "noncentre_design_conditions",
                 "centre_design_conditions", "centre_observations", "model_degrees_of_freedom",
                 "residual_degrees_of_freedom", "pure_error_degrees_of_freedom", "lack_of_fit_degrees_of_freedom"),
    value = c(nrow(reconstructed), length(unique(reconstructed$design_cell_id)),
              sum(!cell_counts$is_centre), sum(cell_counts$is_centre),
              sum(cell_counts$replicate_count[cell_counts$is_centre]), NA, NA,
              nrow(reconstructed) - length(unique(reconstructed$design_cell_id)), NA),
    stringsAsFactors = FALSE
  )

  saveRDS(reconstructed, file.path(paths$models, "reconstructed_model_data.rds"))
  write.csv(reconstructed, file.path(paths$tables, "reconstructed_model_data.csv"), row.names = FALSE)
  saveRDS(list(minitab = minitab, trial_raw = trial_raw, trial_transformed = trial_transformed_full, cleaned = cleaned),
          file.path(paths$models, "input_data_snapshots.rds"))

  list(data = reconstructed, design_summary = design_summary, cell_counts = cell_counts)
}

fit_models <- function(reconstructed, design_summary) {
  models <- list()
  model_summaries <- data.frame()
  coefficient_rows <- data.frame()
  fitted_rows <- data.frame()
  reconciliation <- data.frame()

  for (resp in response_names) {
    dat <- reconstructed
    dat$response <- dat[[resp]]
    contrasts(dat$recipe) <- recipe_contrasts()
    fit <- lm(model_formula, data = dat)
    models[[resp]] <- fit
    saveRDS(fit, file.path(paths$models, paste0(tolower(resp), "_full_quadratic_rsm.rds")))

    summ <- summary(fit)
    press <- sum((residuals(fit) / (1 - hatvalues(fit)))^2)
    pred_r2 <- 1 - press / sum((dat$response - mean(dat$response))^2)
    model_stats <- data.frame(
      response = resp,
      observations = nobs(fit),
      model_df = fit$rank - 1,
      residual_df = df.residual(fit),
      model_matrix_rank = fit$rank,
      residual_standard_error = summ$sigma,
      r_squared_pct = 100 * summ$r.squared,
      adjusted_r_squared_pct = 100 * summ$adj.r.squared,
      PRESS = press,
      predicted_r_squared_pct = 100 * pred_r2,
      stringsAsFactors = FALSE
    )
    model_summaries <- rbind(model_summaries, model_stats)

    coef_table <- as.data.frame(coef(summ))
    coef_table$term <- row.names(coef_table)
    row.names(coef_table) <- NULL
    names(coef_table)[1:4] <- c("estimate", "std_error", "t_value", "p_value")
    coef_table$response <- resp
    coefficient_rows <- rbind(coefficient_rows, coef_table[, c("response", "term", "estimate", "std_error", "t_value", "p_value")])

    eq_fit <- reported_equation_predict(dat, resp)
    fitted_rows <- rbind(fitted_rows, data.frame(
      response = resp,
      minitab_row_id = dat$minitab_row_id,
      design_cell_id = dat$design_cell_id,
      observed = dat$response,
      fitted = fitted(fit),
      residual = residuals(fit),
      reported_rounded_equation_fitted = eq_fit,
      r_minus_reported_rounded_equation = fitted(fit) - eq_fit,
      stringsAsFactors = FALSE
    ))

    ms <- minitab_summary[minitab_summary$response == resp, ]
    reconciliation <- append_recon(reconciliation, resp, "number_of_observations", ms$n_obs, nobs(fit), 0,
                                   "Exact count reported by Minitab Word output.")
    reconciliation <- append_recon(reconciliation, resp, "model_degrees_of_freedom", ms$model_df, fit$rank - 1, 0,
                                   "Model DF calculated from model-matrix rank minus intercept.")
    reconciliation <- append_recon(reconciliation, resp, "residual_degrees_of_freedom", ms$residual_df, df.residual(fit), 0,
                                   "Residual DF calculated from n minus model-matrix rank.")
    reconciliation <- append_recon(reconciliation, resp, "model_matrix_rank", 18, fit$rank, 0,
                                   "Full-rank model matrix should contain 18 parameters including intercept.")
    reconciliation <- append_recon(reconciliation, resp, "residual_standard_error", ms$residual_standard_error, summ$sigma,
                                   if (resp == "Methane") 0.005 else 0.0005,
                                   "Tolerance allows final displayed precision in Word report.")
    reconciliation <- append_recon(reconciliation, resp, "r_squared_pct", ms$r_squared_pct, 100 * summ$r.squared, 0.06,
                                   "Tolerance allows percent rounding in Word report.")
    reconciliation <- append_recon(reconciliation, resp, "adjusted_r_squared_pct", ms$adjusted_r_squared_pct, 100 * summ$adj.r.squared, 0.06,
                                   "Tolerance allows percent rounding in Word report.")
    reconciliation <- append_recon(reconciliation, resp, "PRESS", NA_real_, press, NA_real_,
                                   "PRESS is not displayed in the available Word/Minitab output; R value recorded.")
    reconciliation <- append_recon(reconciliation, resp, "predicted_r_squared_pct", ms$predicted_r_squared_pct, 100 * pred_r2, 0.08,
                                   "Minitab predicted R-sq is PRESS based; tolerance allows displayed percent rounding.")
    max_eq_diff <- max(abs(fitted(fit) - eq_fit), na.rm = TRUE)
    reconciliation <- append_recon(reconciliation, resp, "fitted_values_vs_reported_rounded_equations_max_abs_diff",
                                   0, max_eq_diff, if (resp == "Methane") 0.35 else 0.04,
                                   "Full-precision Minitab fitted values are unavailable; this checks consistency against rounded printed equations.")
    reconciliation <- append_recon(reconciliation, resp, "residuals_available_for_full_precision_minitab_comparison",
                                   NA_real_, NA_real_, NA_real_,
                                   "Full-precision Minitab fitted values/residuals were not present in readable outputs.")

    rc <- reported_coefficients[[resp]]
    est <- coef(fit)
    for (term in names(rc)) {
      tolerance <- if (resp == "Methane") 0.006 else 0.0008
      if (term %in% c("(Intercept)", "I(temp_c^2)", "I(ol_c^2)", "I(time_c^2)") && resp != "Methane") tolerance <- 0.006
      reconciliation <- append_recon(reconciliation, resp, paste0("coefficient_", term), rc[[term]], est[[term]], tolerance,
                                     "Coefficient compared with displayed coded coefficient; tolerance follows displayed precision.")
    }
  }

  common_rank <- unique(model_summaries$model_matrix_rank)
  if (length(common_rank) != 1L) stop("Model matrix rank differs across responses.")
  design_summary$value[design_summary$quantity == "model_degrees_of_freedom"] <- common_rank - 1
  design_summary$value[design_summary$quantity == "residual_degrees_of_freedom"] <- nrow(reconstructed) - common_rank
  design_summary$value[design_summary$quantity == "lack_of_fit_degrees_of_freedom"] <-
    length(unique(reconstructed$design_cell_id)) - common_rank

  write.csv(design_summary, file.path(paths$tables, "design_structure_summary.csv"), row.names = FALSE)
  write.csv(model_summaries, file.path(paths$tables, "model_summary_r.csv"), row.names = FALSE)
  write.csv(coefficient_rows, file.path(paths$tables, "model_coefficients_r.csv"), row.names = FALSE)
  write.csv(fitted_rows, file.path(paths$tables, "model_fitted_residuals.csv"), row.names = FALSE)
  write.csv(reconciliation, file.path(paths$tables, "minitab_r_reconciliation.csv"), row.names = FALSE)
  write_simple_xlsx(reconciliation, file.path(paths$tables, "minitab_r_reconciliation.xlsx"))
  saveRDS(models, file.path(paths$models, "all_full_quadratic_rsm_models.rds"))

  list(models = models, model_summaries = model_summaries, reconciliation = reconciliation, design_summary = design_summary)
}

phase2a_reconstruct_main <- function() {
  log_path <- file.path(paths$logs, "01_reconstruct_rsm_models.log")
  con <- file(log_path, open = "wt")
  sink(con, split = TRUE)
  sink(con, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(con)
  }, add = TRUE)

  cat("Phase 2A reconstruction started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
  cat("Seed: 20260717\n")
  cat("R version:", R.version.string, "\n\n")

  input_files <- unlist(paths[c("minitab_data", "trial_data", "cleaned_data")])
  checksum_rows <- data.frame(
    file = names(input_files),
    path = file.path("data", "input", basename(unname(input_files))),
    md5 = unname(tools::md5sum(input_files)),
    stringsAsFactors = FALSE
  )
  write.csv(checksum_rows, file.path(paths$tables, "input_file_checksums.csv"), row.names = FALSE)

  reconstructed <- build_reconstructed_data()
  fit_result <- fit_models(reconstructed$data, reconstructed$design_summary)

  sink_path <- file.path(paths$logs, "session_info.txt")
  capture.output({
    cat("R version\n")
    cat(R.version.string, "\n\n")
    cat("Package versions\n")
    print(installed.packages()[, c("Package", "Version")])
    cat("\nSession info\n")
    print(sessionInfo())
  }, file = sink_path)

  warnings_path <- file.path(paths$logs, "warnings_01_reconstruct_rsm_models.txt")
  capture.output(warnings(), file = warnings_path)

  expected <- c(
    observations = 180,
    unique_design_conditions = 39,
    centre_design_conditions = 3,
    centre_observations = 36,
    model_degrees_of_freedom = 17,
    residual_degrees_of_freedom = 162,
    pure_error_degrees_of_freedom = 141,
    lack_of_fit_degrees_of_freedom = 21
  )
  ds <- fit_result$design_summary
  got <- setNames(ds$value, ds$quantity)
  bad <- names(expected)[abs(got[names(expected)] - expected) > 1e-12]
  if (length(bad) > 0L) stop("Calculated design/model structure differs for: ", paste(bad, collapse = ", "))

  cat("Phase 2A reconstruction completed successfully.\n")
  invisible(fit_result)
}

if (sys.nframe() == 0L) {
  phase2a_reconstruct_main()
}
