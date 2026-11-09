## Quick exploratory PCA on the four RSM response variables.
## Uses the same 180-observation reconstructed model data as the RSM work.
## This PCA is exploratory only: it was not used to fit or validate the RSM
## models, and its role is optional supplementary material.

set.seed(20260717)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
RV <- ROOT
OUT <- file.path(RV, "outputs", "pca_exploratory")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

data_path <- file.path(RV, "outputs", "tables", "reconstructed_model_data.csv")
d <- read.csv(data_path, stringsAsFactors = FALSE)
stopifnot(nrow(d) == 180)

response_vars <- c("Coliforms", "Ecoli", "Enterococci", "Methane")
stopifnot(all(response_vars %in% colnames(d)))

pca_input <- d[, response_vars]
stopifnot(all(complete.cases(pca_input)))

## --- PCA: centre and scale all four response variables -------------------
pca <- prcomp(pca_input, center = TRUE, scale. = TRUE)

var_explained <- (pca$sdev^2) / sum(pca$sdev^2)
pc1_pct <- round(100 * var_explained[1], 1)
pc2_pct <- round(100 * var_explained[2], 1)

variance_df <- data.frame(
  PC = paste0("PC", seq_along(pca$sdev)),
  sdev = pca$sdev,
  variance = pca$sdev^2,
  proportion_variance_pct = round(100 * var_explained, 4),
  cumulative_variance_pct = round(100 * cumsum(var_explained), 4)
)
write.csv(variance_df, file.path(OUT, "pca_variance_explained.csv"), row.names = FALSE)

loadings_df <- as.data.frame(pca$rotation)
loadings_df <- cbind(variable = rownames(loadings_df), loadings_df)
rownames(loadings_df) <- NULL
write.csv(loadings_df, file.path(OUT, "pca_loadings.csv"), row.names = FALSE)

scores_df <- data.frame(
  trial_row_id = d$trial_row_id,
  minitab_row_id = d$minitab_row_id,
  design_cell_id = d$design_cell_id,
  Temp = d$Temp,
  OL = d$OL,
  Time = d$Time,
  recipe = d$recipe,
  as.data.frame(pca$x)
)
write.csv(scores_df, file.path(OUT, "pca_scores.csv"), row.names = FALSE)

## --- Plot data prep --------------------------------------------------------
plot_df <- scores_df
plot_df$Temp_f <- factor(plot_df$Temp, levels = sort(unique(plot_df$Temp)),
                          labels = paste0(sort(unique(plot_df$Temp)), " °C"))
plot_df$recipe_f <- factor(plot_df$recipe, levels = c("1S:3FW", "2S:1FW", "3S:1FW"))

temp_group_n <- table(plot_df$Temp_f)
## Ellipses are appropriate here: every temperature level has >= 48
## observations, which is ample for a stable 95% normal-theory ellipse.
ellipse_ok <- all(temp_group_n >= 20)

axis_x <- sprintf("PC1 (%.1f%% variance explained)", pc1_pct)
axis_y <- sprintf("PC2 (%.1f%% variance explained)", pc2_pct)

temp_colours <- c("#3B4CC0", "#8C8C8C", "#B40426")
names(temp_colours) <- levels(plot_df$Temp_f)

theme_pca <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 13)
  )

## --- S6: PCA scores plot ---------------------------------------------------
p_scores <- ggplot(plot_df, aes(x = PC1, y = PC2)) +
  { if (ellipse_ok) stat_ellipse(aes(color = Temp_f, group = Temp_f),
                                  type = "norm", level = 0.95,
                                  linetype = 2, linewidth = 0.6, show.legend = FALSE) } +
  geom_point(aes(color = Temp_f, shape = recipe_f), size = 2.6, alpha = 0.85) +
  scale_color_manual(values = temp_colours, name = "Temperature") +
  scale_shape_manual(values = c(16, 17, 15), name = "Feedstock recipe") +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  labs(
    title = "PCA scores: Coliforms, E. coli, Enterococci and Methane",
    x = axis_x, y = axis_y
  ) +
  theme_pca

ggsave(file.path(OUT, "S6_PCA_scores_plot.png"), p_scores,
       width = 7.2, height = 5.6, dpi = 600, bg = "white")

## --- S7: PCA loadings plot --------------------------------------------------
loadings_plot_df <- loadings_df[, c("variable", "PC1", "PC2")]
loadings_plot_df$label <- loadings_plot_df$variable
loadings_plot_df$label[loadings_plot_df$variable == "Ecoli"] <- "E. coli"
## Anchor labels just beyond the arrow tip, and justify text so it extends
## further away from the origin rather than being centred on the anchor -
## this keeps label text off the arrow line itself.
loadings_plot_df$lx <- loadings_plot_df$PC1 * 1.10
loadings_plot_df$ly <- loadings_plot_df$PC2 * 1.10
loadings_plot_df$hjust <- ifelse(loadings_plot_df$PC1 >= 0, 0, 1)
loadings_plot_df$vjust <- ifelse(loadings_plot_df$PC2 >= 0, 0, 1)

loading_scale <- 1.0
p_loadings <- ggplot(loadings_plot_df, aes(x = PC1 * loading_scale, y = PC2 * loading_scale)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_segment(aes(x = 0, y = 0, xend = PC1 * loading_scale, yend = PC2 * loading_scale),
               arrow = arrow(length = unit(0.22, "cm")), colour = "#1F4D78", linewidth = 0.7) +
  ggrepel::geom_text_repel(aes(x = lx, y = ly, label = label, hjust = hjust, vjust = vjust),
                            size = 4.2, fontface = "bold",
                            colour = "#1F4D78", seed = 20260717,
                            box.padding = 0.4, point.padding = 0.2, min.segment.length = 100) +
  coord_equal() +
  xlim(-1.15, 1.15) + ylim(-1.15, 1.15) +
  labs(
    title = "PCA loadings: Coliforms, E. coli, Enterococci and Methane",
    x = axis_x, y = axis_y
  ) +
  theme_pca

ggsave(file.path(OUT, "S7_PCA_loadings_plot.png"), p_loadings,
       width = 6.4, height = 6.0, dpi = 600, bg = "white")

## --- S8: PCA biplot ---------------------------------------------------------
score_range <- max(abs(c(plot_df$PC1, plot_df$PC2)))
loading_range <- max(abs(c(loadings_df$PC1, loadings_df$PC2)))
biplot_scale <- 0.85 * score_range / loading_range

biplot_loadings <- loadings_df[, c("variable", "PC1", "PC2")]
biplot_loadings$label <- biplot_loadings$variable
biplot_loadings$label[biplot_loadings$variable == "Ecoli"] <- "E. coli"
biplot_loadings$xend <- biplot_loadings$PC1 * biplot_scale
biplot_loadings$yend <- biplot_loadings$PC2 * biplot_scale
## Anchor labels just beyond the arrow tip and justify text outward (away
## from the origin) so it does not sit on top of the arrow line.
biplot_loadings$lx <- biplot_loadings$xend * 1.10
biplot_loadings$ly <- biplot_loadings$yend * 1.10
biplot_loadings$hjust <- ifelse(biplot_loadings$PC1 >= 0, 0, 1)
biplot_loadings$vjust <- ifelse(biplot_loadings$PC2 >= 0, 0, 1)

p_biplot <- ggplot() +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_point(data = plot_df, aes(x = PC1, y = PC2, color = Temp_f, shape = recipe_f),
             size = 2.4, alpha = 0.8) +
  geom_segment(data = biplot_loadings,
               aes(x = 0, y = 0, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.22, "cm")), colour = "#B02418", linewidth = 0.8) +
  ggrepel::geom_text_repel(data = biplot_loadings, aes(x = lx, y = ly, label = label, hjust = hjust, vjust = vjust),
                            size = 4.0, fontface = "bold", colour = "#B02418", seed = 20260717,
                            box.padding = 0.4, point.padding = 0.2, min.segment.length = 100) +
  scale_color_manual(values = temp_colours, name = "Temperature") +
  scale_shape_manual(values = c(16, 17, 15), name = "Feedstock recipe") +
  labs(
    title = "PCA biplot: sample scores and response loadings",
    x = axis_x, y = axis_y
  ) +
  theme_pca

ggsave(file.path(OUT, "S8_PCA_biplot.png"), p_biplot,
       width = 7.4, height = 5.8, dpi = 600, bg = "white")

## --- Console summary for downstream reporting ------------------------------
cat("PCA quick exploratory analysis complete.\n")
cat(sprintf("PC1 variance explained: %.2f%%\n", 100 * var_explained[1]))
cat(sprintf("PC2 variance explained: %.2f%%\n", 100 * var_explained[2]))
cat("Loadings (PC1, PC2):\n")
print(loadings_df[, c("variable", "PC1", "PC2")])
cat("Ellipse drawn per temperature group:", ellipse_ok, "\n")
cat("Files written to:", OUT, "\n")
