## =============================================================================
## v_plants — next steps after variable screening
## Author:  Valtýr  (script assembled with Claude)
## Purpose: Ordination visualisation, indicator-species analysis, and a
##          univariate follow-up for the strongest environmental drivers
##          identified in v_plants_variable_screening.html
##
## Assumes the following objects are already in your environment (source
## `v_plants_prep.R` first; `v_plants_variable_screening.Rmd` is an audit trail
## describing how these were chosen, but does not need to be sourced):
##   community_long_paired  – long-format cover table (plot × sample × taxon)
##   community_wide_strict  – wide-format community matrix with plot/sample cols
##   meta_strict            – metadata aligned to community_wide_strict
##   mat_strict             – numeric species matrix (rows = nrow(mat_strict))
##   nmds_strict            – metaMDS result on mat_strict
##   fit_strict             – envfit result on nmds_strict
##   vars_driver_strict_final, driver_final – variable vectors from prep/screening
## =============================================================================

library(tidyverse)
library(vegan)
library(broom)
library(ggrepel)    # install.packages("ggrepel")  if missing
library(indicspecies)  # install.packages("indicspecies")  if missing

set.seed(123)


## ── SECTION 1 ─────────────────────────────────────────────────────────────────
## NMDS ordination plot with environmental vectors
## Purpose: visualise temporal shift in community composition and the drivers
##          that explain most variation in ordination space.
## ─────────────────────────────────────────────────────────────────────────────

# 1a. Extract site scores
nmds_scores <- as.data.frame(scores(nmds_strict, display = "sites")) %>%
  bind_cols(meta_strict %>% select(plot_number, sample))

# 1b. Extract significant vector/factor arrows from envfit
sig_thresh <- 0.05

# continuous vectors: robustly map first two ordination axes to NMDS1/NMDS2
vec_df <- tibble(variable = character(), NMDS1 = numeric(), NMDS2 = numeric(),
                 r2 = numeric(), p_value = numeric(), kind = character())
if (!is.null(fit_strict$vectors) && !is.null(fit_strict$vectors$arrows)) {
  vec_raw <- as.data.frame(fit_strict$vectors$arrows) %>%
    rownames_to_column("variable")
  vec_axes <- setdiff(names(vec_raw), "variable")
  if (length(vec_axes) >= 2) {
    vec_df <- vec_raw %>%
      rename(NMDS1 = all_of(vec_axes[1]), NMDS2 = all_of(vec_axes[2])) %>%
      select(variable, NMDS1, NMDS2) %>%
      mutate(
        r2 = as.numeric(fit_strict$vectors$r[variable]),
        p_value = as.numeric(fit_strict$vectors$pvals[variable]),
        kind = "numeric"
      )
  }
}

# factor centroids: robustly map first two ordination axes to NMDS1/NMDS2
fac_df <- tibble(variable = character(), NMDS1 = numeric(), NMDS2 = numeric(),
                 r2 = numeric(), p_value = numeric(), kind = character())
if (!is.null(fit_strict$factors) && !is.null(fit_strict$factors$centroids)) {
  fac_raw <- as.data.frame(fit_strict$factors$centroids) %>%
    rownames_to_column("variable")
  fac_axes <- setdiff(names(fac_raw), "variable")
  if (length(fac_axes) >= 2) {
    fac_df <- fac_raw %>%
      rename(NMDS1 = all_of(fac_axes[1]), NMDS2 = all_of(fac_axes[2])) %>%
      select(variable, NMDS1, NMDS2) %>%
      mutate(
        r2 = NA_real_,
        p_value = NA_real_,
        kind = "factor"
      )
  }
}

arrow_df <- bind_rows(vec_df, fac_df) %>%
  filter(is.na(p_value) | p_value <= sig_thresh)

# Scaling arrows for display (multiply by 0.7 × ordination radius)
ordi_radius <- max(abs(c(nmds_scores$NMDS1, nmds_scores$NMDS2))) * 0.7

arrow_df <- arrow_df %>%
  mutate(
    xend = NMDS1 * ordi_radius,
    yend = NMDS2 * ordi_radius
  )

# 1c. Build the plot
p_nmds <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  # paired sites connected by thin lines
  geom_line(
    aes(group = plot_number),
    colour = "grey70", linewidth = 0.3, alpha = 0.6
  ) +
  # site points
  geom_point(
    aes(colour = sample, shape = sample),
    size = 2.2, alpha = 0.85
  ) +
  # environmental arrows
  geom_segment(
    data = arrow_df %>% filter(kind == "numeric"),
    aes(x = 0, y = 0, xend = xend, yend = yend),
    arrow = arrow(length = unit(0.2, "cm")),
    colour = "#c0392b", linewidth = 0.7
  ) +
  geom_label(
    data = arrow_df %>% filter(kind == "numeric"),
    aes(x = xend * 1.1, y = yend * 1.1, label = variable),
    size = 2.8, colour = "#c0392b", fill = alpha("white", 0.7),
    label.size = 0
  ) +
  scale_colour_manual(
    values = c(first = "#2980b9", last = "#e67e22"),
    labels = c(first = "First survey", last = "Last survey")
  ) +
  scale_shape_manual(
    values = c(first = 16, last = 17),
    labels = c(first = "First survey", last = "Last survey")
  ) +
  labs(
    title   = "NMDS of plant community composition (Bray–Curtis)",
    subtitle = sprintf(
      "%d revisited plots · stress = %.3f · paired PERMANOVA R² = 2.7%%, p = 0.001",
      length(unique(nmds_scores$plot_number)),
      nmds_strict$stress
    ),
    x = "NMDS1", y = "NMDS2",
    colour = NULL, shape = NULL,
    caption = "Arrows: envfit vectors (p ≤ 0.05). Lines connect first and last survey of same plot."
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")

# ggsave("fig_nmds_paired.pdf", p_nmds, width = 8, height = 7)
ggsave("fig_nmds_paired.png", p_nmds, width = 8, height = 7, dpi = 300)
message("Saved: fig_nmds_paired.png")


## ── SECTION 2 ─────────────────────────────────────────────────────────────────
## Paired change-in-composition: per-plot Bray–Curtis distance
## Purpose: visualise magnitude of compositional change per plot; test whether
##          soil type or moisture class predicts the degree of change.
## ─────────────────────────────────────────────────────────────────────────────

# 2a. Compute per-plot Bray–Curtis distance between first and last survey
bc_matrix <- as.matrix(vegdist(mat_strict, method = "bray"))

plot_bc <- meta_strict %>%
  mutate(row_idx = row_number()) %>%
  group_by(plot_number) %>%
  summarise(
    idx_first = row_idx[sample == "first"],
    idx_last  = row_idx[sample == "last"],
    .groups   = "drop"
  ) %>%
  mutate(bc_dist = map2_dbl(idx_first, idx_last, ~ bc_matrix[.x, .y]))

# 2b. Attach environmental metadata (from first survey for context)
env_first <- meta_strict %>%
  filter(sample == "first") %>%
  select(
    plot_number,
    any_of("soil_type_name"),
    moisture_name,
    háplöntuþekja,
    soil_depth,
    any_of("vegetation_height_mean"),
    total_cover
  )

plot_bc_env <- plot_bc %>%
  left_join(env_first, by = "plot_number")

# 2c. Summary: median Bray–Curtis by soil type (if available)
if ("soil_type_name" %in% names(plot_bc_env)) {
  plot_bc_env %>%
    group_by(soil_type_name) %>%
    summarise(
      n           = n(),
      median_bc   = median(bc_dist, na.rm = TRUE),
      mean_bc     = mean(bc_dist, na.rm = TRUE),
      .groups     = "drop"
    ) %>%
    arrange(desc(median_bc)) %>%
    print()
} else {
  message("Soil type not available in plot_bc_env; skipping soil-type summary.")
}

# 2d. Kruskal–Wallis: does soil type / moisture predict magnitude of change?
if ("soil_type_name" %in% names(plot_bc_env)) {
  kw_soil <- kruskal.test(bc_dist ~ soil_type_name, data = plot_bc_env)
  message(sprintf(
    "Kruskal–Wallis soil type:  χ²=%.2f, df=%d, p=%.4f",
    kw_soil$statistic, kw_soil$parameter, kw_soil$p.value
  ))
} else {
  kw_soil <- NULL
}

kw_moist <- kruskal.test(bc_dist ~ moisture_name, data = plot_bc_env)
message(sprintf(
  "Kruskal–Wallis moisture:   χ²=%.2f, df=%d, p=%.4f",
  kw_moist$statistic, kw_moist$parameter, kw_moist$p.value
))

# NOTE: if soil-type p < 0.05, follow up with Dunn test (dunn.test package)
# if (!is.null(kw_soil) && kw_soil$p.value < 0.05) {
#   dunn.test::dunn.test(plot_bc_env$bc_dist, plot_bc_env$soil_type_name, method = "bh")
# }

# 2e. Boxplot: per-plot BC distance by soil type (if available)
if ("soil_type_name" %in% names(plot_bc_env)) {
  p_bc_soil <- ggplot(
    plot_bc_env %>% filter(!is.na(soil_type_name)),
    aes(x = reorder(soil_type_name, bc_dist, FUN = median),
        y = bc_dist, fill = soil_type_name)
  ) +
    geom_boxplot(alpha = 0.7, outlier.size = 1.2) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
    labs(
      title = "Compositional turnover by soil type",
      subtitle = "Per-plot Bray–Curtis distance (first → last survey)",
      x = "Soil type", y = "Bray–Curtis distance",
      fill = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x  = element_text(angle = 35, hjust = 1),
      legend.position = "none"
    )

  # ggsave("fig_bc_by_soil.pdf", p_bc_soil, width = 8, height = 5)
  ggsave("fig_bc_by_soil.png", p_bc_soil, width = 8, height = 5, dpi = 300)
  message("Saved: fig_bc_by_soil.png")
} else {
  message("Soil type not available; skipping fig_bc_by_soil.png.")
}


## ── SECTION 3 ─────────────────────────────────────────────────────────────────
## Indicator species analysis (IndVal) – first vs last survey
## Purpose: identify taxa that significantly increased or decreased across the
##          paired resurvey, using a permutation-based approach.
## ─────────────────────────────────────────────────────────────────────────────

# indicspecies::multipatt() needs a grouping vector aligned with mat_strict rows
group_vec <- ifelse(meta_strict$sample == "first", 1L, 2L)

indval_result <- multipatt(
  mat_strict,
  cluster      = group_vec,
  control      = how(nperm = 999),
  func         = "IndVal.g",
  duleg        = TRUE
)

# Tidy the output
indval_tbl <- indval_result$sign %>%
  rownames_to_column("taxon") %>%
  as_tibble() %>%
  rename(group = index) %>%
  mutate(
    survey = if_else(group == 1L, "first", "last"),
    p_adj  = p.adjust(p.value, method = "BH")  # Benjamini–Hochberg FDR
  ) %>%
  filter(p_adj <= 0.05) %>%
  arrange(p_adj)

# Top 20 indicators for each survey
cat("\n── Top indicators for FIRST survey ──\n")
indval_tbl %>% filter(survey == "first") %>% slice_head(n = 20) %>% print()

cat("\n── Top indicators for LAST survey ──\n")
indval_tbl %>% filter(survey == "last") %>% slice_head(n = 20) %>% print()

# Save full table
write_csv(indval_tbl, "indval_first_vs_last.csv")
message("Saved: indval_first_vs_last.csv")

# 3b. Quick dot-plot of top 15 indicators per survey
top_indval <- indval_tbl %>%
  group_by(survey) %>%
  slice_min(order_by = p_adj, n = 15) %>%
  ungroup()

p_indval <- ggplot(
  top_indval,
  aes(x = stat, y = reorder(taxon, stat), fill = survey)
) +
  geom_col(alpha = 0.85) +
  facet_wrap(~survey, scales = "free_y") +
  scale_fill_manual(values = c(first = "#2980b9", last = "#e67e22")) +
  labs(
    title    = "Indicator taxa — first vs last survey",
    subtitle = "IndVal.g statistic; top 15 per group (FDR ≤ 5%)",
    x        = "IndVal statistic",
    y        = NULL,
    fill     = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(legend.position = "none")

# ggsave("fig_indval.pdf", p_indval, width = 10, height = 7)
ggsave("fig_indval.png", p_indval, width = 10, height = 7, dpi = 300)
message("Saved: fig_indval.png")


## ── SECTION 4 ─────────────────────────────────────────────────────────────────
## Univariate follow-up: change in key continuous community descriptors
## Purpose: test whether háplöntuþekja (vascular plant cover), vegetation
##          height, total cover, and soil depth changed between surveys.
## Using paired Wilcoxon signed-rank test (non-parametric, matched by plot).
## ─────────────────────────────────────────────────────────────────────────────

univar_vars <- c("háplöntuþekja", "vegetation_height_1", "total_cover", "soil_depth")

# Pivot to wide so each plot has one first-row and one last-row
env_paired_wide <- meta_strict %>%
  select(plot_number, sample, all_of(univar_vars)) %>%
  pivot_wider(
    names_from  = sample,
    values_from = all_of(univar_vars),
    names_sep   = "_"
  )

# Run paired Wilcoxon for each variable
wilcox_results <- map_dfr(univar_vars, function(v) {
  first_col <- paste0(v, "_first")
  last_col  <- paste0(v, "_last")
  x <- env_paired_wide[[first_col]]
  y <- env_paired_wide[[last_col]]
  complete_idx <- complete.cases(x, y)
  n_pairs <- sum(complete_idx)
  test <- wilcox.test(y[complete_idx], x[complete_idx],
                      paired = TRUE, conf.int = TRUE)
  tibble(
    variable    = v,
    n_pairs     = n_pairs,
    median_first = median(x[complete_idx], na.rm = TRUE),
    median_last  = median(y[complete_idx], na.rm = TRUE),
    median_change = median(y[complete_idx] - x[complete_idx], na.rm = TRUE),
    ci_low       = test$conf.int[1],
    ci_high      = test$conf.int[2],
    W            = test$statistic,
    p_value      = test$p.value
  )
}) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

cat("\n── Paired Wilcoxon: change between first and last survey ──\n")
print(wilcox_results)

write_csv(wilcox_results, "univar_paired_wilcox.csv")
message("Saved: univar_paired_wilcox.csv")

# 4b. Strip/dot-plot of pairwise change distributions
env_paired_long_change <- map_dfr(univar_vars, function(v) {
  first_col <- paste0(v, "_first")
  last_col  <- paste0(v, "_last")
  env_paired_wide %>%
    transmute(
      plot_number = plot_number,
      variable    = v,
      change      = .data[[last_col]] - .data[[first_col]]
    )
}) %>%
  filter(!is.na(change)) %>%
  left_join(
    wilcox_results %>% select(variable, p_adj),
    by = "variable"
  ) %>%
  mutate(
    sig_label = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    )
  )

p_change <- ggplot(
  env_paired_long_change,
  aes(x = variable, y = change)
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_violin(fill = "#bdc3c7", alpha = 0.5, trim = TRUE) +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 0.8) +
  geom_text(
    data = env_paired_long_change %>%
      group_by(variable, sig_label) %>% slice(1),
    aes(y = max(env_paired_long_change$change, na.rm = TRUE) * 1.05,
        label = sig_label),
    size = 5, colour = "#c0392b"
  ) +
  labs(
    title    = "Change in continuous community descriptors (last − first)",
    subtitle = "Paired Wilcoxon; * p<0.05, ** p<0.01, *** p<0.001 (BH-adjusted)",
    x        = NULL,
    y        = "Change (last − first)"
  ) +
  theme_classic(base_size = 11)

# ggsave("fig_univar_change.pdf", p_change, width = 8, height = 5)
ggsave("fig_univar_change.png", p_change, width = 8, height = 5, dpi = 300)
message("Saved: fig_univar_change.png")


## ── SECTION 5 ─────────────────────────────────────────────────────────────────
## Next decision point — what to do with these results
## ─────────────────────────────────────────────────────────────────────────────

# After running Sections 1–4, you will have:
#
#  fig_nmds_paired.png       — Ordination showing temporal shift per plot
#  fig_bc_by_soil.png        — Turnover magnitude by soil type
#  indval_first_vs_last.csv  — Indicator taxa table
#  fig_indval.png            — Top indicator taxa plot
#  univar_paired_wilcox.csv  — Univariate change results
#  fig_univar_change.png     — Change distributions for key covariates
#
# Suggested next steps depending on what you find:
#
#  A) If háplöntuþekja (vascular cover) increased significantly → consider a
#     directional dbRDA to formally partition variance explained by sample
#     (time) vs. environment:
#       vegan::dbrda(mat_strict ~ sample + háplöntuþekja + Condition(plot_number), ...)
#
#  B) If soil type or moisture has a large BC effect (Kruskal p < 0.05) →
#     split the dataset by these classes and run stratum-specific PERMANOVAs
#     to ask: "is the temporal signal consistent across soil types?"
#
#  C) If certain taxa appear as consistent gainers/losers in IndVal →
#     check whether these are habitat specialists with known climate or
#     grazing sensitivities (Icelandic vegetation context: Racomitrium heath,
#     Dryas heaths, wetland sedges, etc.)
#
#  D) If temporal shift is spatially structured → explore a Mantel test or
#     spatial autocorrelation of bc_dist vs. geographic distance using
#     longitude_1/latitude_1 from v_plants_use.

message("\nAll done. Review outputs and decide on follow-up based on results.")
sessionInfo()


## =============================================================================
## v_plants — NMDS revision after outlier diagnosis
## Excludes RN-02, RN-03, RN-65-01 from ordination (beach/dune + marsh plots
## that sit ~17 NMDS units from the main cloud and distort axis scaling).
## All three are retained in PERMANOVA results from the screening script.
##
## Assumes in environment:
##   community_wide_strict  – wide matrix with plot_number + sample cols
##   meta_strict            – metadata aligned row-for-row with community_wide_strict
##   mat_strict             – numeric species matrix (318 rows × 747 cols)
##   vars_driver_strict_final, driver_final – from variable screening
## =============================================================================

library(tidyverse)
library(vegan)
library(ggrepel)   # install.packages("ggrepel") if missing

set.seed(123)

## ── 0. Define and remove outlier plots ────────────────────────────────────────
outlier_plots <- c("RN-02", "RN-03", "RN-65-01")

# Boolean index for rows to KEEP
keep_rows <- !meta_strict$plot_number %in% outlier_plots

meta_core   <- meta_strict[keep_rows, ]
mat_core    <- mat_strict[keep_rows, ]

# Drop any taxa that are now all-zero (only present in removed plots)
taxa_present <- colSums(mat_core) > 0
mat_core     <- mat_core[, taxa_present]

message(sprintf(
  "Removed %d plots (%d rows). Matrix now: %d rows × %d taxa.",
  length(outlier_plots), sum(!keep_rows),
  nrow(mat_core), ncol(mat_core)
))


## ── 1. Re-run NMDS on core dataset ────────────────────────────────────────────
# trymax = 200 for a cleaner solution; trace = FALSE keeps console quiet
nmds_core <- metaMDS(
  mat_core,
  distance = "bray",
  k        = 2,
  trymax   = 200,
  trace    = FALSE
)

message(sprintf("NMDS stress (core, k=2): %.4f", nmds_core$stress))
# Stress guide:
#   < 0.10  excellent   |  0.10–0.15  good
#   0.15–0.20  acceptable  |  > 0.20  use with caution (consider k=3)

if (nmds_core$stress > 0.20) {
  warning("Stress > 0.20. Consider running with k = 3.")
}


## ── 2. Re-fit environmental vectors on core dataset ───────────────────────────
# Align env data to the core rows
env_core <- meta_core %>%
  select(all_of(vars_driver_strict_final)) %>%
  mutate(
    across(where(is.character), as.factor),
    across(where(is.logical),   as.factor)
  )

fit_core <- envfit(nmds_core, env_core, permutations = 999)

# Build tidy arrow table (numeric vectors only for biplot arrows)
sig_thresh <- 0.05

arrow_df <- as.data.frame(fit_core$vectors$arrows) %>%
  rownames_to_column("variable") %>%
  rename(ax1 = NMDS1, ax2 = NMDS2) %>%
  mutate(
    r2      = fit_core$vectors$r,
    p_value = fit_core$vectors$pvals
  ) %>%
  filter(p_value <= sig_thresh)

# Factor centroids (plotted as text labels, not arrows)
factor_df <- as.data.frame(scores(fit_core, display = "factors")) %>%
  rownames_to_column("label") %>%
  rename(fax1 = NMDS1, fax2 = NMDS2) %>%
  mutate(
    p_value = fit_core$factors$pvals[
      # match factor name prefix to pvals names
      sapply(rownames(as.data.frame(scores(fit_core, display = "factors"))),
             function(nm) {
               matched <- names(fit_core$factors$pvals)[
                 startsWith(nm, names(fit_core$factors$pvals))
               ]
               if (length(matched)) matched[1] else NA_character_
             })
    ]
  ) %>%
  filter(!is.na(p_value), p_value <= sig_thresh)


## ── 3. Extract NMDS scores ────────────────────────────────────────────────────
site_scores <- as.data.frame(scores(nmds_core, display = "sites")) %>%
  bind_cols(meta_core %>% select(plot_number, sample))

# Scale arrows to ~70 % of ordination radius
ordi_radius <- max(abs(c(site_scores$NMDS1, site_scores$NMDS2))) * 0.70

arrow_df <- arrow_df %>%
  mutate(
    xend = ax1 * ordi_radius,
    yend = ax2 * ordi_radius
  )


## ── 4. Build clean biplot ─────────────────────────────────────────────────────
p_nmds_clean <- ggplot(site_scores, aes(x = NMDS1, y = NMDS2)) +

  # ── paired change lines (lightest layer, drawn first)
  geom_line(
    aes(group = plot_number),
    colour   = "grey60",
    linewidth = 0.30,
    alpha    = 0.55
  ) +

  # ── site points
  geom_point(
    aes(colour = sample, shape = sample),
    size  = 2.0,
    alpha = 0.85
  ) +

  # ── environmental arrows
  geom_segment(
    data = arrow_df,
    aes(x = 0, y = 0, xend = xend, yend = yend),
    arrow     = arrow(length = unit(0.22, "cm"), type = "closed"),
    colour    = "#c0392b",
    linewidth  = 0.75,
    inherit.aes = FALSE
  ) +

  # ── arrow labels (repelled so they don't pile up)
  geom_label_repel(
    data        = arrow_df,
    aes(x = xend * 1.12, y = yend * 1.12, label = variable),
    size        = 3.0,
    colour      = "#c0392b",
    fill        = alpha("white", 0.80),
    label.size  = 0,
    max.overlaps = 20,
    seed        = 42,
    inherit.aes = FALSE
  ) +

  # ── factor centroids (if any passed significance)
  {
    if (nrow(factor_df) > 0)
      geom_label_repel(
        data        = factor_df,
        aes(x = fax1, y = fax2, label = label),
        size        = 2.5,
        colour      = "#27ae60",
        fill        = alpha("white", 0.80),
        label.size  = 0,
        max.overlaps = 15,
        seed        = 43,
        inherit.aes = FALSE
      )
    else
      NULL
  } +

  # ── scales and labels
  scale_colour_manual(
    values = c(first = "#2980b9", last = "#e67e22"),
    labels = c(first = "First survey", last = "Last survey")
  ) +
  scale_shape_manual(
    values = c(first = 16, last = 17),
    labels = c(first = "First survey", last = "Last survey")
  ) +
  coord_equal() +   # equal axes so arrow angles are honest
  labs(
    title    = "NMDS of plant community composition (Bray–Curtis)",
    subtitle = sprintf(
      "%d plots · stress = %.3f · PERMANOVA R² = 2.7%%, p = 0.001  |  excl. RN-02, RN-03, RN-65-01",
      n_distinct(site_scores$plot_number),
      nmds_core$stress
    ),
    x       = "NMDS1",
    y       = "NMDS2",
    colour  = NULL,
    shape   = NULL,
    caption = paste0(
      "Red arrows: envfit vectors (p ≤ 0.05); length ∝ r².\n",
      "Grey lines connect first and last survey of same plot."
    )
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position  = "bottom",
    plot.subtitle    = element_text(size = 8, colour = "grey40"),
    plot.caption     = element_text(size = 7.5, colour = "grey50", hjust = 0)
  )

# ggsave("fig_nmds_core.pdf",  p_nmds_clean, width = 8,  height = 7.5)
ggsave("fig_nmds_core.png",  p_nmds_clean, width = 8,  height = 7.5, dpi = 300)
message("Saved: fig_nmds_core.png")


## ── 5. Quick stress diagnostic (Shepard plot) ─────────────────────────────────
# Useful to include as a supplementary figure or to justify k choice.
# pdf("fig_shepard_core.pdf", width = 5, height = 5)
#   stressplot(nmds_core, main = sprintf("Shepard plot — stress %.3f", nmds_core$stress))
# dev.off()
# message("Saved: fig_shepard_core.pdf")
png("fig_shepard_core.png", width = 1500, height = 1500, res = 300)
stressplot(nmds_core, main = sprintf("Shepard plot — stress %.3f", nmds_core$stress))
dev.off()
message("Saved: fig_shepard_core.png")


## ── 6. What changed? Report updated envfit r² for core dataset ────────────────
fit_tbl_core <- bind_rows(
  tibble(
    variable = names(fit_core$vectors$r),
    r2       = as.numeric(fit_core$vectors$r),
    p_value  = as.numeric(fit_core$vectors$pvals),
    kind     = "numeric"
  ),
  if (!is.null(fit_core$factors)) tibble(
    variable = names(fit_core$factors$r),
    r2       = as.numeric(fit_core$factors$r),
    p_value  = as.numeric(fit_core$factors$pvals),
    kind     = "factor"
  )
) %>%
  arrange(p_value)

cat("\n── envfit results on CORE dataset (outliers removed) ──\n")
print(fit_tbl_core)

write_csv(fit_tbl_core, "envfit_core.csv")
message("Saved: envfit_core.csv")


## ── 7. Separate coastal NMDS (the excluded plots as their own ordination) ──────
# Even though only 3 plots (6 rows), it is worth documenting what those
# communities look like — and whether RN-02/RN-03 shifted notably.

outlier_rows  <- !keep_rows   # rows corresponding to the 3 excluded plots
mat_coastal   <- mat_strict[outlier_rows, ]
meta_coastal  <- meta_strict[outlier_rows, ]

# Drop all-zero taxa
mat_coastal <- mat_coastal[, colSums(mat_coastal) > 0]

message(sprintf(
  "\nCoastal sub-dataset: %d rows × %d taxa.",
  nrow(mat_coastal), ncol(mat_coastal)
))

# With only 6 rows, NMDS is unstable — use PCoA instead
pcoa_coastal <- cmdscale(
  vegdist(mat_coastal, method = "bray"),
  k   = 2,
  eig = TRUE
)

pcoa_df <- as.data.frame(pcoa_coastal$points) %>%
  setNames(c("Axis1", "Axis2")) %>%
  bind_cols(meta_coastal %>% select(plot_number, sample))

var_explained <- round(
  100 * pcoa_coastal$eig[1:2] / sum(abs(pcoa_coastal$eig)), 1
)

p_coastal <- ggplot(pcoa_df, aes(x = Axis1, y = Axis2)) +
  geom_line(aes(group = plot_number), colour = "grey60", linewidth = 0.4) +
  geom_point(aes(colour = sample, shape = sample), size = 3) +
  geom_label_repel(
    aes(label = paste(plot_number, sample)),
    size = 3, max.overlaps = 20, seed = 1
  ) +
  scale_colour_manual(
    values = c(first = "#2980b9", last = "#e67e22"),
    labels = c(first = "First survey", last = "Last survey")
  ) +
  scale_shape_manual(
    values = c(first = 16, last = 17),
    labels = c(first = "First survey", last = "Last survey")
  ) +
  coord_equal() +
  labs(
    title    = "PCoA — excluded coastal/outlier plots",
    subtitle = sprintf(
      "Bray–Curtis · Axis 1 = %.1f%% · Axis 2 = %.1f%%",
      var_explained[1], var_explained[2]
    ),
    x      = sprintf("PCoA Axis 1 (%.1f%%)", var_explained[1]),
    y      = sprintf("PCoA Axis 2 (%.1f%%)", var_explained[2]),
    colour = NULL, shape = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")

# ggsave("fig_pcoa_coastal.pdf", p_coastal, width = 6, height = 6)
ggsave("fig_pcoa_coastal.png", p_coastal, width = 6, height = 6, dpi = 300)
message("Saved: fig_pcoa_coastal.png")

message("\nDone. Next: review fig_nmds_core.png and check stress.")
message("If stress improved below 0.15, the NMDS is solid for publication.")
message("If still >0.15, run with k=3 and use a 3D rotation or pairs plot.")
