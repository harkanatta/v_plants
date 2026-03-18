# ── 0. Build filtered matrix ──────────────────────────────────────────────────
coastal_types <- c("Sandstrandarvist", "Strandmelhólavist",
                   "Sjávarfitjungsvist", "Sjávarkletta- og eyjavist")

# Use rass (already filtered) to rebuild meta and matrix
meta_inland <- meta_strict %>%
  filter(!paste(plot_number, sample, sep = "_") %in% 
           (v_plants_use %>%
              filter(habitat_type_name %in% coastal_types) %>%
              mutate(pv = paste(plot_number, sample, sep = "_")) %>%
              pull(pv)))

mat_inland <- mat_strict[rownames(mat_strict) %in% 
                           paste(meta_inland$plot_number, 
                                 meta_inland$sample, sep = "_"), ]
mat_inland <- mat_inland[, colSums(mat_inland) > 0]

# Sanity check
cat("Plots:", nrow(mat_inland), "| Species:", ncol(mat_inland), "\n")
c("Mertensia maritima", "Honckenya peploides", "Cakile maritima") %in% colnames(mat_inland)

# ── 1. Run DCA ────────────────────────────────────────────────────────────────
dca <- decorana(mat_inland)
print(dca)

# ── 2. Site scores ────────────────────────────────────────────────────────────
site_sc <- scores(dca, display = "sites") %>%
  as_tibble(rownames = "plot_visit")

# ── 3. Species scores ─────────────────────────────────────────────────────────
sp_sc <- scores(dca, display = "species") %>%
  as_tibble(rownames = "taxon_name")

# ── 4. Join metadata ──────────────────────────────────────────────────────────
extra_env_vars <- c(
  "vegetation_height_mean", "mosaþekja", "total_cover",
  "háplöntuþekja", "fléttuþekja", "Mosa- og fléttuskán"
)

env_extra <- env_candidates_ps %>%
  filter(!paste(plot_number, sample, sep = "_") %in%
           paste(meta_strict$plot_number, meta_strict$sample, sep = "_")[
             !paste(meta_strict$plot_number, meta_strict$sample, sep = "_") %in%
               rownames(mat_inland)]) %>%
  select(plot_number, sample, all_of(extra_env_vars))

site_env <- site_sc %>%
  mutate(
    plot_number = str_remove(plot_visit, "_(?:first|last)$"),
    sample      = str_extract(plot_visit, "(?:first|last)$")
  ) %>%
  left_join(meta_inland, by = c("plot_number", "sample")) %>%
  left_join(env_extra,   by = c("plot_number", "sample"))

cat("NAs in haeddypimetrar:", sum(is.na(site_env$haeddypimetrar)), "\n")
cat("NAs in vegetation_height_mean:", sum(is.na(site_env$vegetation_height_mean)), "\n")
cat("Rows:", nrow(site_env), "— expected:", nrow(mat_inland), "\n")

# ── 5. envfit — continuous ────────────────────────────────────────────────────
env_num_mat <- site_env %>%
  select(vegetation_height_mean,
         mosaþekja, total_cover, háplöntuþekja,
         fléttuþekja, `Mosa- og fléttuskán`, grýtniþekja)

set.seed(123)
ef_num  <- envfit(dca, env_num_mat, permutations = 999, na.rm = TRUE)
print(ef_num)

# ── 6. envfit — categorical ───────────────────────────────────────────────────
ef_topo  <- envfit(dca, site_env["topography_name"], permutations = 999, na.rm = TRUE)
ef_moist <- envfit(dca, site_env["moisture_name"],   permutations = 999, na.rm = TRUE)
print(ef_topo)
print(ef_moist)

# ── 7. Arrow and centroid data ────────────────────────────────────────────────
ef_arrows <- scores(ef_num, display = "vectors") %>%
  as_tibble(rownames = "variable") %>%
  mutate(r2 = ef_num$vectors$r, p = ef_num$vectors$pvals) %>%
  filter(p < 0.05)

topo_centroids <- scores(ef_topo, display = "factors") %>%
  as_tibble(rownames = "label")

# ── 8. Plot ───────────────────────────────────────────────────────────────────
arrow_scalar <- 2.0

p_dca <- ggplot() +
  geom_point(
    data = site_env,
    aes(x = DCA1, y = DCA2, colour = moisture_name, shape = sample),
    size = 2.5, alpha = 0.8
  ) +
  geom_segment(
    data = ef_arrows,
    aes(x = 0, y = 0,
        xend = DCA1 * arrow_scalar,
        yend = DCA2 * arrow_scalar),
    arrow = arrow(length = unit(0.22, "cm"), type = "closed"),
    colour = "grey30", linewidth = 0.55
  ) +
  geom_text(
    data = ef_arrows,
    aes(x = DCA1 * arrow_scalar * 1.13,
        y = DCA2 * arrow_scalar * 1.13,
        label = variable),
    size = 2.9, colour = "grey20"
  ) +
  geom_point(
    data = topo_centroids,
    aes(x = DCA1, y = DCA2),
    shape = 17, size = 3.5, colour = "black"
  ) +
  ggrepel::geom_text_repel(
    data = topo_centroids,
    aes(x = DCA1, y = DCA2, label = label),
    size = 3, colour = "black", fontface = "bold",
    box.padding = 0.4, max.overlaps = 20
  ) +
  scale_colour_brewer(palette = "Set2", name = "Moisture") +
  scale_shape_manual(values = c(first = 16, last = 17), name = "Survey") +
  labs(
    x = paste0("DCA1 (length = ", round(dca$evals.decorana[1], 2), " SD)"),
    y = paste0("DCA2 (length = ", round(dca$evals.decorana[2], 2), " SD)"),
    title    = "DCA — Icelandic vascular plant communities (coastal plots removed)",
    subtitle = "▲ Topography centroids  |  Arrows: sig. env. drivers (p < 0.05, envfit)"
  ) +
  theme_classic() +
  theme(legend.position = "right")

ggsave("outputs/p_dca_inland.png",
       plot   = p_dca,
       width  = 16, height = 17, units = "cm",
       dpi    = 300, bg = "white")

# ── 9. Combined plot: paired arrows + centroids (points behind) ─────────────────
site_pairs <- site_env %>%
  select(plot_number, sample, DCA1, DCA2, moisture_name) %>%
  pivot_wider(names_from = sample,
              values_from = c(DCA1, DCA2),
              names_sep = "_") %>%
  filter(!is.na(DCA1_first) & !is.na(DCA1_last))

p_dca_combined <- ggplot() +
  # Points and segments first (behind)
  geom_segment(
    data = site_pairs,
    aes(x = DCA1_first, y = DCA2_first,
        xend = DCA1_last, yend = DCA2_last,
        colour = moisture_name),
    alpha = 0.4, linewidth = 0.4
  ) +
  geom_point(
    data = site_env,
    aes(x = DCA1, y = DCA2, colour = moisture_name, shape = sample),
    size = 2.5, alpha = 0.8
  ) +
  # Env arrows and labels on top
  geom_segment(
    data = ef_arrows,
    aes(x = 0, y = 0,
        xend = DCA1 * arrow_scalar,
        yend = DCA2 * arrow_scalar),
    arrow = arrow(length = unit(0.22, "cm"), type = "closed"),
    colour = "grey30", linewidth = 0.55
  ) +
  geom_text(
    data = ef_arrows,
    aes(x = DCA1 * arrow_scalar * 1.13,
        y = DCA2 * arrow_scalar * 1.13,
        label = variable),
    size = 2.9, colour = "grey20"
  ) +
  geom_point(
    data = topo_centroids,
    aes(x = DCA1, y = DCA2),
    shape = 17, size = 3.5, colour = "black"
  ) +
  ggrepel::geom_text_repel(
    data = topo_centroids,
    aes(x = DCA1, y = DCA2, label = label),
    size = 3, colour = "black", fontface = "bold",
    box.padding = 0.4, max.overlaps = 20
  ) +
  scale_colour_brewer(palette = "Set2", name = "Moisture") +
  scale_shape_manual(values = c(first = 16, last = 17), name = "Survey") +
  labs(
    x = paste0("DCA1 (length = ", round(dca$evals.decorana[1], 2), " SD)"),
    y = paste0("DCA2 (length = ", round(dca$evals.decorana[2], 2), " SD)"),
    title    = "DCA — Icelandic vascular plant communities (coastal plots removed)",
    subtitle = "▲ Topography centroids  |  Arrows: sig. env. drivers (p < 0.05)  |  Lines: first (●) → last (▲) per plot"
  ) +
  theme_classic() +
  theme(legend.position = "right")

ggsave("outputs/p_dca_combined.png",
       plot = p_dca_combined,
       width = 16, height = 17, units = "cm",
       dpi = 300, bg = "white")
















# ── Species richness per plot-visit ───────────────────────────────────────────

# Add first/last via join (paired_long has plot_number, year, sample)
sample_lookup <- paired_long %>% distinct(plot_number, year, sample)

# Mean subplot richness per plot-visit (alpha diversity flavour)
richness_env <- v_plants_use %>%
  inner_join(sample_lookup, by = c("plot_number", "year")) %>%
  group_by(plot_number, subplot_number, sample) %>%
  summarise(n_taxa_subplot = n_distinct(taxon_name), .groups = "drop") %>%
  group_by(plot_number, sample) %>%
  summarise(
    mean_subplot_richness = mean(n_taxa_subplot),
    max_subplot_richness  = max(n_taxa_subplot),
    n_subplots            = n(),
    .groups = "drop"
  )

# Total richness per plot-visit (gamma at plot level)
total_richness <- v_plants_use %>%
  inner_join(sample_lookup, by = c("plot_number", "year")) %>%
  group_by(plot_number, sample) %>%
  summarise(total_plot_richness = n_distinct(taxon_name), .groups = "drop")

richness_env <- richness_env %>%
  left_join(total_richness, by = c("plot_number", "sample"))

richness_env


site_env <- site_env %>%
  left_join(richness_env, by = c("plot_number", "sample"))

# Sanity check
cat("NAs in mean_subplot_richness:", sum(is.na(site_env$mean_subplot_richness)), "\n")
cat("NAs in total_plot_richness:",   sum(is.na(site_env$total_plot_richness)), "\n")

# Updated envfit matrix
env_num_mat2 <- site_env %>%
  select(vegetation_height_mean,
         mosaþekja, total_cover, háplöntuþekja,
         fléttuþekja, `Mosa- og fléttuskán`, grýtniþekja,
         mean_subplot_richness, total_plot_richness)

set.seed(123)
ef_num2 <- envfit(dca, env_num_mat2, permutations = 999, na.rm = TRUE)
print(ef_num2)








