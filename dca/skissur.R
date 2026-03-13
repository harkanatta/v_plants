dca_vars <- v_plants_use %>%
  select(
    taxon_name,
    plot_number,
    subplot_number,
    haeddypimetrar,
    year,
    vegetation_height_mean,
    mosaþekja,
    total_cover,
    háplöntuþekja,
    fléttuþekja,
    `Mosa- og fléttuskán`,
    grýtniþekja,
    topography_name,
    field_date_parsed,
    moisture_name,
    soil_depth,
    soil_type_name,
    habitat_type_name,
    permafrost
  )

dca_vars


n_taxa <- v_plants_use %>%
  group_by(plot_number, subplot_number) %>%
  summarise(n_taxa = n_distinct(taxon_name), .groups = "drop")





# ── 0. Stamp rownames (meta_strict is row-aligned with mat_strict by construction)
rownames(mat_strict) <- paste(meta_strict$plot_number, meta_strict$sample, sep = "_")

# ── 1. Run DCA ────────────────────────────────────────────────────────────────
dca <- decorana(mat_strict)
print(dca)

# ── 2. Site scores ────────────────────────────────────────────────────────────
site_sc <- scores(dca, display = "sites") %>%
  as_tibble(rownames = "plot_visit")

# ── 3. Species scores ─────────────────────────────────────────────────────────
sp_sc <- scores(dca, display = "species") %>%
  as_tibble(rownames = "taxon_name")

# ── 4. Join metadata ──────────────────────────────────────────────────────────
# meta_strict has: plot_number, sample, vars_driver_strict_final, 
#                  habitat_type_name, grýtniþekja, gap_years, haeddypimetrar
# env_candidates_ps has the remaining display vars we want
extra_env_vars <- c(
  "vegetation_height_mean", "mosaþekja", "total_cover",
  "háplöntuþekja", "fléttuþekja", "Mosa- og fléttuskán"
)

env_extra <- env_candidates_ps %>%
  select(plot_number, sample, all_of(extra_env_vars))

site_env <- site_sc %>%
  mutate(
    plot_number = str_remove(plot_visit, "_(?:first|last)$"),
    sample      = str_extract(plot_visit, "(?:first|last)$")
  ) %>%
  left_join(meta_strict, by = c("plot_number", "sample")) %>%
  left_join(env_extra,   by = c("plot_number", "sample"))

# Sanity check
cat("NAs in haeddypimetrar:", sum(is.na(site_env$haeddypimetrar)), "\n")
cat("NAs in vegetation_height_mean:", sum(is.na(site_env$vegetation_height_mean)), "\n")
cat("Rows:", nrow(site_env), "— expected:", nrow(mat_strict), "\n")

# ── 5. envfit — continuous variables ──────────────────────────────────────────
env_num_mat <- site_env %>%
  select(
    haeddypimetrar, vegetation_height_mean,
    mosaþekja, total_cover, háplöntuþekja,
    fléttuþekja, `Mosa- og fléttuskán`, grýtniþekja
  )

set.seed(123)
ef_num  <- envfit(dca, env_num_mat, permutations = 999, na.rm = TRUE)
print(ef_num)

# ── 6. envfit — categorical ───────────────────────────────────────────────────
ef_topo <- envfit(dca, site_env["topography_name"], permutations = 999, na.rm = TRUE)
ef_moist <- envfit(dca, site_env["moisture_name"],  permutations = 999, na.rm = TRUE)
print(ef_topo)
print(ef_moist)

# ── 7. Arrow and centroid data for ggplot ─────────────────────────────────────
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
    title  = "DCA — Icelandic vascular plant communities",
    subtitle = "▲ Topography centroids  |  Arrows: sig. env. drivers (p < 0.05, envfit)"
  ) +
  theme_classic() +
  theme(legend.position = "right")


ggsave("outputs/p_dca.png",
       plot = p_dca,
       width  = 16,   # cm — square-ish canvas to match coord_equal() panel
       height = 17,   # cm — extra height for title + legend + caption
       units  = "cm",
       dpi    = 300,
       bg     = "white")
