# Within-plot change: first vs last survey (two-wave panel)
# Requires: v_plants2, plot_years. Uses vegan for diversity and distance.
library(tidyverse)
library(vegan)

# ------------------------------------------------------------------------------
# 1. Restrict to revisited plots and parse first/last year
# ------------------------------------------------------------------------------
revisited <- plot_years %>%
  filter(n_years == 2) %>%
  separate(years, into = c("year_first", "year_last"), sep = ", ", convert = TRUE)

# ------------------------------------------------------------------------------
# 2. Cover per plot × year × taxon (pool subplots: mean cover)
# ------------------------------------------------------------------------------
dat <- v_plants2 %>%
  filter(plot_number %in% revisited$plot_number) %>%
  mutate(cover_num = coalesce(scale_value_value,
    case_when(scale_value_code %in% c("•", "·", "\u2022") ~ 0.25,
              scale_value_code == "+" ~ 0.75,
              scale_value_code == "1" ~ 3, scale_value_code == "2" ~ 15,
              scale_value_code == "3" ~ 37.5, scale_value_code == "4" ~ 62.5,
              scale_value_code == "5" ~ 87.5, TRUE ~ NA_real_)))

pooled_plot_year <- dat %>%
  filter(!is.na(cover_num), cover_num > 0) %>%
  group_by(plot_number, year, taxon_name) %>%
  summarise(cover = mean(cover_num, na.rm = TRUE), .groups = "drop")

# ------------------------------------------------------------------------------
# 3. Diversity per plot × year (S, Shannon, Simpson, Evenness)
# ------------------------------------------------------------------------------
wide_plot_year <- pooled_plot_year %>%
  pivot_wider(names_from = taxon_name, values_from = cover, values_fill = 0)

split_by_plot_year <- wide_plot_year %>%
  group_by(plot_number, year) %>%
  group_split()

species_matrices <- split_by_plot_year %>%
  map(~ .x %>% select(-plot_number, -year) %>% as.matrix() %>% `[`(1, , drop = FALSE))
ids <- split_by_plot_year %>% map(~ .x %>% select(plot_number, year) %>% slice(1))

diversity_plot_year <- tibble(
  plot_number = map_chr(ids, ~ .x$plot_number),
  year        = map_dbl(ids, ~ .x$year),
  N           = map_dbl(species_matrices, ~ sum(.x)),
  S           = map_dbl(species_matrices, ~ specnumber(.x)),
  Shannon     = map_dbl(species_matrices, ~ diversity(.x, index = "shannon")),
  Simpson     = map_dbl(species_matrices, ~ diversity(.x, index = "simpson")),
  Evenness    = NA_real_
) %>%
  mutate(Evenness = Shannon / log(pmax(S, 1)))

# ------------------------------------------------------------------------------
# 4. Paired change: first vs last (per plot)
# ------------------------------------------------------------------------------
paired <- diversity_plot_year %>%
  inner_join(revisited %>% select(plot_number, year_first, year_last), by = "plot_number") %>%
  filter(year %in% c(year_first, year_last)) %>%
  select(plot_number, year, S, Shannon, Simpson, Evenness, year_first, year_last)

first_survey  <- paired %>% filter(year == year_first) %>% select(plot_number, S_first = S, Shannon_first = Shannon, Simpson_first = Simpson, Evenness_first = Evenness)
last_survey   <- paired %>% filter(year == year_last)  %>% select(plot_number, S_last = S, Shannon_last = Shannon, Simpson_last = Simpson, Evenness_last = Evenness)

change <- first_survey %>%
  full_join(last_survey, by = "plot_number") %>%
  mutate(
    delta_S       = S_last - S_first,
    delta_Shannon = Shannon_last - Shannon_first,
    delta_Simpson = Simpson_last - Simpson_first,
    delta_Evenness = Evenness_last - Evenness_first
  )

# Mean and median within-plot change
change_summary <- change %>%
  summarise(
    n_plots = n(),
    mean_delta_S       = mean(delta_S, na.rm = TRUE),
    median_delta_S     = median(delta_S, na.rm = TRUE),
    mean_delta_Shannon = mean(delta_Shannon, na.rm = TRUE),
    median_delta_Shannon = median(delta_Shannon, na.rm = TRUE),
    mean_delta_Simpson = mean(delta_Simpson, na.rm = TRUE),
    median_delta_Simpson = median(delta_Simpson, na.rm = TRUE),
    mean_delta_Evenness = mean(delta_Evenness, na.rm = TRUE),
    median_delta_Evenness = median(delta_Evenness, na.rm = TRUE)
  )
print(change_summary)

# ------------------------------------------------------------------------------
# 5. Change in distribution tails (e.g. loss of high-cover taxa)
# ------------------------------------------------------------------------------
# Per plot × year: proportion of total cover in "high cover" taxa (e.g. single taxon > 50% or top decile)
pooled_with_total <- pooled_plot_year %>%
  group_by(plot_number, year) %>%
  mutate(
    total_cover = sum(cover),
    pct_cover  = cover / total_cover,
    rank_cover = row_number(desc(cover))
  ) %>%
  ungroup()

# Option A: share of cover in dominant taxon (max pct per plot-year)
dominance <- pooled_with_total %>%
  group_by(plot_number, year) %>%
  summarise(
    max_pct = max(pct_cover, na.rm = TRUE),
    n_high  = sum(pct_cover >= 0.25, na.rm = TRUE),  # taxa with ≥25% cover
    .groups = "drop"
  )

revisited_years <- revisited %>% select(plot_number, year_first, year_last)
dom_first <- dominance %>%
  inner_join(revisited_years, by = "plot_number") %>%
  filter(year == year_first) %>%
  select(plot_number, max_pct_first = max_pct, n_high_first = n_high)
dom_last <- dominance %>%
  inner_join(revisited_years, by = "plot_number") %>%
  filter(year == year_last) %>%
  select(plot_number, max_pct_last = max_pct, n_high_last = n_high)

tail_change <- dom_first %>%
  full_join(dom_last, by = "plot_number") %>%
  mutate(
    delta_max_pct = max_pct_last - max_pct_first,
    delta_n_high  = n_high_last - n_high_first
  )

tail_summary <- tail_change %>%
  summarise(
    mean_delta_max_pct = mean(delta_max_pct, na.rm = TRUE),
    median_delta_max_pct = median(delta_max_pct, na.rm = TRUE),
    mean_delta_n_high  = mean(delta_n_high, na.rm = TRUE),
    median_delta_n_high  = median(delta_n_high, na.rm = TRUE)
  )
print(tail_summary)

# ------------------------------------------------------------------------------
# 5b. Which taxa drive the increase in dominance? (dominant taxon per plot × year)
# ------------------------------------------------------------------------------
# Dominant taxon = taxon with max cover share in each plot × year (rank_cover == 1)
dominant_taxon_per_plot_year <- pooled_with_total %>%
  filter(rank_cover == 1L) %>%
  select(plot_number, year, dominant_taxon = taxon_name, max_pct = pct_cover)

dom_taxon_first <- dominant_taxon_per_plot_year %>%
  inner_join(revisited_years, by = "plot_number") %>%
  filter(year == year_first) %>%
  select(plot_number, dominant_first = dominant_taxon, max_pct_first = max_pct)
dom_taxon_last <- dominant_taxon_per_plot_year %>%
  inner_join(revisited_years, by = "plot_number") %>%
  filter(year == year_last) %>%
  select(plot_number, dominant_last = dominant_taxon, max_pct_last = max_pct)

# Per-plot: who was dominant at first vs last (and change in max share)
dominant_taxon_change <- dom_taxon_first %>%
  full_join(dom_taxon_last, by = "plot_number") %>%
  mutate(
    same_dominant = (dominant_first == dominant_last),
    delta_max_pct = max_pct_last - max_pct_first
  )

# Summary: which taxa are dominant at first vs last (frequency)
dominant_first_freq <- dominant_taxon_change %>%
  count(dominant_first, name = "n_plots_first") %>%
  filter(!is.na(dominant_first)) %>%
  arrange(desc(n_plots_first))
dominant_last_freq <- dominant_taxon_change %>%
  count(dominant_last, name = "n_plots_last") %>%
  filter(!is.na(dominant_last)) %>%
  arrange(desc(n_plots_last))

# Culprits: taxa that are dominant at last survey (especially where dominance increased)
dominant_last_freq
# Where dominance increased (delta_max_pct > 0): who is the new/larger dominant?
culprits_increased_dominance <- dominant_taxon_change %>%
  filter(delta_max_pct > 0) %>%
  count(dominant_last, name = "n_plots_increased") %>%
  arrange(desc(n_plots_increased))
# Same but with mean delta_max_pct per dominant_at_last
culprits_with_delta <- dominant_taxon_change %>%
  filter(delta_max_pct > 0) %>%
  group_by(dominant_last) %>%
  summarise(
    n_plots = n(),
    mean_delta_max_pct = mean(delta_max_pct, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_plots))

# Optional: taxa that *gained* dominance (not dominant at first, dominant at last)
gained_dominance <- dominant_taxon_change %>%
  filter(dominant_first != dominant_last | is.na(dominant_first)) %>%
  count(dominant_last, name = "n_plots_gained") %>%
  filter(!is.na(dominant_last)) %>%
  arrange(desc(n_plots_gained))

# Print key tables (pair with dominance plot)
message("Dominant taxon at FIRST survey (top 15):")
print(dominant_first_freq %>% head(15))
message("Dominant taxon at LAST survey (top 15):")
print(dominant_last_freq %>% head(15))
message("Taxa driving INCREASE in dominance (dominant at last when delta_max_pct > 0):")
print(culprits_with_delta)
message("Taxa that GAINED dominance (different taxon dominant at last):")
print(gained_dominance %>% head(15))

# ------------------------------------------------------------------------------
# 5c. Prepare slope data (kept for downstream plot objects)
# ------------------------------------------------------------------------------
# Build flow table: same taxon on both sides, ranked by n_plots in each period
n_slope <- 10L   # number of taxa to show (rows)
first_wide <- dominant_first_freq %>% rename(taxon = dominant_first, n_first = n_plots_first)
last_wide  <- dominant_last_freq  %>% rename(taxon = dominant_last,  n_last = n_plots_last)
slope_taxa <- full_join(first_wide, last_wide, by = "taxon") %>%
  mutate(
    n_first = replace_na(n_first, 0L),
    n_last  = replace_na(n_last,  0L)
  ) %>%
  filter(n_first > 0 | n_last > 0) %>%
  mutate(
    rank_first = min_rank(desc(n_first)),
    rank_last  = min_rank(desc(n_last))
  ) %>%
  filter(rank_first <= n_slope | rank_last <= n_slope) %>%
  slice_max(n_first + n_last, n = n_slope, with_ties = FALSE) %>%
  mutate(
    rank_first = min_rank(desc(n_first)),
    rank_last  = min_rank(desc(n_last))
  ) %>%
  arrange(rank_first)

# ------------------------------------------------------------------------------
# 5d. Slope graph in ggplot2 (same idea, no extra package)
#     Optional: CGPfunctions::newggslopegraph() or ggalluvial for ribbon-style.
# ------------------------------------------------------------------------------
# Long format: one row per taxon per period; include rank so y-axis = rank (no overlap)
slope_taxa_long <- slope_taxa %>%
  select(taxon, n_first, n_last, rank_first, rank_last) %>%
  pivot_longer(
    cols = c(n_first, n_last),
    names_to = "period_n",
    values_to = "n_plots"
  ) %>%
  mutate(
    period = if_else(period_n == "n_first", "First survey", "Last survey"),
    rank   = if_else(period_n == "n_first", rank_first, rank_last)
  ) %>%
  select(taxon, period, n_plots, rank)

# Slope graph: y = rank (1 at top) so each side has clear order; colour = taxon (line + point match)
p_slope_gg <- ggplot(slope_taxa_long, aes(x = period, y = rank, group = taxon, colour = taxon)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.5) +
  geom_text(
    data = subset(slope_taxa_long, period == "First survey"),
    aes(label = str_trunc(taxon, 16)), hjust = 1.08, size = 3, show.legend = FALSE
  ) +
  geom_text(
    data = subset(slope_taxa_long, period == "Last survey"),
    aes(label = str_trunc(taxon, 16)), hjust = -0.08, size = 3, show.legend = FALSE
  ) +
  scale_y_reverse(breaks = 1:max(slope_taxa_long$rank), minor_breaks = NULL) +
  scale_x_discrete(expand = expansion(mult = 0.2)) +
  scale_colour_hue(l = 45, c = 120) +
  labs(x = NULL, y = "Rank (1 = most often dominant)", title = "Dominant taxon: First vs Last survey (slope graph)") +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    panel.grid.major.x = element_blank()
  )
# print(p_slope_gg)


# Optional: networkD3 Sankey (2 columns only)
# Left  = first-survey dominant species (ordered by first-survey importance, top 10)
# Right = rank 1..10 in last survey
# Link  = number of plots that move from each first-survey species to each last-survey rank
n_sankey <- 10L
if (requireNamespace("networkD3", quietly = TRUE)) {
  first_importance <- dominant_taxon_change %>%
    filter(!is.na(dominant_first)) %>%
    count(dominant_first, name = "n_first") %>%
    arrange(desc(n_first)) %>%
    slice_head(n = n_sankey)
  last_importance <- dominant_taxon_change %>%
    filter(!is.na(dominant_last)) %>%
    count(dominant_last, name = "n_last") %>%
    arrange(desc(n_last)) %>%
    slice_head(n = n_sankey)
  first_species <- first_importance$dominant_first
  last_species  <- last_importance$dominant_last
  rank_last_vec <- setNames(seq_along(last_species), last_species)

  flow_sankey <- dominant_taxon_change %>%
    filter(dominant_first %in% first_species, !is.na(dominant_last)) %>%
    mutate(rank_last = rank_last_vec[dominant_last]) %>%
    filter(!is.na(rank_last)) %>%
    count(dominant_first, rank_last, name = "n_plots") %>%
    filter(n_plots > 0)

  if (nrow(flow_sankey) > 0L) {
    left_nodes <- tibble(
      name = paste0(seq_along(first_species), ". ", first_species)
    )
    right_nodes <- tibble(
      name = paste0("Rank ", seq_len(n_sankey))
    )
    nodes_df <- bind_rows(left_nodes, right_nodes)

    links_df <- flow_sankey %>%
      mutate(
        source = match(dominant_first, first_species) - 1L,
        target = n_sankey + as.integer(rank_last) - 1L,
        value  = as.numeric(n_plots)
      ) %>%
      select(source, target, value)

    sankey_dominant <- networkD3::sankeyNetwork(
      Links = links_df,
      Nodes = nodes_df,
      Source = "source",
      Target = "target",
      Value = "value",
      NodeID = "name",
      units = "",
      fontSize = 14,
      nodeWidth = 28,
      sinksRight = TRUE,
      iterations = 0
    )
    if (requireNamespace("htmlwidgets", quietly = TRUE)) {
      sankey_dominant <- htmlwidgets::onRender(
        sankey_dominant,
        "
        function(el, x) {
          d3.select(el).selectAll('.link title')
            .text(function(d) {
              return d.source.name + ' -> ' + d.target.name + '\\n' + d.value + ' plots';
            });
          d3.select(el).selectAll('.node title')
            .text(function(d) {
              return d.name + '\\n' + d.value + ' plots';
            });
        }
        "
      )
    }
    print(sankey_dominant)
  }
}

# ------------------------------------------------------------------------------
# 6. Change in community composition (distance within plot, first vs last)
# ------------------------------------------------------------------------------
all_taxa <- wide_plot_year %>% select(-plot_number, -year) %>% names()
get_vec <- function(plot_num, yr) {
  r <- wide_plot_year %>% filter(plot_number == plot_num, year == yr)
  if (nrow(r) == 0) return(setNames(rep(0, length(all_taxa)), all_taxa))
  v <- r %>% select(-plot_number, -year) %>% as.matrix() %>% `[`(1, )
  full_vec <- setNames(rep(0, length(all_taxa)), all_taxa)
  full_vec[names(v)] <- v
  as.numeric(full_vec[all_taxa])
}

distances <- revisited %>%
  mutate(
    vec_first = map2(plot_number, year_first, get_vec),
    vec_last  = map2(plot_number, year_last, get_vec),
    mat_both  = map2(vec_first, vec_last, ~ rbind(.x, .y)),
    bc_dist   = map_dbl(mat_both, ~ as.numeric(vegdist(.x, method = "bray")))
  ) %>%
  select(plot_number, year_first, year_last, bc_dist)

comp_summary <- distances %>%
  summarise(
    n_plots = n(),
    mean_bc_dist   = mean(bc_dist, na.rm = TRUE),
    median_bc_dist = median(bc_dist, na.rm = TRUE),
    sd_bc_dist     = sd(bc_dist, na.rm = TRUE)
  )
print(comp_summary)

# ------------------------------------------------------------------------------
# 7. Visuals (pair with readout)
# ------------------------------------------------------------------------------
# Long format for diversity deltas
change_long <- change %>%
  select(plot_number, delta_S, delta_Shannon, delta_Simpson, delta_Evenness) %>%
  pivot_longer(cols = c(delta_S, delta_Shannon, delta_Simpson, delta_Evenness),
               names_to = "index", values_to = "delta") %>%
  mutate(index = factor(index, levels = c("delta_S", "delta_Shannon", "delta_Simpson", "delta_Evenness"),
                        labels = c("Δ Species richness (S)", "Δ Shannon", "Δ Simpson", "Δ Evenness")))

# 1) Distribution of within-plot change (diversity)
p_delta_diversity <- ggplot(change_long, aes(x = delta)) +
  geom_histogram(bins = 25, fill = "steelblue", alpha = 0.8, colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  facet_wrap(vars(index), scales = "free", ncol = 2) +
  labs(x = "Within-plot change (last − first)", y = "Number of plots",
       title = "Distribution of diversity change (n = 167 revisited plots)") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"))
# print(p_delta_diversity)

# 2) First vs last (paired) – S and Shannon
p_first_last_S <- ggplot(change, aes(x = S_first, y = S_last)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(alpha = 0.6, size = 2) +
  labs(x = "Species richness (first survey)", y = "Species richness (last survey)",
       title = "Paired plots: first vs last survey (1:1 = no change)") +
  theme_minimal()
# print(p_first_last_S)

p_first_last_Shannon <- ggplot(change, aes(x = Shannon_first, y = Shannon_last)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(alpha = 0.6, size = 2) +
  labs(x = "Shannon (first survey)", y = "Shannon (last survey)",
       title = "Paired plots: first vs last survey (1:1 = no change)") +
  theme_minimal()
# print(p_first_last_Shannon)

# 3) Tail change: dominance
p_delta_dominance <- ggplot(tail_change, aes(x = delta_max_pct)) +
  geom_histogram(bins = 25, fill = "darkorange", alpha = 0.8, colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  labs(x = "Change in dominant taxon share (last − first)", y = "Number of plots",
       title = "Increase in dominance (max % cover per plot)") +
  theme_minimal()
# print(p_delta_dominance)

# 3b) Which taxa are driving the increase? (dominant at last when delta_max_pct > 0)
p_culprits <- culprits_with_delta %>%
  mutate(dominant_last = str_trunc(dominant_last, 35)) %>%
  slice_max(n_plots, n = 20) %>%
  ggplot(aes(x = fct_reorder(dominant_last, n_plots), y = n_plots)) +
  geom_col(fill = "darkorange", alpha = 0.85) +
  geom_text(aes(label = n_plots), hjust = -0.2, size = 3) +
  coord_flip() +
  labs(x = NULL, y = "Number of plots (dominance increased)",
       title = "Taxa dominant at last survey when dominance increased (top 20)") +
  theme_minimal(base_size = 10)
# print(p_culprits)

# 4) Composition distance (Bray–Curtis)
p_bc_dist <- ggplot(distances, aes(x = bc_dist)) +
  geom_histogram(bins = 25, fill = "forestgreen", alpha = 0.8, colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = median(distances$bc_dist), linetype = "dashed", colour = "darkgreen", linewidth = 0.8) +
  labs(x = "Bray–Curtis distance (first vs last)", y = "Number of plots",
       title = "Compositional change per plot (0 = identical, 1 = no overlap)") +
  theme_minimal()
# print(p_bc_dist)

# Optional: one-page summary (patchwork if available)
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  p_readout <- p_delta_diversity /
    (p_first_last_S + p_first_last_Shannon) /
    (p_delta_dominance + p_bc_dist) +
    plot_annotation(title = "Paired change: first vs last survey (167 plots)")
  # print(p_readout)
  # ggsave("paired_change_readout.png", p_readout, width = 10, height = 11, dpi = 150)
}

# ------------------------------------------------------------------------------
# Outputs
# ------------------------------------------------------------------------------
# change         = per-plot deltas (delta_S, delta_Shannon, etc.)
# change_summary = mean/median of those deltas
# tail_change    = per-plot change in dominance (max_pct, n_high)
# tail_summary   = mean/median tail change
# distances      = per-plot Bray-Curtis between first and last survey
# comp_summary   = mean/median composition distance
# dominant_taxon_change, culprits_with_delta, gained_dominance (section 5b)
# Plots: p_delta_diversity, p_first_last_S, p_first_last_Shannon, p_delta_dominance, p_culprits, p_slope_gg, p_bc_dist
