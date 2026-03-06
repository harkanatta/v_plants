# Diversity indices from vegetation data (tegundir.txt, subplot.txt)
# Shannon, Margalef, Simpson, Evenness (Pielou), Berger–Parker, inverse Simpson
# Requires: tidyverse, vegan (optional but recommended)

library(tidyverse)

# ------------------------------------------------------------------------------
# 1. Read data
# ------------------------------------------------------------------------------
# Paths: run with setwd() to project dir, or set base_dir explicitly
base_dir <- getwd()

tegundir_path <- file.path(base_dir, "tegundir.txt")
subplot_path  <- file.path(base_dir, "subplot.txt")

stopifnot("tegundir.txt not found" = file.exists(tegundir_path))

# Read Cover as character so "3"/"4" are not parsed as numeric (we map them via Braun-Blanquet)
tegundir <- read_tsv(tegundir_path, locale = locale(encoding = "UTF-8"),
                     col_types = cols(Cover = col_character(), .default = col_guess()),
                     show_col_types = FALSE)
subplot  <- if (file.exists(subplot_path)) {
  read_tsv(subplot_path, locale = locale(encoding = "UTF-8"), show_col_types = FALSE)
} else {
  tibble()
}

# ------------------------------------------------------------------------------
# 2. Parse Cover to a single numeric (%; mixed input: ranges, Braun-Blanquet, numbers)
# Braun-Blanquet: • 0–0.5% → 0.25; + 0.5–1% → 0.75; 1→3, 2→15, 3→37.5, 4→62.5, 5→87.5
# Ranges "low-high" → midpoint; other numbers → as-is (percent).
# ------------------------------------------------------------------------------
cover_to_numeric <- function(x) {
  # Braun-Blanquet: symbol/code -> midpoint of % range
  bb <- c(
    "•" = 0.25,   # 0–0.5%  (U+2022 bullet; other bullets handled below)
    "+" = 0.75,   # 0.5–1%
    "1" = 3,      # 1–5%
    "2" = 15,     # 5–25%
    "3" = 37.5,   # 25–50%
    "4" = 62.5,   # 50–75%
    "5" = 87.5    # 75–100%
  )
  out <- rep(NA_real_, length(x))
  for (i in seq_along(x)) {
    s <- as.character(x[i])
    s <- trimws(s)
    # Normalise "3.0" -> "3" so BB classes 1-5 are recognised when read as numeric
    s <- sub("\\.0+$", "", s)
    if (is.na(s) || !nzchar(s)) next
    # Range "number-number"
    if (grepl("^[0-9.]+-[0-9.]+$", s)) {
      r <- as.numeric(strsplit(s, "-", fixed = TRUE)[[1]])
      out[i] <- mean(r, na.rm = TRUE)
      next
    }
    # Braun-Blanquet: bullet (common Unicode variants)
    if (s %in% c("\u2022", "\u2023", "\u2219", "\u00b7", "•", "·")) {
      out[i] <- 0.25
      next
    }
    # Braun-Blanquet: + or digit 1–5 only
    if (s == "+") {
      out[i] <- 0.75
      next
    }
    if (s %in% c("1", "2", "3", "4", "5")) {
      out[i] <- as.numeric(bb[s])   # BB class -> midpoint %
      next
    }
    # Single numeric (percent, e.g. 25 or 7.5)
    if (grepl("^[0-9.]+$", s)) {
      out[i] <- as.numeric(s)
      next
    }
  }
  out
}

tegundir <- tegundir %>%
  mutate(cover_num = cover_to_numeric(Cover))

# Drop rows with no usable cover (optional: treat as presence = 1 for richness-only indices)
tegundir <- tegundir %>%
  filter(!is.na(cover_num), cover_num > 0)

# ------------------------------------------------------------------------------
# 3. Diversity index functions
# ------------------------------------------------------------------------------
# x = vector of abundances (or cover) per species in one sample
# All indices use proportions p_i = x_i / sum(x); N = sum(x), S = number of species

shannon <- function(x) {
  p <- x / sum(x)
  p <- p[p > 0]
  -sum(p * log(p))
}

simpson_lambda <- function(x) {
  p <- x / sum(x)
  sum(p^2)
}

simpson_inverse <- function(x) 1 / simpson_lambda(x)   # 1/lambda
simpson_1_minus  <- function(x) 1 - simpson_lambda(x)  # 1 - lambda (Gini–Simpson)

margalef <- function(x) {
  S <- length(x)
  N <- sum(x)
  if (N <= 0) return(NA_real_)
  (S - 1) / log(N)
}

pielou_evenness <- function(x) {
  S <- length(x)
  if (S <= 0) return(NA_real_)
  H <- shannon(x)
  H / log(S)
}

berger_parker <- function(x) {
  p <- x / sum(x)
  max(p)
}

# ------------------------------------------------------------------------------
# 4. Compute indices per subplot
# ------------------------------------------------------------------------------
by_subplot <- tegundir %>%
  group_by(Subplot, Taxon) %>%
  summarise(cover = sum(cover_num), .groups = "drop") %>%
  group_by(Subplot) %>%
  summarise(
    abundances = list(cover),
    .groups = "drop"
  )

div_calc <- function(ab) {
  ab <- ab[[1]]
  if (length(ab) == 0 || sum(ab) <= 0) {
    return(tibble(
      shannon = NA_real_, margalef = NA_real_, simpson_1minus = NA_real_,
      simpson_inverse = NA_real_, evenness = NA_real_, berger_parker = NA_real_,
      richness = 0L, total_cover = 0
    ))
  }
  tibble(
    shannon         = shannon(ab),
    margalef        = margalef(ab),
    simpson_1minus  = simpson_1_minus(ab),
    simpson_inverse = simpson_inverse(ab),
    evenness        = pielou_evenness(ab),
    berger_parker   = berger_parker(ab),
    richness        = length(ab),
    total_cover     = sum(ab)
  )
}

diversity <- by_subplot %>%
  mutate(indices = map(abundances, div_calc)) %>%
  unnest(indices) %>%
  select(-abundances)

# Optional: derive a short site/plot label from Subplot (e.g. "GD-S1", "Lh-S1")
diversity <- diversity %>%
  mutate(
    site = str_trim(str_replace(Subplot, "^[0-9]+ - [^-]+ - ([^-]+) - .*", "\\1"))
  )

# ------------------------------------------------------------------------------
# 5. Join subplot metadata if available
# ------------------------------------------------------------------------------
if (nrow(subplot) > 0L && "Plot" %in% names(subplot)) {
  # Match Subplot to Plot + subplot number if needed for extra grouping
  subplot_key <- subplot %>%
    transmute(
      plot_id = Plot,
      subplot_num = `Subplot number`,
      total_cover_pct = `Total cover (%)`
    )
  # Simple join by matching start of Subplot to Plot (if structure allows)
  diversity <- diversity %>%
    mutate(plot_id = str_replace(Subplot, " - [0-9]+$", "")) %>%
    left_join(
      subplot_key %>% distinct(plot_id = plot_id, .keep_all = TRUE),
      by = "plot_id"
    )
}

# ------------------------------------------------------------------------------
# 6. Visualisations (tidyverse / ggplot2)
# ------------------------------------------------------------------------------
# Long format for faceted plots
diversity_long <- diversity %>%
  select(Subplot, site, shannon, margalef, simpson_1minus, simpson_inverse,
         evenness, berger_parker, richness) %>%
  pivot_longer(
    cols = c(shannon, margalef, simpson_1minus, simpson_inverse, evenness, berger_parker, richness),
    names_to = "index",
    values_to = "value"
  ) %>%
  filter(!is.na(value))

# Labels for plots
index_lab <- c(
  shannon = "Shannon H'",
  margalef = "Margalef D",
  simpson_1minus = "Simpson (1−λ)",
  simpson_inverse = "Inverse Simpson",
  evenness = "Pielou J (evenness)",
  berger_parker = "Berger–Parker (dominance)",
  richness = "Species richness"
)
diversity_long <- diversity_long %>%
  mutate(index_lab = factor(index, levels = names(index_lab), labels = index_lab))

# --- Plot 1: Boxplots of each index by site ---
p_by_site <- ggplot(diversity_long %>% filter(index != "richness"),
                    aes(x = fct_infreq(site), y = value, fill = site)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.6) +
  facet_wrap(vars(index_lab), scales = "free_y", ncol = 2) +
  labs(x = "Site", y = "Index value", title = "Diversity indices by site") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )
ggsave(file.path(base_dir, "diversity_by_site.png"), p_by_site,
       width = 10, height = 9, dpi = 150)
print(p_by_site)

# --- Plot 2: Richness by site (bar) ---
p_richness <- diversity %>%
  group_by(site) %>%
  summarise(mean_richness = mean(richness), .groups = "drop") %>%
  ggplot(aes(x = fct_reorder(site, mean_richness), y = mean_richness, fill = site)) +
  geom_col(alpha = 0.8) +
  labs(x = "Site", y = "Mean species richness", title = "Mean species richness per subplot by site") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
ggsave(file.path(base_dir, "richness_by_site.png"), p_richness, width = 8, height = 5, dpi = 150)
print(p_richness)

# --- Plot 3: Shannon vs Evenness (scatter) ---
p_shannon_evenness <- ggplot(diversity, aes(x = shannon, y = evenness, colour = site)) +
  geom_point(alpha = 0.8, size = 2) +
  labs(x = "Shannon H'", y = "Pielou evenness J'",
       title = "Shannon diversity vs evenness by subplot") +
  theme_minimal()
ggsave(file.path(base_dir, "shannon_vs_evenness.png"), p_shannon_evenness,
       width = 7, height = 5, dpi = 150)
print(p_shannon_evenness)

# --- Plot 4: All indices distributions (histograms / density) ---
p_dist <- ggplot(diversity_long %>% filter(!index %in% c("richness")),
                 aes(x = value)) +
  geom_histogram(bins = 25, fill = "steelblue", alpha = 0.7) +
  facet_wrap(vars(index_lab), scales = "free", ncol = 2) +
  labs(x = "Value", y = "Count", title = "Distribution of diversity indices") +
  theme_minimal(base_size = 11)
ggsave(file.path(base_dir, "diversity_distributions.png"), p_dist,
       width = 10, height = 8, dpi = 150)
print(p_dist)

# --- Plot 5: Heatmap-style tile (sites x indices, mean value) ---
heat_data <- diversity_long %>%
  filter(index != "richness") %>%
  group_by(site, index_lab) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")

p_heat <- ggplot(heat_data, aes(x = index_lab, y = fct_infreq(site), fill = mean_value)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  scale_fill_viridis_c(option = "viridis", na.value = "grey90") +
  labs(x = "Index", y = "Site", fill = "Mean value",
       title = "Mean diversity indices by site (heatmap)") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
ggsave(file.path(base_dir, "diversity_heatmap.png"), p_heat, width = 9, height = 6, dpi = 150)
print(p_heat)

# ------------------------------------------------------------------------------
# 7. Optional: use vegan for Shannon/Simpson (same results)
# ------------------------------------------------------------------------------
if (requireNamespace("vegan", quietly = TRUE)) {
  # Verify: vegan::diversity(ab, "shannon") == shannon(ab)
  # diversity already computed above; vegan can be used in div_calc if preferred
  message("vegan is installed; diversity() and diversity(x, 'simpson') match our functions.")
}

# ------------------------------------------------------------------------------
# 8. Write summary table
# ------------------------------------------------------------------------------
diversity %>%
  select(Subplot, site, shannon, margalef, simpson_1minus, evenness, berger_parker, richness) %>%
  write_csv(file.path(base_dir, "diversity_summary.csv"), na = "")

message("Done. Outputs: diversity_summary.csv, diversity_by_site.png, richness_by_site.png, ",
        "shannon_vs_evenness.png, diversity_distributions.png, diversity_heatmap.png")
