# Keep only rows where Subplot contains any of the listed Heiti (site codes).
# Run from project dir, or set base_dir. Expects tegundir in workspace or reads from file.

library(tidyverse)

base_dir <- getwd()

# ------------------------------------------------------------------------------
# Heiti to include (Subplot must contain at least one of these)
# ------------------------------------------------------------------------------
heiti_include <- c(
  "AL-10-01", "AL-24-03", "AL-40-01", "AL-40-03", "AL-42-04", "AL-43-01",
  "EY-127-02", "EY-127-04", "FV01-1", "HF-65-01", "HF-65-02", "HF-65-03",
  "HF-65-05", "HF-65-06", "HF-65-07", "HF-65-08", "HF-65-09", "HF-65-13",
  "HRM01", "HRM02", "HRM04", "HRM06", "ID-12", "ID-19", "KSM17",
  "MM-10-01", "MM-110-01", "MS-42-01", "MS-42-02", "MS-65-01", "MS-65-03",
  "MS-65-04", "NV-40-01", "NV-41-01", "NV-41-06", "NV-41-08", "NV-42-05",
  "OF-01", "OF-02", "OF-04", "OX-01-2022", "OX-02-2022", "RN-04",
  "RN-65-01", "RN-65-02", "RN-65-02+", "RN-65-02X", "RN-65-03",
  "S30-1", "SF-66-01", "SK-EY-02", "SL-41-01", "SL-41-02", "SL-42-05",
  "SL42-06", "SL43-04", "T02-3", "T05-1", "T05-2", "T05-3", "T05-4",
  "TH-40-02", "TH-43-01", "TH-MYV-01-2022", "TH-MYV-02-2022", "U16",
  "V03-2", "V05-4", "V13-2", "VF 40-01", "VF-03", "VF-06", "VF-10-08",
  "VF-11", "VF-40-02", "VF41-09", "VF59-03", "VF-65-02", "VF-65-03",
  "VF65-05", "VF65-05AA", "VF-65-06", "AEY-01", "AEY-02", "AEY03",
  "G03-1", "G03-3", "G15-3", "G20-3", "G20-4", "HEY-01", "HMM07", "HMM08",
  "HRM07", "KKB39", "KMM35", "KMM68", "KST52", "KVV84", "MS-30-02",
  "MS-40-01", "MS-40-03", "MS-42-03", "MS-42-04", "MY-78-01", "RN-05-02",
  "RN-07-02", "RN-07-03", "RN-24-03", "SN-42-05", "SN-59-04", "SN-59-02Y",
  "SN-59-05X", "V-02-3", "V06-1", "V06-4", "V07-2", "V07-3", "V09-1",
  "V09-2", "V14-3", "V15-3", "V16-2", "V16-4", "V18-1", "VF 59-04",
  "VF08", "VF59-01", "E07-1", "E07-3", "G03-2", "G03-5", "G09-1",
  "HRM09", "HRM10", "HST06", "KBA05", "KKB08", "KMM14", "KSD39", "KVV00",
  "MY-100-01", "MY-650-03X", "NV-42-01", "OF-03", "RN-01", "RN-02",
  "RN-03", "RN-04-13", "RN-04-14", "RN-15-03", "RN-23-02", "RN-24-02",
  "RN-25-04", "RN-25-07", "S21-1", "S22-1", "S27-2", "SL-80-04",
  "SN-41-04B", "SN-53-02", "T02-2", "TO5-5", "*VF16", "VF20-01", "VF21-10",
  "VF22-02", "VF30-01", "MY-05-03B", "MY-31-03", "MY-37-01", "MY-69-04",
  "RN-17-03", "RN-23-04", "RN-28-03", "S21-2", "S21-3", "S25-2", "S26-1",
  "S27-1", "SK-31-03", "SK-32-04", "SK-58-01", "SN-40-01", "FG2-1",
  "FA9-1", "FT10-3", "FG2-5", "FD1-8", "L10-1", "L13-5", "V02-1",
  "V02-2", "V02-4", "V06-2", "V07-1", "V14-1", "V16-3", "VF-12", "VF-04",
  "VF-05", "VF-41-07", "VF-40-06", "NV-40-02", "NV-40-03", "NV-40-04",
  "NV-40-05", "NV-40-08"
) %>% str_trim()

# ------------------------------------------------------------------------------
# Load data (use existing tegundir if present, else read from file)
# ------------------------------------------------------------------------------
if (!exists("tegundir")) {
  p <- file.path(base_dir, "tegundir.txt")
  if (!file.exists(p)) stop("tegundir.txt not found and tegundir not in workspace.")
  tegundir <- read_tsv(p, locale = locale(encoding = "UTF-8"),
                      col_types = cols(Cover = col_character(), .default = col_guess()),
                      show_col_types = FALSE)
}

# ------------------------------------------------------------------------------
# Build include pattern (literal substring match; escape regex metacharacters)
# ------------------------------------------------------------------------------
# Escape regex metacharacters so heiti are matched literally ( ] and [ must be escaped in class)
escape_regex <- function(x) str_replace_all(x, "([.+*?^$()\\[\\]|\\\\])", "\\\\\\1")
pattern_include <- heiti_include %>% escape_regex() %>% str_c(collapse = "|")

# ------------------------------------------------------------------------------
# Cover -> single numeric % (Braun-Blanquet: • 0.25, + 0.75, 1→3, 2→15, 3→37.5, 4→62.5, 5→87.5)
# ------------------------------------------------------------------------------
cover_to_numeric <- function(x) {
  bb <- c("•" = 0.25, "+" = 0.75, "1" = 3, "2" = 15, "3" = 37.5, "4" = 62.5, "5" = 87.5)
  out <- rep(NA_real_, length(x))
  for (i in seq_along(x)) {
    s <- trimws(sub("\\.0+$", "", as.character(x[i])))
    if (is.na(s) || !nzchar(s)) next
    if (grepl("^[0-9.]+-[0-9.]+$", s)) {
      r <- as.numeric(strsplit(s, "-", fixed = TRUE)[[1]])
      out[i] <- mean(r, na.rm = TRUE)
    } else if (s %in% c("\u2022", "\u2023", "\u2219", "\u00b7", "•", "·")) {
      out[i] <- 0.25
    } else if (s == "+") {
      out[i] <- 0.75
    } else if (s %in% c("1", "2", "3", "4", "5")) {
      out[i] <- as.numeric(bb[s])
    } else if (grepl("^[0-9.]+$", s)) {
      out[i] <- as.numeric(s)
    }
  }
  out
}

# ------------------------------------------------------------------------------
# Filter: keep only rows where Subplot contains at least one heiti
# ------------------------------------------------------------------------------
tegundir_filtered <- tegundir %>%
  filter(str_detect(Subplot, pattern_include)) %>%
  mutate(cover_num = cover_to_numeric(Cover)) %>%
  mutate(
    # Subplot format: "ID - Name - SITECODE - date - N" -> site = 3rd segment, subplot_no = last
    parts = str_split(Subplot, " - "),
    site = map_chr(parts, ~ if (length(.x) >= 3) .x[3] else NA_character_),
    subplot_no = map_chr(parts, ~ .x[length(.x)]),
    parts = NULL   # drop helper
  )

# ------------------------------------------------------------------------------
# Report
# ------------------------------------------------------------------------------
n_before <- nrow(tegundir)
n_after  <- nrow(tegundir_filtered)
message("Subplot filter (include heiti): ", n_before, " -> ", n_after, " rows (excluded ", n_before - n_after, ")")

matched <- heiti_include[map_lgl(heiti_include, ~ any(str_detect(tegundir$Subplot, fixed(.x))))]
message("Heiti with at least one match: ", length(matched), " of ", length(heiti_include))

# Optional: write filtered table
# write_tsv(tegundir_filtered, file.path(base_dir, "tegundir_filtered.txt"), na = "")


# ==============================================================================
# Diversity indices (after pooling subplots by site: mean cover per site × Taxon)
# ==============================================================================
# Pool subplots first (mean cover per site × Taxon), then one species matrix row per site.
library(vegan)

# Pool: one row per site × Taxon, cover = mean across subplots
pooled <- tegundir_filtered %>%
  filter(!is.na(cover_num), cover_num > 0) %>%
  group_by(site, Taxon) %>%
  summarise(cover = mean(cover_num), .groups = "drop")

# Species matrix: one row per site, columns = Taxon, values = cover
species_wide <- pooled %>%
  pivot_wider(names_from = Taxon, values_from = cover, values_fill = 0)

# Row identifiers and matrix (extract site first, then matrix from remaining columns)
sites <- species_wide$site
species_matrix <- species_wide %>%
  select(-site) %>%
  as.matrix()
rownames(species_matrix) <- sites

# Diversity indices per site
diversity_by_site <- tibble(
  site = sites,
  N = rowSums(species_matrix),
  S = specnumber(species_matrix),
  Shannon = diversity(species_matrix, index = "shannon"),
  Simpson = diversity(species_matrix, index = "simpson"),
  Margalef = (S - 1) / log(pmax(N, 1)),
  Evenness = Shannon / log(pmax(S, 1))
)

# Long format (index, Skor) like gagnavinnsla tafla
diversity_by_site_long <- diversity_by_site %>%
  select(site, S, Shannon, Simpson, Margalef, Evenness) %>%
  pivot_longer(cols = c(S, Shannon, Simpson, Margalef, Evenness),
               names_to = "index", values_to = "Skor")

# Plots
p_diversity_site <- ggplot(diversity_by_site_long, aes(x = fct_reorder(site, Skor, .fun = mean), y = Skor, fill = site)) +
  geom_col(alpha = 0.8) +
  facet_wrap(vars(index), scales = "free_y", ncol = 2) +
  labs(x = "Site", y = "Skor", title = "Diversity indices by site (subplots pooled, mean cover)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
# print(p_diversity_site)
# ggsave(file.path(base_dir, "diversity_by_site_veg.png"), p_diversity_site, width = 10, height = 8, dpi = 150)

p_richness_site <- ggplot(diversity_by_site_long %>% filter(index == "S"),
                          aes(x = fct_reorder(site, Skor), y = Skor, fill = site)) +
  geom_col(alpha = 0.8) +
  labs(x = "Site", y = "Species richness (S)", title = "Species richness by site") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
# print(p_richness_site)

message("Diversity objects: pooled, species_wide, species_matrix, diversity_by_site, diversity_by_site_long, p_diversity_site, p_richness_site")
