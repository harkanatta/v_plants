# Pre-analysis checklist for two-wave paired vegetation data
# Purpose: run integrity checks BEFORE inferential models.
# Output style: PASS / FLAG tables + optional core tests.

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(vegan)
})

# ------------------------------------------------------------------------------
# 0) Settings
# ------------------------------------------------------------------------------
base_dir <- getwd()
plants_path <- file.path(base_dir, "plants_export", "v_plants_allt.csv")
stopifnot(file.exists(plants_path))

# If TRUE, subset by heiti_include. Supply heiti_include vector below.
use_heiti_include <- FALSE
heiti_include <- character(0)

# Core variables to audit for missingness and model readiness
core_vars <- c(
  "plot_number", "field_date", "taxon_name",
  "scale_value_code", "scale_value_value",
  "habitat_type_name", "moisture_name", "soil_type_name",
  "altitude", "slope", "precipitation", "temperature_year"
)

# ------------------------------------------------------------------------------
# 1) Load + date/year parsing
# ------------------------------------------------------------------------------
v_plants <- readr::read_csv(
  plants_path,
  show_col_types = FALSE,
  locale = readr::locale(encoding = "UTF-8")
)

if (use_heiti_include) {
  if (length(heiti_include) == 0) stop("use_heiti_include=TRUE but heiti_include is empty.")
  v_plants <- v_plants %>% filter(plot_number %in% stringr::str_trim(heiti_include))
}

v_plants2 <- v_plants %>%
  mutate(
    field_date_parsed = parse_date_time(
      field_date,
      orders = c("Y-m-d", "d-m-Y", "d/m/Y", "Y/m/d", "d.m.Y", "Y.m.d"),
      quiet = TRUE
    ),
    year = year(field_date_parsed)
  ) %>%
  filter(!is.na(plot_number), !is.na(taxon_name), !is.na(year))

# ------------------------------------------------------------------------------
# 2) Pairing integrity checks
# ------------------------------------------------------------------------------
plot_years <- v_plants2 %>%
  distinct(plot_number, year) %>%
  arrange(plot_number, year) %>%
  group_by(plot_number) %>%
  summarise(
    n_years = n(),
    year_first = first(year),
    year_last = last(year),
    years = str_c(year, collapse = ", "),
    .groups = "drop"
  )

revisited <- plot_years %>% filter(n_years == 2)

check_pairing <- tibble(
  check = c(
    "Input rows > 0",
    "Rows retained after year parsing > 0",
    "At least one revisited plot (n_years==2)",
    "No revisited plot with identical year_first/year_last"
  ),
  value = c(
    nrow(v_plants),
    nrow(v_plants2),
    nrow(revisited),
    sum(revisited$year_first == revisited$year_last, na.rm = TRUE)
  ),
  pass = c(
    nrow(v_plants) > 0,
    nrow(v_plants2) > 0,
    nrow(revisited) > 0,
    sum(revisited$year_first == revisited$year_last, na.rm = TRUE) == 0
  )
)

# Duplicate check after aggregation level used for community matrix
cover_dat <- v_plants2 %>%
  inner_join(revisited %>% select(plot_number, year_first, year_last), by = "plot_number") %>%
  mutate(
    sample = case_when(
      year == year_first ~ "first",
      year == year_last ~ "last",
      TRUE ~ NA_character_
    ),
    cover_num = coalesce(
      scale_value_value,
      case_when(
        scale_value_code %in% c("•", "·", "\u2022") ~ 0.25,
        scale_value_code == "+" ~ 0.75,
        scale_value_code == "1" ~ 3,
        scale_value_code == "2" ~ 15,
        scale_value_code == "3" ~ 37.5,
        scale_value_code == "4" ~ 62.5,
        scale_value_code == "5" ~ 87.5,
        TRUE ~ NA_real_
      )
    )
  ) %>%
  filter(!is.na(sample), !is.na(taxon_name), !is.na(cover_num), cover_num > 0)

dup_count <- cover_dat %>%
  count(plot_number, sample, taxon_name) %>%
  filter(n > 1) %>%
  nrow()

check_duplicates <- tibble(
  check = "No duplicates at plot_number x sample x taxon_name level",
  value = dup_count,
  pass = dup_count == 0
)

# ------------------------------------------------------------------------------
# 3) Missingness audit
# ------------------------------------------------------------------------------
present_core <- intersect(core_vars, names(v_plants2))
missing_overall <- v_plants2 %>%
  summarise(across(all_of(present_core), ~ mean(is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "prop_missing") %>%
  arrange(desc(prop_missing))

# Missingness shift between first and last sample among revisited plots
sample_meta <- v_plants2 %>%
  inner_join(revisited %>% select(plot_number, year_first, year_last), by = "plot_number") %>%
  mutate(
    sample = case_when(
      year == year_first ~ "first",
      year == year_last ~ "last",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(sample))

missing_by_sample <- map_dfr(present_core, function(v) {
  x <- sample_meta %>%
    mutate(is_miss = is.na(.data[[v]])) %>%
    group_by(sample) %>%
    summarise(prop_miss = mean(is_miss), n = n(), .groups = "drop") %>%
    pivot_wider(names_from = sample, values_from = c(prop_miss, n), values_fill = 0)
  tibble(
    variable = v,
    miss_first = x$prop_miss_first %||% NA_real_,
    miss_last = x$prop_miss_last %||% NA_real_,
    delta_pp = 100 * ((x$prop_miss_last %||% 0) - (x$prop_miss_first %||% 0))
  )
}) %>% arrange(desc(abs(delta_pp)))

# Broader variable-level profile over all columns
var_profile <- v_plants2 %>%
  summarise(across(
    everything(),
    list(
      prop_missing = ~mean(is.na(.x)),
      n_unique = ~n_distinct(.x, na.rm = TRUE)
    ),
    .names = "{.col}__{.fn}"
  )) %>%
  pivot_longer(
    everything(),
    names_to = c("variable", ".value"),
    names_sep = "__"
  ) %>%
  mutate(
    type = map_chr(variable, ~ class(v_plants2[[.x]])[1]),
    prop_missing = as.numeric(prop_missing),
    n_unique = as.numeric(n_unique),
    recommended_role = case_when(
      variable %in% c(
        "id", "plot_id", "subplot_id", "subplot_number",
        "group_type_subplot_id", "project_id", "project_number",
        "analyser_id", "tag_id"
      ) ~ "id / index (not model term)",
      variable %in% c("plot_number") ~ "plot identifier (pairing / strata)",
      variable %in% c("field_date", "field_date_parsed", "year") ~ "time / pairing (not covariate)",
      variable %in% c(
        "taxon_id", "taxon_name",
        "scale_value_code", "scale_value_value",
        "tag_name"
      ) ~ "species / cover core",
      variable %in% c(
        "comments", "subplot_comments",
        "project_descrtiption", "images", "subplot_images"
      ) ~ "free text / notes",
      variable %in% c("longitude_1", "latitude_1", "longitude_2", "latitude_2") ~ "spatial coordinate",
      variable %in% c(
        "precipitation", "temperature_year",
        "temperature_january", "temperature_july", "slope"
      ) & prop_missing >= 0.9 ~ "exclude: effectively not measured",
      prop_missing >= 0.95 ~ "exclude: ~all missing",
      prop_missing >= 0.5 ~ "sparse: only use in targeted subsets",
      type %in% c("numeric", "integer", "double") ~ "candidate numeric covariate / PCA",
      type %in% c("character", "factor") & n_unique <= 20 ~ "candidate factor covariate / stratifier",
      type %in% c("character", "factor") & n_unique > 20 ~ "high-cardinality label / not a direct model term",
      type %in% c("logical") ~ "logical indicator (candidate factor)",
      TRUE ~ "inspect manually"
    )
  ) %>%
  arrange(desc(prop_missing), variable)

# ------------------------------------------------------------------------------
# 4) Subplot comparability sanity (within plot: first vs last overlap)
# ------------------------------------------------------------------------------
subplots_long <- v_plants2 %>%
  inner_join(revisited %>% select(plot_number, year_first, year_last), by = "plot_number") %>%
  mutate(
    sample = case_when(
      year == year_first ~ "first",
      year == year_last ~ "last",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(sample), !is.na(subplot_id)) %>%
  distinct(plot_number, sample, subplot_id)

subplot_sets <- subplots_long %>%
  group_by(plot_number, sample) %>%
  summarise(ids = list(unique(subplot_id)), .groups = "drop") %>%
  pivot_wider(names_from = sample, values_from = ids)

subplot_overlap <- subplot_sets %>%
  mutate(
    n_first = map_int(first, length),
    n_last = map_int(last, length),
    n_intersect = map2_int(first, last, ~ length(intersect(.x, .y))),
    n_union = map2_int(first, last, ~ length(union(.x, .y))),
    jaccard = if_else(n_union > 0, n_intersect / n_union, NA_real_)
  ) %>%
  select(plot_number, n_first, n_last, n_intersect, n_union, jaccard)

check_subplot <- tibble(
  check = c(
    "Median subplot Jaccard overlap >= 0.5",
    "Plots with zero subplot overlap <= 10%"
  ),
  value = c(
    median(subplot_overlap$jaccard, na.rm = TRUE),
    mean(subplot_overlap$n_intersect == 0, na.rm = TRUE)
  ),
  pass = c(
    median(subplot_overlap$jaccard, na.rm = TRUE) >= 0.5,
    mean(subplot_overlap$n_intersect == 0, na.rm = TRUE) <= 0.10
  )
)

# ------------------------------------------------------------------------------
# 5) Build paired community matrix for multivariate tests
# ------------------------------------------------------------------------------
community <- cover_dat %>%
  group_by(plot_number, sample, taxon_name) %>%
  summarise(cover = mean(cover_num), .groups = "drop") %>%
  pivot_wider(names_from = taxon_name, values_from = cover, values_fill = 0)

# Keep complete pairs only (exactly first+last rows)
pair_counts <- community %>% count(plot_number, name = "n_samples")
community <- community %>%
  inner_join(pair_counts %>% filter(n_samples == 2), by = "plot_number") %>%
  select(-n_samples)

meta <- community %>% select(plot_number, sample)
mat <- community %>% select(-plot_number, -sample) %>% as.matrix()

check_multivar <- tibble(
  check = c(
    "Community matrix has rows",
    "Community matrix has columns (taxa)",
    "Exactly two rows per plot_number in matrix"
  ),
  value = c(
    nrow(mat),
    ncol(mat),
    all((meta %>% count(plot_number) %>% pull(n)) == 2)
  ),
  pass = c(
    nrow(mat) > 0,
    ncol(mat) > 0,
    all((meta %>% count(plot_number) %>% pull(n)) == 2)
  )
)

# ------------------------------------------------------------------------------
# 6) Pre-tests for PERMANOVA assumptions + core paired PERMANOVA
# ------------------------------------------------------------------------------
dist_bray <- vegdist(mat, method = "bray")

bd <- betadisper(dist_bray, group = meta$sample)
bd_test <- permutest(bd, permutations = 999)

perm <- adonis2(
  mat ~ sample,
  data = meta,
  method = "bray",
  permutations = 999,
  strata = meta$plot_number
)

# ------------------------------------------------------------------------------
# 7) Print checklist output
# ------------------------------------------------------------------------------
cat("\n=== PRE-ANALYSIS CHECKLIST ===\n")
print(bind_rows(check_pairing, check_duplicates, check_subplot, check_multivar) %>%
        mutate(status = if_else(pass, "PASS", "FLAG")) %>%
        select(status, check, value))

cat("\n=== MISSINGNESS (OVERALL) ===\n")
print(missing_overall, n = min(30, nrow(missing_overall)))

cat("\n=== MISSINGNESS SHIFT FIRST vs LAST (percentage points) ===\n")
print(missing_by_sample, n = min(30, nrow(missing_by_sample)))

cat("\n=== SUBPLOT OVERLAP SUMMARY ===\n")
print(subplot_overlap %>%
        summarise(
          n_plots = n(),
          median_jaccard = median(jaccard, na.rm = TRUE),
          mean_jaccard = mean(jaccard, na.rm = TRUE),
          prop_zero_overlap = mean(n_intersect == 0, na.rm = TRUE)
        ))

cat("\n=== BETADISPER TEST (homogeneity of dispersion by sample) ===\n")
print(bd_test)

cat("\n=== PAIRED PERMANOVA (adonis2 with strata = plot_number) ===\n")
print(perm)

cat("\nDone.\n")

