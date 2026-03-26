suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(vegan)
})

set.seed(123)

# ------------------------------------------------------------------------------
# 0) Settings
# ------------------------------------------------------------------------------
base_dir_candidates <- c(
  getwd(),
  normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = FALSE)
)
plants_candidates <- file.path(base_dir_candidates, "plants_export", "v_plants_allt.csv")
plants_exists <- file.exists(plants_candidates)
if (!any(plants_exists)) {
  stop(
    "Could not locate plants_export/v_plants_allt.csv.\nTried:\n",
    paste0(" - ", plants_candidates, collapse = "\n")
  )
}
plants_path <- plants_candidates[which(plants_exists)[1]]
base_dir <- dirname(dirname(plants_path))

# Optional species-exclusion sheet
# - CSV: v_plants_taxa_exclude_flags.csv (two columns, no header)
# - Column 1: exclusion flag (e.g. "x" = drop this taxon from all analyses)
# - Column 2: taxon_name as it appears in v_plants_raw$taxon_name
taxon_exclude_path <- file.path(base_dir, "v_plants_taxa_exclude_flags.csv")
taxa_exclude <- character(0)
if (file.exists(taxon_exclude_path)) {
  taxon_exclude_tbl <- readr::read_csv(
    taxon_exclude_path,
    col_names = c("exclude_flag", "taxon_name"),
    show_col_types = FALSE
  )
  
  # Build a simple character vector of taxon names to drop
  # (any non-empty flag in column 1 is treated as "exclude")
  taxa_exclude <- taxon_exclude_tbl %>%
    filter(!is.na(exclude_flag), exclude_flag != "") %>%
    pull(taxon_name) %>%
    unique()
}

# Optional column selection for environmental variables
# - CSV: v_plants_envvar_selection.csv (two columns, no header)
# - Column 1: flag (e.g. "x" = include this column as an environmental candidate)
# - Column 2: variable name as it appears in v_plants_raw / v_plants2
raw_cols_sheet_path <- file.path(base_dir, "v_plants_envvar_selection.csv")

vars_env_from_sheet <- character(0)
if (file.exists(raw_cols_sheet_path)) {
  raw_cols_tbl <- readr::read_csv(
    raw_cols_sheet_path,
    col_names = c("flag", "variable"),
    show_col_types = FALSE
  )
  
  # Take all variables where the flag column is non-empty
  # (this allows using "x", "1", etc. as a generic "selected" marker)
  vars_env_from_sheet <- raw_cols_tbl %>%
    filter(!is.na(flag), flag != "") %>%
    pull(variable) %>%
    unique()
}

# Special handling for vegetation height:
# if any of the four raw height columns are flagged in the sheet, we
# replace them with a single derived variable name "vegetation_height_mean"
# (the actual mean is computed later when building v_plants2).
height_vars_sheet <- c("vegetation_height_1", "vegetation_height_2",
                       "vegetation_height_3", "vegetation_height_4")
if (any(height_vars_sheet %in% vars_env_from_sheet)) {
  vars_env_from_sheet <- c(
    setdiff(vars_env_from_sheet, height_vars_sheet),
    "vegetation_height_mean"
  )
}

use_heiti_include <- TRUE
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
)

vars_exclude_missing <- c(
  "biomass", "precipitation", "stone_class", "temperature_january",
  "temperature_july", "temperature_year", "heildarthekja", "subplot_carbon",
  "subplot_nitrogen", "subplot_ph", "plant_community_other",
  "plant_community_3_code", "plant_community_3_id", "plant_community_3_name",
  "slope_aspect", "slope", "autt", "surface_type_code", "surface_type_id",
  "surface_type_name", "soil_temperature_2", "plant_community_cover", "sinuthekja"
)
vars_exclude_text <- c(
  "comments", "subplot_comments", "images", "subplot_images", "project_description"
)
vars_exclude_strict <- c(vars_exclude_missing, vars_exclude_text)

# ------------------------------------------------------------------------------
# 1) Load + filter + parse date
# ------------------------------------------------------------------------------
v_plants_raw <- readr::read_csv(
  plants_path,
  show_col_types = FALSE,
  locale = readr::locale(encoding = "UTF-8")
)

if (use_heiti_include) {
  if (length(heiti_include) == 0) {
    stop("use_heiti_include = TRUE but heiti_include is empty.")
  }
  v_plants_raw <- v_plants_raw %>%
    filter(plot_number %in% stringr::str_trim(heiti_include))
}

if (length(taxa_exclude) > 0) {
  v_plants_raw <- v_plants_raw %>%
    filter(!taxon_name %in% taxa_exclude)
}

v_plants_raw_all <- v_plants_raw

valid_subplots_main <- as.character(1:8)
outside_subplot_codes <- c("0", "Utan reita", "utan reita")

v_plants_raw_main <- v_plants_raw_all %>%
  filter(is.na(subplot_number) | subplot_number %in% valid_subplots_main)

v_plants2 <- v_plants_raw_main %>%
  mutate(
    field_date_parsed = parse_date_time(
      field_date,
      orders = c("Y-m-d", "d-m-Y", "d/m/Y", "Y/m/d", "d.m.Y", "Y.m.d"),
      quiet = TRUE
    ),
    year = year(field_date_parsed),
    vegetation_height_mean = {
      h <- cbind(
        vegetation_height_1,
        vegetation_height_2,
        vegetation_height_3,
        vegetation_height_4
      )
      ifelse(rowSums(!is.na(h)) == 0, NA_real_, rowMeans(h, na.rm = TRUE))
    },
    grytnithekja = replace_na(grytnithekja, 0)
  ) %>%
  filter(!is.na(plot_number), !is.na(taxon_name), !is.na(year))

v_plants_use <- v_plants2 %>%
  select(-any_of(vars_exclude_strict))


# Optional: table from QGIS (CSV) ??? NA/NULL replaced with 0, joined to v_plants_use and env/meta
my_table_path <- file.path(base_dir, "outputs", "qgis_share", "my_table.csv")
my_table <- NULL
if (file.exists(my_table_path)) {
  my_table <- readr::read_csv(
    my_table_path,
    show_col_types = FALSE,
    na = c("", "NA", "NULL")
  ) %>%
    mutate(
      across(where(is.numeric), ~ replace_na(.x, 0)),
      across(where(is.character), ~ replace_na(.x, ""))
    )
  if ("plot_number" %in% names(my_table)) {
    my_table <- my_table %>% distinct(plot_number, .keep_all = TRUE)
    v_plants_use <- v_plants_use %>%
      left_join(my_table, by = "plot_number")
  }
}

# ------------------------------------------------------------------------------
# 2) Paired design + community matrix
# ------------------------------------------------------------------------------
plot_years <- v_plants2 %>%
  distinct(plot_number, year) %>%
  arrange(plot_number, year) %>%
  group_by(plot_number) %>%
  summarise(
    n_years = n(),
    year_first = first(year),
    year_last = last(year),
    .groups = "drop"
  )

revisited <- plot_years %>% filter(n_years >= 2)

paired_long <- v_plants2 %>%
  inner_join(revisited %>% select(plot_number, year_first, year_last), by = "plot_number") %>%
  mutate(
    sample = case_when(
      year == year_first ~ "first",
      year == year_last ~ "last",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(sample))

paired_cover <- paired_long %>%
  mutate(
    cover_num = coalesce(
      scale_value_value,
      case_when(
        scale_value_code %in% c("???", "??", "\u2022") ~ 0.25,
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
  filter(!is.na(taxon_name), !is.na(cover_num), cover_num > 0)

community_long <- paired_cover %>%
  group_by(plot_number, sample, taxon_name) %>%
  summarise(cover = mean(cover_num), .groups = "drop")

pair_counts <- community_long %>%
  group_by(plot_number) %>%
  summarise(n_samples = n_distinct(sample), .groups = "drop")

community_long_paired <- community_long %>%
  inner_join(pair_counts %>% filter(n_samples == 2), by = "plot_number") %>%
  select(-n_samples)

community_wide <- community_long_paired %>%
  pivot_wider(names_from = taxon_name, values_from = cover, values_fill = 0)

meta <- community_wide %>% select(plot_number, sample)
mat <- community_wide %>% select(-plot_number, -sample) %>% as.matrix()

# ------------------------------------------------------------------------------
# 2a) Per-plot Bray???Curtis distance (first vs last) + gap in years
# ------------------------------------------------------------------------------
bc_per_plot <- community_wide %>%
  group_by(plot_number) %>%
  summarise(
    bc = {
      mat2 <- pick(where(is.numeric)) %>% as.matrix()
      if (nrow(mat2) == 2) {
        as.numeric(vegan::vegdist(mat2, method = "bray")[1])
      } else {
        NA_real_
      }
    },
    .groups = "drop"
  )

gap_tbl <- plot_years %>%
  transmute(plot_number, gap_years = year_last - year_first)

bc_with_gap <- bc_per_plot %>%
  left_join(gap_tbl, by = "plot_number")

# ------------------------------------------------------------------------------
# 2b) Plot-level diversity (all indices from combined subplot + outside cover)
# ------------------------------------------------------------------------------
# Design:
#   - Subplot species (subplot_number 1–8): measured BB cover, treated as the
#     best estimate of plot-level cover density.
#   - Outside-subplot species (subplot_number in outside_subplot_codes): recorded
#     as present in the plot but not within any subplot. Assigned MIN_BB_COVER
#     (the lowest Braun–Blanquet value, 0.25 %) to enter the abundance vector.
#   - S, Shannon H', Simpson 1–D, Pielou's J: all computed from this combined
#     cover set in a single pass.
#
# Geometry (for documentation / sanity):
#   Subplot:  33 cm × 100 cm = 0.33 m²
#   Plot:     400 m²  (always)
#   Subplots per plot: 6 or 8  → total subplot area 1.98 or 2.64 m²
#   Subplot cover is taken as representative of plot-level density; outside-
#   subplot species receive MIN_BB_COVER as a conservative presence indicator.

PLOT_AREA_M2    <- 400
SUBPLOT_W_CM    <- 33
SUBPLOT_L_CM    <- 100
SUBPLOT_AREA_M2 <- (SUBPLOT_W_CM / 100) * (SUBPLOT_L_CM / 100)   # 0.33 m²
MIN_BB_COVER    <- 0.25   # •/??/... codes → lowest Braun–Blanquet value

bb_to_cover <- function(code, value) {
  coalesce(
    value,
    case_when(
      code %in% c("...", "??", "\u2022") ~ 0.25,
      code == "+"  ~ 0.75,
      code == "1"  ~ 3,
      code == "2"  ~ 15,
      code == "3"  ~ 37.5,
      code == "4"  ~ 62.5,
      code == "5"  ~ 87.5,
      TRUE         ~ NA_real_
    )
  )
}

build_diversity_plot <- function(dat_all) {
  # dat_all: plant records filtered to subplot_number %in%
  #          c(valid_subplots_main, outside_subplot_codes) or NA,
  #          with columns: plot_number, year, taxon_name, subplot_number,
  #                        scale_value_code, scale_value_value.
  
  # --- n_subplots per plot-year (6 or 8 depending on plot design) ---
  n_subplots_py <- dat_all %>%
    filter(subplot_number %in% valid_subplots_main) %>%
    distinct(plot_number, year, subplot_number) %>%
    count(plot_number, year, name = "n_subplots")
  
  # --- Subplot species: measured BB cover ---
  subplot_covers <- dat_all %>%
    filter(subplot_number %in% valid_subplots_main) %>%
    mutate(cover_num = bb_to_cover(scale_value_code, scale_value_value)) %>%
    filter(!is.na(cover_num), cover_num > 0) %>%
    group_by(plot_number, year, taxon_name) %>%
    summarise(cover = mean(cover_num, na.rm = TRUE), .groups = "drop")
  
  # --- Outside-subplot species: present in plot but absent from subplots ---
  # Only species not already captured by subplot_covers receive MIN_BB_COVER.
  outside_covers <- dat_all %>%
    filter(subplot_number %in% outside_subplot_codes) %>%
    filter(!is.na(taxon_name)) %>%
    distinct(plot_number, year, taxon_name) %>%
    anti_join(subplot_covers, by = c("plot_number", "year", "taxon_name")) %>%
    mutate(cover = MIN_BB_COVER)
  
  # --- Combined cover table → diversity indices ---
  all_covers <- bind_rows(subplot_covers, outside_covers)
  
  all_wide <- all_covers %>%
    tidyr::pivot_wider(
      names_from  = taxon_name,
      values_from = cover,
      values_fill = 0
    )
  
  mats <- all_wide %>%
    select(-plot_number, -year) %>%
    as.matrix()
  
  tibble(
    plot_number = all_wide$plot_number,
    year        = all_wide$year,
    S           = vegan::specnumber(mats),
    Shannon     = vegan::diversity(mats, index = "shannon"),
    Simpson     = vegan::diversity(mats, index = "simpson")
  ) %>%
    mutate(Pielou_J = Shannon / log(pmax(S, 1))) %>%
    left_join(n_subplots_py, by = c("plot_number", "year"))
}

# Input: all records with subplot_number in subplots OR outside codes (or NA)
v_plants_div_input <- v_plants_raw_all %>%
  filter(
    is.na(subplot_number) |
      subplot_number %in% c(valid_subplots_main, outside_subplot_codes)
  ) %>%
  mutate(
    field_date_parsed = parse_date_time(
      field_date,
      orders = c("Y-m-d", "d-m-Y", "d/m/Y", "Y/m/d", "d.m.Y", "Y.m.d"),
      quiet = TRUE
    ),
    year = year(field_date_parsed)
  ) %>%
  filter(!is.na(plot_number), !is.na(taxon_name), !is.na(year))

diversity_py        <- build_diversity_plot(v_plants_div_input)
diversity_py_hybrid <- diversity_py   # backward-compatible alias

revisited_yrs <- plot_years %>%
  select(plot_number, year_first, year_last) %>%
  filter(year_last > year_first)  # must have at least 1 year between surveys

pair_div <- function(div_tbl) {
  paired <- div_tbl %>%
    inner_join(revisited_yrs, by = "plot_number") %>%
    filter(year %in% c(year_first, year_last))
  
  first_div <- paired %>%
    filter(year == year_first) %>%
    select(
      plot_number,
      S_first        = S,
      Shannon_first  = Shannon,
      Simpson_first  = Simpson,
      Pielou_J_first = Pielou_J
    )
  
  last_div <- paired %>%
    filter(year == year_last) %>%
    select(
      plot_number,
      S_last        = S,
      Shannon_last  = Shannon,
      Simpson_last  = Simpson,
      Pielou_J_last = Pielou_J
    )
  
  first_div %>%
    full_join(last_div, by = "plot_number") %>%
    mutate(
      delta_S        = S_last        - S_first,
      delta_Shannon  = Shannon_last  - Shannon_first,
      delta_Simpson  = Simpson_last  - Simpson_first,
      delta_Pielou_J = Pielou_J_last - Pielou_J_first
    )
}

change_core <- pair_div(diversity_py_hybrid)

# S_total_py: total plot richness (subplots + outside) per plot-year.
# All indices in diversity_py already include outside-subplot species,
# so S_total equals S here. Kept for ordination-arrow backward compatibility.
S_total_py <- diversity_py %>% select(plot_number, year, S_total = S)

# ------------------------------------------------------------------------------
# 3) Variable roles and strict env set
# ------------------------------------------------------------------------------
id_vars <- c(
  "id", "plot_id", "subplot_id", "group_type_subplot_id",
  "project_id", "project_number", "analyser_id", "tag_id"
)
time_pairing_vars <- c("field_date", "field_date_parsed", "year", "plot_number")
species_cover_core_vars <- c("taxon_id", "taxon_name", "scale_value_code", "scale_value_value", "tag_name")
spatial_vars <- c("longitude_1", "latitude_1", "longitude_2", "latitude_2")
design_method_vars <- c(
  "plot_size_id", "plot_size_width", "plot_size_length",
  "subplot_size_id", "subplot_size_width", "subplot_size_length",
  "subplot_group_scale_id", "group_scale_name", "subplot_taxon_scale_id",
  "taxon_scale_name", "project_name", "analyser_name", "area"
)

var_catalog <- tibble(variable = names(v_plants_use)) %>%
  mutate(
    type = map_chr(variable, ~ class(v_plants_use[[.x]])[1]),
    prop_missing = map_dbl(variable, ~ mean(is.na(v_plants_use[[.x]]))),
    n_unique = map_dbl(variable, ~ n_distinct(v_plants_use[[.x]], na.rm = TRUE)),
    analysis_role = case_when(
      variable %in% id_vars ~ "identifier",
      variable %in% time_pairing_vars ~ "time_or_pairing",
      variable %in% species_cover_core_vars ~ "species_cover_core",
      variable %in% spatial_vars ~ "spatial",
      variable %in% design_method_vars ~ "design_or_metadata",
      type %in% c("numeric", "integer", "double") & prop_missing <= 0.40 ~ "env_numeric_candidate",
      type %in% c("character", "factor") & n_unique <= 20 & prop_missing <= 0.40 ~ "env_factor_candidate",
      type %in% c("logical") & prop_missing <= 0.40 ~ "env_factor_candidate",
      type %in% c("numeric", "integer", "double") & prop_missing > 0.40 ~ "deferred_sparse_numeric",
      type %in% c("character", "factor", "logical") & prop_missing > 0.40 ~ "deferred_sparse_categorical",
      TRUE ~ "review_manually"
    )
  )

vars_env_numeric <- var_catalog %>%
  filter(
    variable %in% vars_env_from_sheet,
    type %in% c("numeric", "integer", "double")
  ) %>%
  pull(variable)

vars_env_factor <- var_catalog %>%
  filter(
    variable %in% vars_env_from_sheet,
    type %in% c("character", "factor", "logical")
  ) %>%
  pull(variable)

mode_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

# For numeric env vars: treat NA as 0 when the variable was measured (has any non-NA).
# Columns that are entirely NA stay NA (variable never measured).
env_candidates_ps <- paired_long %>%
  mutate(across(
    any_of(vars_env_numeric),
    ~ if (all(is.na(.))) . else replace_na(., 0)
  )) %>%
  group_by(plot_number, sample) %>%
  summarise(
    across(any_of(vars_env_numeric), ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)),
    across(any_of(vars_env_factor), ~ mode_value(as.character(.x))),
    .groups = "drop"
  ) %>%
  left_join(gap_tbl, by = "plot_number")

# Elevation 
if (!is.null(my_table) && "plot_number" %in% names(my_table)) {
  env_candidates_ps <- env_candidates_ps %>%
    left_join(my_table, by = "plot_number")
}

vars_driver_strict_final <- c(
  "haeddypimetrar",
  "soil_depth",
  #"gap_years",
  "soil_type_name",
  "moisture_name",
  "topography_name",
  "permafrost"
)
#vars_driver_strict_final <- intersect(vars_driver_strict_final, names(env_candidates_ps))

# Add extra exploratory vars here (do NOT add to vars_driver_strict_final)
extra_env_exploratory <- c("habitat_type_name", "gap_years")

env_strict <- env_candidates_ps %>%
  select(
    plot_number, sample,
    all_of(vars_driver_strict_final),
    any_of(extra_env_exploratory)   # included but not used for drop_na()
  ) %>%
  distinct()


env_strict_cc <- env_strict %>%
  drop_na(all_of(vars_driver_strict_final))

plots_with_complete_pairs <- env_strict_cc %>%
  count(plot_number) %>%
  filter(n == 2) %>%
  pull(plot_number)

community_wide_strict <- community_wide %>%
  filter(plot_number %in% plots_with_complete_pairs)

meta_strict <- community_wide_strict %>%
  select(plot_number, sample) %>%
  left_join(env_strict_cc, by = c("plot_number", "sample")) %>%
  left_join(gap_tbl, by = "plot_number")

# Add area (county) to meta_strict
area_lookup <- v_plants_use %>%
  group_by(plot_number) %>%
  summarise(area = first(na.omit(area)), .groups = "drop")

meta_strict <- meta_strict %>%
  left_join(area_lookup, by = "plot_number")

mat_strict <- community_wide_strict %>%
  mutate(rn = paste(plot_number, sample, sep = "_")) %>%
  select(-plot_number, -sample) %>%
  tibble::column_to_rownames("rn") %>%
  as.matrix()
mat_strict <- mat_strict[, colSums(mat_strict) > 0]
# ------------------------------------------------------------------------------
# 4) Strict model objects used by downstream reports
# ------------------------------------------------------------------------------
#nmds_strict <- vegan::metaMDS(mat_strict, distance = "bray", trymax = 100, trace = FALSE)

nmds_cache <- file.path(base_dir, "outputs", "nmds_strict_cache.rds")
if (file.exists(nmds_cache)) {
  nmds_strict <- readRDS(nmds_cache)
} else {
  nmds_strict <- vegan::metaMDS(mat_strict, distance = "bray", trymax = 100, trace = FALSE)
  saveRDS(nmds_strict, nmds_cache)
}

envfit_formula_df <- meta_strict %>%
  dplyr::select(all_of(vars_driver_strict_final)) %>%
  mutate(
    across(where(is.character), as.factor),
    across(where(is.logical), as.factor)
  )

fit_strict <- vegan::envfit(nmds_strict, envfit_formula_df, permutations = 999)

perm_strict <- vegan::adonis2(
  mat_strict ~ sample,
  data = meta_strict,
  method = "bray",
  permutations = 999,
  strata = meta_strict$plot_number
)

bd_strict <- vegan::betadisper(vegan::vegdist(mat_strict, method = "bray"), group = meta_strict$sample)
bd_test_strict <- vegan::permutest(bd_strict, permutations = 999)


# ------------------------------------------------------------------------------
# 4b) Inland core matrix + NMDS + envfit (coastal habitat types excluded)
# ------------------------------------------------------------------------------
coastal_types <- c("Sandstrandarvist", "Strandmelhólavist",
                   "Sjávarfitjungsvist", "Sjávarkletta- og eyjavist")

# Detect outlier plots (same logic as Rmd)
nmds_scores_tmp <- as.data.frame(vegan::scores(nmds_strict, display = "sites")) %>%
  bind_cols(meta_strict %>% select(plot_number, sample))
centroid_tmp <- colMeans(nmds_scores_tmp[, c("NMDS1", "NMDS2")])
outlier_plots <- nmds_scores_tmp %>%
  mutate(dist_c = sqrt((NMDS1 - centroid_tmp[1])^2 + (NMDS2 - centroid_tmp[2])^2)) %>%
  arrange(desc(dist_c)) %>%
  distinct(plot_number, .keep_all = TRUE) %>%
  slice_head(n = 3) %>%
  pull(plot_number)

# Exclude outlier plots AND coastal habitat type plots
coastal_plots <- meta_strict %>%
  filter(habitat_type_name %in% coastal_types) %>%
  distinct(plot_number) %>%
  pull(plot_number)

exclude_inland <- union(outlier_plots, coastal_plots)
keep_inland    <- !meta_strict$plot_number %in% exclude_inland

meta_core_inland <- meta_strict[keep_inland, ]
mat_core_inland  <- mat_strict[keep_inland, ]
mat_core_inland  <- mat_core_inland[, colSums(mat_core_inland) > 0]

# NMDS — cache separately from nmds_strict
nmds_inland_cache <- file.path(base_dir, "outputs", "nmds_inland_cache.rds")
if (file.exists(nmds_inland_cache)) {
  nmds_core_inland <- readRDS(nmds_inland_cache)
} else {
  nmds_core_inland <- vegan::metaMDS(
    mat_core_inland,
    distance = "bray",
    k        = 2,
    trymax   = 200,
    trace    = FALSE
  )
  saveRDS(nmds_core_inland, nmds_inland_cache)
}

# envfit on inland matrix — full driver set + habitat + area
env_core_inland <- meta_core_inland %>%
  select(all_of(vars_driver_strict_final),
         any_of(c("habitat_type_name", "area"))) %>%
  mutate(
    across(where(is.character), as.factor),
    across(where(is.logical),   as.factor)
  )

fit_core_inland <- vegan::envfit(
  nmds_core_inland,
  env_core_inland,
  permutations = 999,
  na.rm        = TRUE
)

# Keep this explicit for downstream reporting consistency
driver_final <- c(
  "haplontuthekja",
  "soil_depth",
  "total_cover",
  "vegetation_height_mean",
  "soil_type_name",
  "moisture_name"
)


"area" %in% names(meta_strict)
"area" %in% names(meta_core_inland)
exists("fit_core_inland")
nrow(meta_core_inland)

source("v_plants_prep.R")

saveRDS(
  list(
    community_long_paired = community_long_paired,
    meta_strict           = meta_strict,
    meta_core_inland      = meta_core_inland,
    mat_core_inland       = mat_core_inland,
    nmds_core_inland      = nmds_core_inland,
    fit_core_inland       = fit_core_inland,
    outlier_plots         = outlier_plots
  ),
  file = "v_plants_app_data.rds"
)
