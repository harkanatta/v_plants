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
  "temperature_july", "temperature_year", "heildarþekja", "subplot_carbon",
  "subplot_nitrogen", "subplot_ph", "plant_community_other",
  "plant_community_3_code", "plant_community_3_id", "plant_community_3_name",
  "slope_aspect", "slope", "autt", "surface_type_code", "surface_type_id",
  "surface_type_name", "soil_temperature_2", "plant_community_cover"
)
vars_exclude_text <- c(
  "comments", "subplot_comments", "images", "subplot_images", "project_descrtiption"
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
    }
  ) %>%
  filter(!is.na(plot_number), !is.na(taxon_name), !is.na(year))

v_plants_use <- v_plants2 %>%
  select(-any_of(vars_exclude_strict))

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
# 2a) Per-plot Bray–Curtis distance (first vs last) + gap in years
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
# 2b) Plot-level diversity with and without outside-subplot records (optional)
# ------------------------------------------------------------------------------
build_diversity_py <- function(dat) {
  dat %>%
    mutate(
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
    filter(!is.na(cover_num), cover_num > 0) %>%
    group_by(plot_number, year, taxon_name) %>%
    summarise(cover = mean(cover_num, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = taxon_name, values_from = cover, values_fill = 0) %>%
    {
      mats <- as.matrix(dplyr::select(., -plot_number, -year))
      tibble(
        plot_number = .$plot_number,
        year        = .$year,
        S           = vegan::specnumber(mats),
        Shannon     = vegan::diversity(mats, index = "shannon"),
        Simpson     = vegan::diversity(mats, index = "simpson")
      ) %>%
        mutate(Evenness = Shannon / log(pmax(S, 1)))
    }
}

v_plants_div_core <- v_plants_raw_all %>%
  filter(is.na(subplot_number) | subplot_number %in% valid_subplots_main) %>%
  mutate(
    field_date_parsed = parse_date_time(
      field_date,
      orders = c("Y-m-d", "d-m-Y", "d/m/Y", "Y/m/d", "d.m.Y", "Y.m.d"),
      quiet = TRUE
    ),
    year = year(field_date_parsed)
  ) %>%
  filter(!is.na(plot_number), !is.na(taxon_name), !is.na(year))

v_plants_div_all <- v_plants_raw_all %>%
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

diversity_py_core <- build_diversity_py(v_plants_div_core)
diversity_py_all  <- build_diversity_py(v_plants_div_all)

revisited_yrs <- plot_years %>% select(plot_number, year_first, year_last)

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
      Evenness_first = Evenness
    )

  last_div <- paired %>%
    filter(year == year_last) %>%
    select(
      plot_number,
      S_last        = S,
      Shannon_last  = Shannon,
      Simpson_last  = Simpson,
      Evenness_last = Evenness
    )

  first_div %>%
    full_join(last_div, by = "plot_number") %>%
    mutate(
      delta_S        = S_last        - S_first,
      delta_Shannon  = Shannon_last  - Shannon_first,
      delta_Simpson  = Simpson_last  - Simpson_first,
      delta_Evenness = Evenness_last - Evenness_first
    )
}

change_core <- pair_div(diversity_py_core)
change_all  <- pair_div(diversity_py_all)

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

env_candidates_ps <- paired_long %>%
  group_by(plot_number, sample) %>%
  summarise(
    across(all_of(vars_env_numeric), ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)),
    across(all_of(vars_env_factor), ~ mode_value(as.character(.x))),
    .groups = "drop"
  ) %>%
  left_join(gap_tbl, by = "plot_number")

# Optional: elevation from gpkg (used for bc_enriched and for meta_strict when available)
elev_tbl <- NULL
if (requireNamespace("sf", quietly = TRUE)) {
  gpkg_candidates <- c(
    file.path(base_dir, "outputs", "qgis_share", "snid_with_elevation.gpkg"),
    file.path(base_dir, "gagnagrunnurNI", "outputs", "qgis_share", "snid_with_elevation.gpkg")
  )
  gpkg_path <- gpkg_candidates[file.exists(gpkg_candidates)][1]
  if (!is.na(gpkg_path)) {
    snid_elev <- sf::st_read(gpkg_path, quiet = TRUE)
    elev_tbl <- snid_elev %>%
      sf::st_drop_geometry() %>%
      group_by(plot_number) %>%
      summarise(haeddypimetrar = first(haeddypimetrar), .groups = "drop")
    env_first <- env_candidates_ps %>%
      filter(sample == "first") %>%
      select(plot_number, moisture_name, habitat_type_name, topography_name, permafrost) %>%
      distinct()
    bc_enriched <- bc_with_gap %>%
      left_join(env_first, by = "plot_number") %>%
      left_join(elev_tbl, by = "plot_number")
  } else {
    bc_enriched <- NULL
  }
} else {
  bc_enriched <- NULL
}
if (!is.null(elev_tbl)) {
  env_candidates_ps <- env_candidates_ps %>%
    left_join(elev_tbl, by = "plot_number")
}

vars_driver_strict_final <- c(
  "háplöntuþekja",
  "soil_depth",
  "total_cover",
  "vegetation_height_mean",
  "soil_type_name",
  "moisture_name",
  "topography_name",
  "permafrost"
)
vars_driver_strict_final <- intersect(vars_driver_strict_final, names(env_candidates_ps))

env_strict <- env_candidates_ps %>%
  select(plot_number, sample, all_of(vars_driver_strict_final)) %>%
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
if (!is.null(elev_tbl)) {
  meta_strict <- meta_strict %>%
    left_join(elev_tbl, by = "plot_number")
}

mat_strict <- community_wide_strict %>%
  select(-plot_number, -sample) %>%
  as.matrix()

# ------------------------------------------------------------------------------
# 4) Strict model objects used by downstream reports
# ------------------------------------------------------------------------------
nmds_strict <- vegan::metaMDS(mat_strict, distance = "bray", trymax = 100, trace = FALSE)

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

# Keep this explicit for downstream reporting consistency
driver_final <- c(
  "háplöntuþekja",
  "soil_depth",
  "total_cover",
  "vegetation_height_mean",
  "soil_type_name",
  "moisture_name"
)

