# -------------------------------------------------------------------------------
# v_plants_glimpse_for_llm.R — Run in R console; copy console output to feed an LLM
# Assumes: setwd() is project root (e.g. d:/gagnagrunnurNI or where plants_export lives)
# -------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

base_candidates <- c(getwd(), normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = FALSE))
plants_path <- NULL
for (b in base_candidates) {
  p <- file.path(b, "plants_export", "v_plants_allt.csv")
  if (file.exists(p)) { plants_path <- p; break }
}
if (is.null(plants_path)) stop("Could not find plants_export/v_plants_allt.csv")
base_dir <- dirname(dirname(plants_path))

cat("\n========== 1) INPUT CSVs ==========\n\n")

# Main data
cat("--- plants_export/v_plants_allt.csv (first as v_plants_raw) ---\n")
v_plants_raw <- readr::read_csv(plants_path, show_col_types = FALSE)
dplyr::glimpse(v_plants_raw)
cat("Rows:", nrow(v_plants_raw), " Cols:", ncol(v_plants_raw), "\n\n")

# Taxa exclude
taxon_path <- file.path(base_dir, "v_plants_taxa_exclude_flags.csv")
if (file.exists(taxon_path)) {
  cat("--- v_plants_taxa_exclude_flags.csv ---\n")
  taxon_exclude_tbl <- readr::read_csv(taxon_path, col_names = c("exclude_flag", "taxon_name"), show_col_types = FALSE)
  dplyr::glimpse(taxon_exclude_tbl)
  cat("Rows:", nrow(taxon_exclude_tbl), "\n\n")
} else {
  cat("--- v_plants_taxa_exclude_flags.csv not found ---\n\n")
}

# Env var selection
envvar_path <- file.path(base_dir, "v_plants_envvar_selection.csv")
if (file.exists(envvar_path)) {
  cat("--- v_plants_envvar_selection.csv ---\n")
  raw_cols_tbl <- readr::read_csv(envvar_path, col_names = c("flag", "variable"), show_col_types = FALSE)
  dplyr::glimpse(raw_cols_tbl)
  cat("Rows:", nrow(raw_cols_tbl), "\n\n")
} else {
  cat("--- v_plants_envvar_selection.csv not found ---\n\n")
}

cat("\n========== 2) SOURCE PREP (creates all outputs) ==========\n\n")
source(file.path(base_dir, "v_plants_prep.R"), echo = FALSE)

cat("\n========== 3) OUTPUTS FROM v_plants_prep.R ==========\n\n")

cat("--- v_plants2 (main filtered long data) ---\n")
dplyr::glimpse(v_plants2)
cat("Rows:", nrow(v_plants2), "\n\n")

cat("--- v_plants_use (plots in heiti_include) ---\n")
dplyr::glimpse(v_plants_use)
cat("Rows:", nrow(v_plants_use), "\n\n")

cat("--- plot_years ---\n")
dplyr::glimpse(plot_years)
cat("Rows:", nrow(plot_years), "\n\n")

cat("--- community_wide (first/last wide) ---\n")
dplyr::glimpse(community_wide)
cat("Rows:", nrow(community_wide), " Cols:", ncol(community_wide), "\n\n")

cat("--- meta, mat (dim only) ---\n")
cat("meta:", nrow(meta), "x", ncol(meta), "  mat:", nrow(mat), "x", ncol(mat), "\n\n")

cat("--- bc_per_plot ---\n")
dplyr::glimpse(bc_per_plot)
cat("Rows:", nrow(bc_per_plot), "\n\n")

cat("--- gap_tbl ---\n")
dplyr::glimpse(gap_tbl)
cat("\n")

cat("--- bc_with_gap ---\n")
dplyr::glimpse(bc_with_gap)
cat("Rows:", nrow(bc_with_gap), "\n\n")

if (!is.null(bc_enriched) && nrow(bc_enriched) > 0) {
  cat("--- bc_enriched ---\n")
  dplyr::glimpse(bc_enriched)
  cat("Rows:", nrow(bc_enriched), "\n\n")
} else {
  cat("--- bc_enriched: NULL or empty (sf/gpkg not used) ---\n\n")
}

cat("--- diversity_py_core ---\n")
dplyr::glimpse(diversity_py_core)
cat("Rows:", nrow(diversity_py_core), "\n\n")

cat("--- diversity_py_hybrid (S from core+outside, H/Simpson/J from core) ---\n")
dplyr::glimpse(diversity_py_hybrid)
cat("Rows:", nrow(diversity_py_hybrid), "\n\n")

cat("--- change_core (paired diversity change) ---\n")
dplyr::glimpse(change_core)
cat("Rows:", nrow(change_core), "\n\n")

cat("--- var_catalog (first 20 rows) ---\n")
print(head(var_catalog, 20))
cat("Total variables:", nrow(var_catalog), "\n\n")

cat("--- env_candidates_ps ---\n")
dplyr::glimpse(env_candidates_ps)
cat("Rows:", nrow(env_candidates_ps), "\n\n")

cat("--- env_strict ---\n")
dplyr::glimpse(env_strict)
cat("Rows:", nrow(env_strict), "\n\n")

cat("--- community_wide_strict, meta_strict, mat_strict (dims) ---\n")
cat("community_wide_strict:", nrow(community_wide_strict), "x", ncol(community_wide_strict), "\n")
cat("meta_strict:", nrow(meta_strict), "x", ncol(meta_strict), "\n")
cat("mat_strict:", nrow(mat_strict), "x", ncol(mat_strict), "\n\n")

cat("--- vars_driver_strict_final ---\n")
print(vars_driver_strict_final)
cat("\n--- driver_final ---\n")
print(driver_final)
cat("\n")

cat("--- perm_strict (adonis2) ---\n")
print(perm_strict)
cat("\n")

cat("--- fit_strict (envfit) vectors ---\n")
if (!is.null(fit_strict$vectors)) print(fit_strict$vectors)
cat("\n")

cat("--- nmds_strict (metaMDS) ---\n")
cat("converged:", nmds_strict$convergence, " stress:", nmds_strict$stress, "\n")
cat("points dim:", nrow(nmds_strict$points), "x", ncol(nmds_strict$points), "\n\n")

cat("--- bd_strict (betadisper) ---\n")
cat("groups:", paste(names(bd_strict$group), collapse = ", "), "\n")
cat("distances length:", length(bd_strict$distances), "\n")
cat("--- bd_test_strict (permutest) ---\n")
print(bd_test_strict)
cat("\n")

cat("========== END GLIMPSE ==========\n")
