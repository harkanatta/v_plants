# v_plants — Iceland vegetation resurvey analysis (R)

This repository contains an R workflow for analysing long-interval paired vegetation resurveys in Icelandic plots.

## What this project does

- Reads the exported vegetation dataset from `plants_export/v_plants_allt.csv`
- Applies two *explicit* filter sheets:
  - `v_plants_taxa_exclude_flags.csv`: taxa to exclude globally (flagged rows)
  - `v_plants_envvar_selection.csv`: environmental variables selected for analysis (flagged rows)
- Builds a paired **first vs last** dataset per plot
- Harmonises cover values (Braun–Blanquet codes + numeric estimates → one numeric cover scale)
- Produces analysis-ready objects (community matrices, NMDS, PERMANOVA, envfit, etc.)

## Repository structure (key files)

- `v_plants_prep.R`: single “source of truth” for data prep + core model objects
- `v_plants_variable_screening.Rmd`: diagnostics / transparency around variables and strict datasets
- `v_plants_analysis.Rmd`: main analysis report (sources `v_plants_prep.R`)
- `v_plants_next_steps.R`: follow-up analyses and figures (assumes prep objects exist)
- `fyrirlestur/`: ioslides presentations + CSS

## Data and Git LFS

`plants_export/v_plants_allt.csv` is larger than GitHub’s normal file-size limit. This repo is configured to track it with **Git LFS**.

After cloning, run:

```bash
git lfs install
git lfs pull
```

## Installing R packages (RStudio / container library path)

In this environment, install packages into:

```r
install.packages("vegan", lib = "/home/rstudio/R/library")
```

To make that library visible in a session:

```r
.libPaths(c("/home/rstudio/R/library", .libPaths()))
```

## How to run

In R:

```r
source("v_plants_prep.R")
rmarkdown::render("v_plants_variable_screening.Rmd")
rmarkdown::render("v_plants_analysis.Rmd")
rmarkdown::render("fyrirlestur/v_plants_presentation.Rmd")
rmarkdown::render("fyrirlestur/v_plants_presentation_is.Rmd")
```

## Notes on encoding (þæöð)

When exporting CSVs intended for Excel, prefer:

```r
readr::write_excel_csv(df, "outputs/file_excel.csv")
```

