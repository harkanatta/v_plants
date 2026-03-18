# One-off script: fix Icelandic column names in v_plants_allt.csv
# Run from project root. Deletes self after use.

library(readr)

csv_path <- "plants_export/v_plants_allt.csv"
if (!file.exists(csv_path)) stop("CSV not found: ", csv_path)

# Old (Icelandic) -> New (ASCII)
name_map <- c(
  "háplöntuþekja" = "haplontuthekja",
  "mosaþekja" = "mosathekja",
  "fléttuþekja" = "flettuthekja",
  "engjaskófir" = "engjaskofir",
  "breiskjufléttur" = "breiskjuflettur",
  "Mosa- og fléttuskán" = "mosa_og_flettuskan",
  "kræðufléttur" = "kraeduflettur",
  "svarðmosi" = "svardmosi",
  "grýtniþekja" = "grytnithekja",
  "sinuþekja" = "sinuthekja",
  "greinaþekja" = "greinathekja",
  "barrþekja" = "barrthekja",
  "heildarþekja" = "heildarthekja"
)

dat <- read_csv(csv_path, show_col_types = FALSE, locale = locale(encoding = "UTF-8"))
nms <- names(dat)
for (i in seq_along(name_map)) {
  idx <- which(nms == names(name_map)[i])
  if (length(idx) > 0) nms[idx] <- name_map[i]
}
names(dat) <- nms
write_csv(dat, csv_path)
message("Updated ", csv_path, " with ASCII column names")
