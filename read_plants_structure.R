# Read and inspect structure of plants_export/v_plants_allt.csv
library(tidyverse)

path <- file.path(getwd(), "plants_export", "v_plants_allt.csv")
if (!file.exists(path)) stop("File not found: ", path)

v_plants <- read_csv(path, locale = locale(encoding = "UTF-8"), show_col_types = FALSE)

# Structure
str(v_plants)

# Column names (numbered)
names(v_plants) %>% enframe(value = "column") %>% print(n = Inf)

# Compact summary
glimpse(v_plants)

# Key columns for diversity / joining: plot identifier, taxon, cover scale
# plot_number = site/plot code (e.g. 2016-A22-01, NV-sk-01, KMM68)
# taxon_name, scale_value_code (Braun-Blanquet: • + 1 2 3 4), scale_value_value (numeric)
v_plants %>%
  select(plot_number, area, field_date, subplot_id, subplot_number, total_cover,
         habitat_type_name, taxon_name, scale_value_code, scale_value_value, tag_name) %>%
  head(20)

# Distinct plot_number (sites) and scale types
v_plants %>% distinct(plot_number) %>% arrange(plot_number) %>% print(n = 50)
v_plants %>% count(taxon_scale_name, scale_value_code, name = "n") %>% print(n = 30)
