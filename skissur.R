# Fj. hápl. snið: distinct vascular plant species per plot×sample
richness_snid <- v_plants_use %>%
  filter(!is.na(subplot_number), subplot_number %in% valid_subplots_main) %>%
  inner_join(
    plot_years %>% select(plot_number, year_first, year_last),
    by = "plot_number"
  ) %>%
  mutate(sample = case_when(
    year == year_first ~ "first",
    year == year_last  ~ "last",
    TRUE               ~ NA_character_
  )) %>%
  filter(!is.na(sample)) %>%
  group_by(plot_number, sample) %>%
  summarise(fj_hapl_snid = n_distinct(taxon_name), .groups = "drop")

# Fj. hápl. í reit: distinct vascular plant species per subplot, then mean across subplots
richness_reit <- v_plants_use %>%
  filter(!is.na(subplot_number), subplot_number %in% valid_subplots_main) %>%
  inner_join(
    plot_years %>% select(plot_number, year_first, year_last),
    by = "plot_number"
  ) %>%
  mutate(sample = case_when(
    year == year_first ~ "first",
    year == year_last  ~ "last",
    TRUE               ~ NA_character_
  )) %>%
  filter(!is.na(sample)) %>%
  group_by(plot_number, sample, subplot_number) %>%
  summarise(n = n_distinct(taxon_name), .groups = "drop") %>%
  group_by(plot_number, sample) %>%
  summarise(fj_hapl_reit = mean(n), .groups = "drop")
