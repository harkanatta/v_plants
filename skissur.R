dca_vars <- v_plants_use %>%
  select(
    plot_number,
    subplot_number,
    year,
    vegetation_height_mean,
    mosaþekja,
    total_cover,
    háplöntuþekja,
    fléttuþekja,
    `Mosa- og fléttuskán`,
    grýtniþekja,
    topography_name
  )

dca_vars
