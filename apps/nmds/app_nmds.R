# =============================================================================
# v_plants ??? NMDS Interactive Biplot  (v4)
# Single-file Shiny app
#
# Dependencies: shiny, tidyverse, vegan, bslib, plotly
# Data:         sources v_plants_prep.R from the working directory
#
# Run with: shiny::runApp("app_nmds.R")
# =============================================================================

library(shiny)
library(tidyverse)
library(vegan)
library(bslib)
library(plotly)

# -----------------------------------------------------------------------------
# 0.  Load data
# -----------------------------------------------------------------------------
if (!exists("community_long_paired")) {
  dat <- readRDS("v_plants_app_data.rds")
  list2env(dat, envir = .GlobalEnv)
}

# -----------------------------------------------------------------------------
# 1.  Static objects
# -----------------------------------------------------------------------------

# Site scores ??? explicit column names, categorical + numeric all in one object
site_scores <- as.data.frame(vegan::scores(nmds_core_inland, display = "sites")) %>%
  bind_cols(meta_core_inland %>% select(
    plot_number, sample,
    habitat_type_name, soil_type_name, moisture_name,
    topography_name, area, permafrost,
    soil_depth, haeddypimetrar,
    any_of("total_cover")
  ))

# Hover numeric ??? ONLY the first???last pivoted numeric columns
# Categoricals are already in site_scores; do NOT include them here
# to avoid duplicate column names after left_join in renderPlotly
hover_num_cols <- intersect(c("soil_depth", "haeddypimetrar", "total_cover"),
                            names(meta_core_inland))

hover_meta <- meta_core_inland %>%
  select(plot_number, sample, all_of(hover_num_cols)) %>%
  pivot_wider(
    id_cols     = plot_number,
    names_from  = sample,
    values_from = all_of(hover_num_cols),
    names_sep   = "_"
  )

# Top 5 species per plot (last survey)
top5 <- community_long_paired %>%
  filter(plot_number %in% unique(site_scores$plot_number),
         sample == "last") %>%
  group_by(plot_number) %>%
  arrange(desc(cover)) %>%
  slice_head(n = 5) %>%
  summarise(top5 = paste(taxon_name, collapse = "<br>&nbsp;&nbsp;"),
            .groups = "drop")

# Trajectory pairs (first ??? last, one row per plot)
traj <- site_scores %>%
  select(plot_number, sample, NMDS1, NMDS2) %>%
  pivot_wider(names_from  = sample,
              values_from = c(NMDS1, NMDS2),
              names_sep   = "_") %>%
  filter(complete.cases(.))

# envfit: continuous vectors
ordi_radius <- max(abs(c(site_scores$NMDS1, site_scores$NMDS2))) * 0.65

envfit_vectors <- tryCatch({
  as.data.frame(fit_core_inland$vectors$arrows) %>%
    rownames_to_column("variable") %>%
    rename(ax1 = NMDS1, ax2 = NMDS2) %>%
    mutate(
      r2      = as.numeric(fit_core_inland$vectors$r),
      p_value = as.numeric(fit_core_inland$vectors$pvals)
    ) %>%
    filter(p_value <= 0.05) %>%
    mutate(xend = ax1 * ordi_radius,
           yend = ax2 * ordi_radius)
}, error = function(e) tibble())

# envfit: factor centroids
envfit_factors <- tryCatch({
  as.data.frame(fit_core_inland$factors$centroids) %>%
    rownames_to_column("label") %>%
    rename(cx = NMDS1, cy = NMDS2) %>%
    mutate(
      p_value = sapply(label, function(nm) {
        matched <- names(fit_core_inland$factors$pvals)[
          startsWith(nm, names(fit_core_inland$factors$pvals))
        ]
        if (length(matched)) fit_core_inland$factors$pvals[matched[1]] else NA_real_
      })
    ) %>%
    filter(!is.na(p_value), p_value <= 0.05)
}, error = function(e) tibble())

# Colour-by and filter choices
colour_choices <- c(
  "Moisture"            = "moisture_name",
  "Soil type"           = "soil_type_name",
  "Topography"          = "topography_name",
  "Habitat type"        = "habitat_type_name",
  "Area / county"       = "area",
  "Sample (first/last)" = "sample"
  )
colour_choices <- colour_choices[colour_choices %in% names(site_scores)]

mk_filter <- function(col) {
  c("All" = "All", sort(unique(na.omit(site_scores[[col]]))))
}
area_choices     <- mk_filter("area")
habitat_choices  <- mk_filter("habitat_type_name")
soil_choices     <- mk_filter("soil_type_name")
moisture_choices <- mk_filter("moisture_name")

# =============================================================================
# 2.  UI
# =============================================================================

ui <- page_sidebar(
  title = NULL,
  theme = bs_theme(
    bootswatch   = "flatly",
    base_font    = font_google("IBM Plex Sans"),
    heading_font = font_google("IBM Plex Mono"),
    primary      = "#2c6e49",
    secondary    = "#52796f",
    "sidebar-bg" = "#f4f6f4"
  ),
  
  tags$head(tags$style(HTML("
    .sec-label {
      font-size:0.7rem; font-weight:600; color:#6c757d;
      text-transform:uppercase; letter-spacing:0.05em;
      margin:0.9rem 0 0.25rem;
    }
    .stress-bar {
      background:#f4f6f4; border:1px solid #dee2e6; border-radius:6px;
      padding:0.4rem 0.9rem; margin-bottom:0.6rem;
      font-size:0.78rem; color:#495057;
      font-family:'IBM Plex Mono',monospace;
    }
    .stress-bar span { color:#2c6e49; font-weight:700; }
  "))),
  
  sidebar = sidebar(
    width = 255,
    
    tags$div(
      style = "margin-bottom:0.4rem;",
      tags$span(
        style = "font-family:'IBM Plex Mono',monospace; color:#2c6e49;
                 font-size:0.85rem; font-weight:600;",
        "v_plants"
      ),
      tags$span(
        style = "font-size:0.72rem; color:#6c757d; margin-left:0.4rem;",
        "NMDS biplot"
      )
    ),
    
    tags$hr(style = "margin:0.3rem 0 0.6rem;"),
    
    div(class = "sec-label", "Colour by"),
    selectInput("colour_by", NULL,
                choices  = colour_choices,
                selected = colour_choices[1]),
    
    tags$hr(style = "margin:0.6rem 0 0.4rem;"),
    
    div(class = "sec-label", "Highlight"),
    selectInput("filter_area",     "Area / county", choices = area_choices,     selected = "All"),
    selectInput("filter_habitat",  "Habitat type",  choices = habitat_choices,  selected = "All"),
    selectInput("filter_soil",     "Soil type",     choices = soil_choices,     selected = "All"),
    selectInput("filter_moisture", "Moisture",      choices = moisture_choices, selected = "All"),
    
    tags$hr(style = "margin:0.6rem 0 0.4rem;"),
    
    div(class = "sec-label", "Overlays"),
    checkboxInput("show_vectors", "envfit vectors (continuous)", value = TRUE),
    checkboxInput("show_factors", "envfit centroids (factors)",  value = FALSE),
    checkboxInput("show_traj",    "First \u2192 last lines",      value = TRUE),
    checkboxInput("show_labels",  "Plot labels",                  value = FALSE)
  ),
  
  div(class = "stress-bar",
      HTML(sprintf(
        "NMDS stress: <span>%.4f</span> &nbsp;&middot;&nbsp; k&nbsp;=&nbsp;%d &nbsp;&middot;&nbsp; %d plots &nbsp;&middot;&nbsp; coastal types + outliers excluded",
        nmds_core_inland$stress,
        nmds_core_inland$ndim,
        length(unique(site_scores$plot_number))
      ))
  ),
  
  card(
    card_body(plotlyOutput("nmds_plot", height = "75vh")),
    full_screen = TRUE,
    fill = FALSE
  )
)

# =============================================================================
# 3.  Server
# =============================================================================

server <- function(input, output, session) {
  
  # ?????? Highlighted plots ??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
  highlighted_plots <- reactive({
    plots <- unique(site_scores$plot_number)
    
    if (input$filter_area != "All")
      plots <- intersect(plots,
                         site_scores %>% filter(area == input$filter_area) %>%
                           pull(plot_number) %>% unique())
    if (input$filter_habitat != "All")
      plots <- intersect(plots,
                         site_scores %>% filter(habitat_type_name == input$filter_habitat) %>%
                           pull(plot_number) %>% unique())
    if (input$filter_soil != "All")
      plots <- intersect(plots,
                         site_scores %>% filter(soil_type_name == input$filter_soil) %>%
                           pull(plot_number) %>% unique())
    if (input$filter_moisture != "All")
      plots <- intersect(plots,
                         site_scores %>% filter(moisture_name == input$filter_moisture) %>%
                           pull(plot_number) %>% unique())
    plots
  })
  
  # ?????? NMDS plot ??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
  output$nmds_plot <- renderPlotly({
    
    hp         <- highlighted_plots()
    cb         <- input$colour_by
    any_filter <- !all(c(input$filter_area, input$filter_habitat,
                         input$filter_soil, input$filter_moisture) == "All")
    
    # Join only hover_meta (numeric only) ??? no duplicate categorical columns
    df <- site_scores %>%
      left_join(top5,       by = "plot_number") %>%
      left_join(hover_meta, by = "plot_number") %>%
      mutate(
        highlighted = plot_number %in% hp,
        pt_opacity  = if (any_filter) if_else(highlighted, 0.90, 0.10) else 0.85,
        symbol_val  = if_else(sample == "first", "circle", "triangle-up")
      )
    
    # Colour column ??? outside mutate, base R extraction
    df$colour_val <- if (cb %in% names(df)) as.character(df[[cb]]) else "all"
    df$colour_val[is.na(df$colour_val)] <- "unknown"
    
    # Hover text ??? built outside mutate, column by column
    ht <- paste0("<b>", df$plot_number, "</b> [", df$sample, "]<br>")
    if ("habitat_type_name" %in% names(df))
      ht <- paste0(ht, "Habitat: ",    df$habitat_type_name, "<br>")
    if ("soil_type_name" %in% names(df))
      ht <- paste0(ht, "Soil: ",       df$soil_type_name,    "<br>")
    if ("moisture_name" %in% names(df))
      ht <- paste0(ht, "Moisture: ",   df$moisture_name,     "<br>")
    if ("permafrost" %in% names(df))
      ht <- paste0(ht, "Permafrost: ", df$permafrost,        "<br>")
    if ("area" %in% names(df))
      ht <- paste0(ht, "Area: ",       df$area,              "<br>")
    if (all(c("soil_depth_first", "soil_depth_last") %in% names(df)))
      ht <- paste0(ht, "Soil depth: ", round(df$soil_depth_first, 1),
                   " \u2192 ", round(df$soil_depth_last, 1), " cm<br>")
    if (all(c("haeddypimetrar_first", "haeddypimetrar_last") %in% names(df)))
      ht <- paste0(ht, "Elevation: ", round(df$haeddypimetrar_first, 0),
                   " \u2192 ", round(df$haeddypimetrar_last, 0), " m<br>")
    if ("top5" %in% names(df))
      ht <- paste0(ht, "<br><i>Top 5 (last):</i><br>&nbsp;&nbsp;", df$top5)
    df$hover_text <- ht
    
    p <- plot_ly() %>% config(displayModeBar = FALSE)
    
    # ?????? Trajectory lines ??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
    if (input$show_traj) {
      traj_df <- traj %>%
        mutate(
          highlighted  = plot_number %in% hp,
          line_col     = if_else(highlighted | !any_filter, "#52796f", "#cccccc"),
          line_opacity = if (any_filter) if_else(highlighted, 0.60, 0.07) else 0.30
        )
      
      for (i in seq_len(nrow(traj_df))) {
        r <- traj_df[i, ]
        p <- p %>% add_lines(
          x = c(r$NMDS1_first, r$NMDS1_last),
          y = c(r$NMDS2_first, r$NMDS2_last),
          line       = list(color = r$line_col, dash = "dot", width = 0.9),
          opacity    = r$line_opacity,
          showlegend = FALSE,
          hoverinfo  = "skip",
          inherit    = FALSE
        )
      }
    }
    
    # ?????? Points by colour group ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
    for (grp in sort(unique(df$colour_val))) {
      sub <- df %>% filter(colour_val == grp)
      p <- p %>% add_trace(
        data          = sub,
        x             = ~NMDS1,
        y             = ~NMDS2,
        type          = "scatter",
        mode          = "markers",
        name          = grp,
        marker        = list(
          symbol  = ~symbol_val,
          size    = 9,
          opacity = ~pt_opacity,
          line    = list(color = "white", width = 0.6)
        ),
        text          = ~hover_text,
        hovertemplate = "%{text}<extra></extra>",
        inherit       = FALSE
      )
    }
    
    # ?????? Plot labels ????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
    if (input$show_labels) {
      lbl_df <- df %>%
        group_by(plot_number) %>%
        summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2),
                  highlighted = any(highlighted), .groups = "drop") %>%
        mutate(lbl_opacity = if (any_filter) if_else(highlighted, 1, 0.12) else 0.55)
      
      p <- p %>% add_text(
        data       = lbl_df,
        x = ~NMDS1, y = ~NMDS2,
        text       = ~plot_number,
        textfont   = list(size = 7.5, family = "IBM Plex Mono",
                          color = "rgba(60,60,60,1)"),
        opacity    = ~lbl_opacity,
        showlegend = FALSE,
        hoverinfo  = "skip",
        inherit    = FALSE
      )
    }
    
    # ?????? envfit continuous vectors ??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
    if (input$show_vectors && nrow(envfit_vectors) > 0) {
      for (i in seq_len(nrow(envfit_vectors))) {
        r <- envfit_vectors[i, ]
        
        p <- p %>% add_annotations(
          x = r$xend, y = r$yend,
          ax = 0, ay = 0,
          xref = "x", yref = "y", axref = "x", ayref = "y",
          showarrow  = TRUE,
          arrowhead  = 2, arrowsize = 0.8,
          arrowwidth = 1.5, arrowcolor = "#c0392b",
          text       = "",
          inherit    = FALSE
        )
        
        p <- p %>% add_annotations(
          x         = r$xend * 1.05,
          y         = r$yend * 1.05,
          xref      = "x", yref = "y",
          showarrow = FALSE,
          text      = paste0("<b>", r$variable, "</b>  r2=",
                             round(r$r2, 2), " p=", round(r$p_value, 3)),
          font      = list(color = "#c0392b", size = 10,
                           family = "IBM Plex Mono"),
          xanchor   = if (r$xend >= 0) "left" else "right",
          yanchor   = if (r$yend >= 0) "bottom" else "top",
          inherit   = FALSE
        )
      }
    }
    
    # ?????? envfit factor centroids ????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
    if (input$show_factors && nrow(envfit_factors) > 0) {
      p <- p %>% add_trace(
        data          = envfit_factors,
        x             = ~cx,
        y             = ~cy,
        type          = "scatter",
        mode          = "markers+text",
        name          = "factor centroids",
        marker        = list(symbol = "diamond", size = 7,
                             color = "#8e44ad", opacity = 0.75,
                             line = list(color = "white", width = 0.5)),
        text          = ~label,
        textfont      = list(size = 8, color = "#8e44ad",
                             family = "IBM Plex Mono"),
        textposition  = "top center",
        hovertemplate = "<b>%{text}</b><extra></extra>",
        showlegend    = TRUE,
        inherit       = FALSE
      )
    }
    
    # ?????? Layout ????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
    p %>% layout(
      xaxis  = list(title = "NMDS1", zeroline = FALSE,
                    gridcolor = "#e9ecef", showgrid = TRUE),
      yaxis  = list(title = "NMDS2", zeroline = FALSE,
                    gridcolor = "#e9ecef", showgrid = TRUE,
                    scaleanchor = "x", scaleratio = 1),
      legend = list(orientation = "v", x = 1.01, y = 0.99,
                    font = list(size = 10)),
      margin = list(l = 10, r = 10, t = 10, b = 40),
      paper_bgcolor = "rgba(0,0,0,0)",
      plot_bgcolor  = "rgba(0,0,0,0)",
      font = list(family = "IBM Plex Sans")
    )
  })
}

# =============================================================================
# 4.  Run
# =============================================================================
shinyApp(ui = ui, server = server)