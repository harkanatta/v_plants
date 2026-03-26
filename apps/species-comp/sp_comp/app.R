# =============================================================================
# v_plants ??? Species Composition Explorer  (v7)
# Single-file Shiny app
#
# Dependencies: shiny, tidyverse, vegan, bslib, DT
# Data:         sources v_plants_prep.R from the working directory
#
# Run with:     shiny::runApp("app_species_explorer.R")
# =============================================================================

library(shiny)
library(tidyverse)
library(vegan)
library(bslib)
library(DT)

# -----------------------------------------------------------------------------
# 0.  Load data
# -----------------------------------------------------------------------------
if (!exists("community_long_paired")) {
  dat <- readRDS("v_plants_app_data.rds")
  list2env(dat, envir = .GlobalEnv)
}

# -----------------------------------------------------------------------------
# 1.  Derived objects
# -----------------------------------------------------------------------------

if (!exists("outlier_plots")) {
  nmds_scores_tmp <- as.data.frame(vegan::scores(nmds_strict, display = "sites")) %>%
    bind_cols(meta_strict %>% select(plot_number, sample))
  centroid_tmp <- colMeans(nmds_scores_tmp[, c("NMDS1", "NMDS2")])
  outlier_plots <- nmds_scores_tmp %>%
    mutate(dist_c = sqrt((NMDS1 - centroid_tmp[1])^2 + (NMDS2 - centroid_tmp[2])^2)) %>%
    arrange(desc(dist_c)) %>%
    distinct(plot_number, .keep_all = TRUE) %>%
    slice_head(n = 3) %>%
    pull(plot_number)
}

meta_lookup <- meta_strict %>%
  filter(sample == "first") %>%
  select(plot_number,
         any_of(c("soil_type_name", "moisture_name", "topography_name",
                  "habitat_type_name", "haeddypimetrar", "gap_years"))) %>%
  distinct(plot_number, .keep_all = TRUE)

all_plots         <- sort(unique(meta_strict$plot_number))
non_outlier_plots <- setdiff(all_plots, outlier_plots)

mk_named <- function(x) c("All" = "All", setNames(x, x))
soil_types     <- mk_named(sort(unique(na.omit(meta_lookup$soil_type_name))))
moisture_types <- mk_named(sort(unique(na.omit(meta_lookup$moisture_name))))
habitat_types  <- mk_named(sort(unique(na.omit(meta_lookup$habitat_type_name))))

bc_per_plot <- community_long_paired %>%
  select(plot_number, sample, taxon_name, cover) %>%
  pivot_wider(names_from = taxon_name, values_from = cover, values_fill = 0) %>%
  group_by(plot_number) %>%
  group_modify(function(df, key) {
    mat <- df %>% select(-sample) %>% as.matrix()
    if (nrow(mat) < 2) return(tibble(bc = NA_real_))
    tibble(bc = as.numeric(vegdist(mat, method = "bray"))[1])
  }) %>%
  ungroup()

bc_median <- median(bc_per_plot$bc, na.rm = TRUE)
bc_q75    <- quantile(bc_per_plot$bc, 0.75, na.rm = TRUE)

bc_label <- function(bc) {
  case_when(is.na(bc)      ~ "\u2014",
            bc < bc_median ~ "low turnover",
            bc < bc_q75    ~ "moderate turnover",
            TRUE           ~ "high turnover")
}

bc_col <- function(bc) {
  case_when(is.na(bc)      ~ "#6c757d",
            bc < bc_median ~ "#2c6e49",
            bc < bc_q75    ~ "#e67e22",
            TRUE           ~ "#c0392b")
}

diversity_stats <- function(cover) {
  x <- cover[cover > 0]
  if (length(x) < 2) return(list(H = NA_real_, J = NA_real_))
  p <- x / sum(x)
  H <- -sum(p * log(p))
  J <- H / log(length(x))
  list(H = H, J = J)
}

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
    .nav-row { display:flex; gap:0.4rem; margin-top:0.5rem; }
    .nav-row .btn { flex:1; font-size:0.8rem; padding:0.25rem 0; }

    .sec-label {
      font-size:0.7rem; font-weight:600; color:#6c757d;
      text-transform:uppercase; letter-spacing:0.05em;
      margin:0.9rem 0 0.25rem;
    }

    .plot-info-bar {
      background:#f4f6f4;
      border:1px solid #dee2e6;
      border-radius:6px;
      padding:0.5rem 0.9rem;
      margin-bottom:0.75rem;
      line-height:1.55;
    }
    .plot-info-bar .line1 {
      display:flex; flex-wrap:wrap; align-items:baseline; gap:0 0.45rem;
      font-size:0.88rem;
    }
    .plot-info-bar .line2 {
      display:flex; flex-wrap:wrap; align-items:baseline; gap:0 0.45rem;
      font-size:0.78rem; color:#495057; margin-top:0.1rem;
    }
    .plot-id {
      font-family:'IBM Plex Mono',monospace;
      font-weight:700; font-size:1rem; color:#2c6e49;
    }
    .sep  { color:#adb5bd; }
    .mono { font-family:'IBM Plex Mono',monospace; font-weight:600; }

    .dataTables_wrapper { font-size:0.82rem; }
    table.dataTable thead th {
      background:#f4f6f4; color:#2c6e49;
      font-family:'IBM Plex Mono',monospace;
    }
  "))),
  
  # ?????? Sidebar ????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
  sidebar = sidebar(
    width = 255,
    
    tags$hr(style = "margin:0.3rem 0 0.5rem;"),
    
    div(class = "sec-label", "Filter"),
    selectInput("filter_habitat",  NULL, choices = habitat_types,  selected = "All"),
    selectInput("filter_soil",     NULL, choices = soil_types,     selected = "All"),
    selectInput("filter_moisture", NULL, choices = moisture_types, selected = "All"),
    checkboxInput("exclude_outliers", "Exclude coastal outliers", value = TRUE),
    
    tags$hr(style = "margin:0.6rem 0 0.4rem;"),
    
    div(class = "sec-label", "Plot"),
    selectInput("plot_id", NULL,
                choices  = non_outlier_plots,
                selected = non_outlier_plots[1]),
    div(class = "nav-row",
        actionButton("btn_prev", "\u2190 prev",
                     class = "btn btn-outline-secondary btn-sm"),
        actionButton("btn_next", "next \u2192",
                     class = "btn btn-outline-secondary btn-sm")
    )
  ),
  
  # ?????? Main panel ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
  uiOutput("plot_info_bar"),
  
  card(
    card_body(
      DT::dataTableOutput("species_dt")
    ),
    full_screen = TRUE,
    fill = FALSE
  )
  
) # end page_sidebar

# =============================================================================
# 3.  Server
# =============================================================================

server <- function(input, output, session) {
  
  # ?????? Filtered plot list ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
  filtered_plots <- reactive({
    plots <- if (input$exclude_outliers) non_outlier_plots else all_plots
    
    if (input$filter_habitat != "All" && "habitat_type_name" %in% names(meta_lookup))
      plots <- intersect(plots,
                         meta_lookup %>% filter(habitat_type_name == input$filter_habitat) %>%
                           pull(plot_number))
    if (input$filter_soil != "All" && "soil_type_name" %in% names(meta_lookup))
      plots <- intersect(plots,
                         meta_lookup %>% filter(soil_type_name == input$filter_soil) %>%
                           pull(plot_number))
    if (input$filter_moisture != "All" && "moisture_name" %in% names(meta_lookup))
      plots <- intersect(plots,
                         meta_lookup %>% filter(moisture_name == input$filter_moisture) %>%
                           pull(plot_number))
    sort(plots)
  })
  
  observeEvent(filtered_plots(), {
    fp  <- filtered_plots()
    cur <- isolate(input$plot_id)
    updateSelectInput(session, "plot_id",
                      choices  = fp,
                      selected = if (cur %in% fp) cur else fp[1])
  })
  
  # ?????? Prev / Next ????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
  observeEvent(input$btn_prev, {
    fp  <- filtered_plots()
    idx <- match(input$plot_id, fp)
    if (!is.na(idx) && idx > 1)
      updateSelectInput(session, "plot_id", choices = fp, selected = fp[idx - 1])
  })
  
  observeEvent(input$btn_next, {
    fp  <- filtered_plots()
    idx <- match(input$plot_id, fp)
    if (!is.na(idx) && idx < length(fp))
      updateSelectInput(session, "plot_id", choices = fp, selected = fp[idx + 1])
  })
  
  # ?????? Plot data ??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
  plot_data <- reactive({
    req(input$plot_id)
    
    community_long_paired %>%
      filter(plot_number == input$plot_id) %>%
      select(sample, taxon_name, cover) %>%
      pivot_wider(names_from = sample, values_from = cover, values_fill = 0) %>%
      mutate(
        first  = coalesce(first, 0),
        last   = coalesce(last,  0),
        change = last - first
      ) %>%
      filter(first > 0 | last > 0) %>%
      arrange(desc(pmax(first, last)))
  })
  
  # ?????? Plot info bar ??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
  output$plot_info_bar <- renderUI({
    req(input$plot_id)
    
    row       <- meta_lookup %>% filter(plot_number == input$plot_id)
    bc        <- bc_per_plot %>% filter(plot_number == input$plot_id) %>% pull(bc)
    bc        <- if (length(bc) == 0) NA_real_ else bc[1]
    pd        <- plot_data()
    rich_first <- pd %>% filter(first > 0) %>% nrow()
    rich_last  <- pd %>% filter(last  > 0) %>% nrow()
    div_first  <- diversity_stats(pd$first)
    div_last   <- diversity_stats(pd$last)
    
    fmt <- function(x) if (is.na(x)) "?" else sprintf("%.2f", x)
    sep <- tags$span(class = "sep", "\u00b7")
    arr <- "\u2192"
    
    # safe field extractor for character cols
    field <- function(col) {
      if (col %in% names(row) && !is.na(row[[col]])) as.character(row[[col]]) else NULL
    }
    
    # formatted numeric fields
    elev_str <- if ("haeddypimetrar" %in% names(row) && !is.na(row$haeddypimetrar))
      paste0(round(row$haeddypimetrar), " m.y.s") else NULL
    gap_str  <- if ("gap_years" %in% names(row) && !is.na(row$gap_years))
      paste0(row$gap_years, " yr gap") else NULL
    
    # line 1: plot ID ?? habitat ?? soil ?? moisture ?? topography ?? elevation ?? gap
    line1_items <- list(tags$span(class = "plot-id", input$plot_id))
    for (val in Filter(Negate(is.null), list(
      field("habitat_type_name"),
      field("soil_type_name"),
      field("moisture_name"),
      field("topography_name"),
      elev_str,
      gap_str
    ))) {
      line1_items <- c(line1_items, list(sep, tags$span(val)))
    }
    
    # line 2: BC ?? turnover ?? Richness ?? Shannon ?? Evenness
    bc_str    <- if (is.na(bc)) "\u2014" else sprintf("%.3f", bc)
    bc_colour <- bc_col(bc)
    
    line2 <- div(class = "line2",
                 tags$span("BC:"),
                 tags$span(class = "mono",
                           style = paste0("color:", bc_colour, ";"),
                           bc_str),
                 sep,
                 tags$span(style = paste0("color:", bc_colour, ";"),
                           bc_label(bc)),
                 sep,
                 tags$span(sprintf("Richness: %d %s %d", rich_first, arr, rich_last)),
                 sep,
                 tags$span(sprintf("Shannon: %s %s %s",
                                   fmt(div_first$H), arr, fmt(div_last$H))),
                 sep,
                 tags$span(sprintf("Evenness: %s %s %s",
                                   fmt(div_first$J), arr, fmt(div_last$J)))
    )
    
    div(class = "plot-info-bar",
        div(class = "line1", line1_items),
        line2
    )
  })
  
  # ?????? DT species table ?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
  output$species_dt <- DT::renderDataTable({
    
    df <- plot_data() %>%
      transmute(
        Taxon       = taxon_name,
        `First (%)` = round(first,  1),
        `Last (%)`  = round(last,   1),
        Change      = round(change, 1)
      )
    
    DT::datatable(
      df,
      rownames  = FALSE,
      selection = "none",
      class     = "stripe hover compact",
      options   = list(
        pageLength     = -1,
        dom            = "ft",
        order          = list(list(2L, "desc")),
        scrollY        = "65vh",
        scrollCollapse = TRUE,
        columnDefs     = list(
          list(className = "dt-right", targets = 1:3),
          list(width = "52%", targets = 0)
        )
      )
    ) %>%
      DT::formatStyle(
        "Change",
        color = DT::styleInterval(
          c(-0.001, 0.001),
          c("#c0392b", "#6c757d", "#2c6e49")
        ),
        fontWeight = DT::styleInterval(
          c(-0.001, 0.001),
          c("600", "400", "600")
        )
      ) %>%
      DT::formatStyle(
        c("First (%)", "Last (%)"),
        background         = DT::styleColorBar(c(0, 100), "#e9f2ec"),
        backgroundSize     = "98% 55%",
        backgroundRepeat   = "no-repeat",
        backgroundPosition = "center"
      )
    
  }, server = FALSE)
}

# =============================================================================
# 4.  Run
# =============================================================================
shinyApp(ui = ui, server = server)