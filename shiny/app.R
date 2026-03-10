library(shiny)
library(ggplot2)
library(dplyr)
library(vegan)
library(ggrepel)

set.seed(123)

## Load prep objects from v_plants_prep.R --------------------------------------

prep_candidates <- c("v_plants_prep.R", "../v_plants_prep.R")
prep_path <- prep_candidates[file.exists(prep_candidates)][1]
if (is.na(prep_path)) {
  stop("Could not find v_plants_prep.R from working directory or its parent.")
}
# Source into the app's environment so base/utils functions are available
source(prep_path, local = TRUE)

## Candidate environmental variables for envfit (from prep)
vars_driver_strict_final <- vars_driver_strict_final
vars_driver_strict_final <- intersect(vars_driver_strict_final, names(meta_strict))

extra_numeric <- c("gap_years", "haeddypimetrar")
extra_numeric <- extra_numeric[extra_numeric %in% names(meta_strict)]

env_choices <- unique(c(vars_driver_strict_final, extra_numeric))

## Site scores aligned with meta_strict
site_scores <- as.data.frame(scores(nmds_strict, display = "sites")) %>%
  bind_cols(meta_strict %>% select(plot_number, sample))

## UI ---------------------------------------------------------------------------

ui <- fluidPage(
  titlePanel("v_plants NMDS explorer"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      helpText("Choose environmental variables to fit onto the NMDS ordination."),
      checkboxGroupInput(
        "env_vars",
        "Environmental variables",
        choices  = env_choices,
        selected = vars_driver_strict_final
      ),
      sliderInput(
        "p_thresh",
        "Show arrows with p ≤",
        min   = 0.001,
        max   = 0.20,
        value = 0.05,
        step  = 0.001
      ),
      checkboxInput(
        "show_factors",
        "Show factor centroids (categorical envfit results)",
        value = TRUE
      )
    ),
    mainPanel(
      width = 9,
      plotOutput(
        "nmds_plot",
        height = "600px",
        brush  = brushOpts(id = "nmds_brush")
      ),
      h4("Brushed points (first + last per plot)"),
      tableOutput("brush_table"),
      h4("Envfit summary for selected variables"),
      tableOutput("envfit_table")
    )
  )
)

## Server -----------------------------------------------------------------------

server <- function(input, output, session) {

  env_data <- reactive({
    vars <- input$env_vars
    if (is.null(vars) || length(vars) == 0) {
      return(NULL)
    }
    df <- meta_strict %>%
      select(all_of(vars)) %>%
      mutate(
        across(where(is.character), as.factor),
        across(where(is.logical), as.factor)
      )
    df
  })

  envfit_result <- reactive({
    df <- env_data()
    if (is.null(df)) {
      return(NULL)
    }
    vegan::envfit(nmds_strict, df, permutations = 999)
  })

  output$nmds_plot <- renderPlot({
    fit <- envfit_result()
    p_thresh <- input$p_thresh

    # Base points and connecting lines
    p <- ggplot(site_scores, aes(x = NMDS1, y = NMDS2)) +
      geom_line(
        aes(group = plot_number),
        colour = "grey70", linewidth = 0.3, alpha = 0.6
      ) +
      geom_point(
        aes(colour = sample, shape = sample),
        size = 2.2, alpha = 0.85
      ) +
      scale_colour_manual(
        values = c(first = "#2980b9", last = "#e67e22"),
        labels = c(first = "First survey", last = "Last survey")
      ) +
      scale_shape_manual(
        values = c(first = 16, last = 17),
        labels = c(first = "First survey", last = "Last survey")
      ) +
      coord_equal() +
      labs(
        title    = "NMDS of plant community composition (Bray–Curtis)",
        subtitle = sprintf(
          "%d plots · stress = %.3f",
          length(unique(site_scores$plot_number)),
          nmds_strict$stress
        ),
        x = "NMDS1", y = "NMDS2",
        colour = NULL, shape = NULL
      ) +
      theme_classic(base_size = 11) +
      theme(legend.position = "bottom")

    if (!is.null(fit)) {
      # Numeric vectors
      vec_df <- tibble(variable = character(), NMDS1 = numeric(), NMDS2 = numeric(),
                       r2 = numeric(), p_value = numeric(), kind = character())
      if (!is.null(fit$vectors) && !is.null(fit$vectors$arrows)) {
        vec_raw <- as.data.frame(fit$vectors$arrows) %>%
          rownames_to_column("variable")
        vec_axes <- setdiff(names(vec_raw), "variable")
        if (length(vec_axes) >= 2) {
          vec_df <- vec_raw %>%
            rename(NMDS1 = all_of(vec_axes[1]), NMDS2 = all_of(vec_axes[2])) %>%
            select(variable, NMDS1, NMDS2) %>%
            mutate(
              r2      = as.numeric(fit$vectors$r[variable]),
              p_value = as.numeric(fit$vectors$pvals[variable]),
              kind    = "numeric"
            )
        }
      }

      # Factor centroids (optional)
      fac_df <- tibble(variable = character(), NMDS1 = numeric(), NMDS2 = numeric(),
                       r2 = numeric(), p_value = numeric(), kind = character())
      if (input$show_factors &&
          !is.null(fit$factors) && !is.null(fit$factors$centroids)) {
        fac_raw <- as.data.frame(fit$factors$centroids) %>%
          rownames_to_column("variable")
        fac_axes <- setdiff(names(fac_raw), "variable")
        if (length(fac_axes) >= 2) {
          fac_df <- fac_raw %>%
            rename(NMDS1 = all_of(fac_axes[1]), NMDS2 = all_of(fac_axes[2])) %>%
            select(variable, NMDS1, NMDS2) %>%
            mutate(
              r2      = NA_real_,
              p_value = as.numeric(fit$factors$pvals[variable]),
              kind    = "factor"
            )
        }
      }

      arrow_df <- bind_rows(vec_df, fac_df) %>%
        filter(is.na(p_value) | p_value <= p_thresh)

      if (nrow(arrow_df) > 0) {
        ordi_radius <- max(abs(c(site_scores$NMDS1, site_scores$NMDS2))) * 0.7
        arrow_df <- arrow_df %>%
          mutate(
            xend = NMDS1 * ordi_radius,
            yend = NMDS2 * ordi_radius
          )

        p <- p +
          geom_segment(
            data = arrow_df %>% filter(kind == "numeric"),
            aes(x = 0, y = 0, xend = xend, yend = yend),
            arrow = arrow(length = unit(0.2, "cm")),
            colour = "#c0392b", linewidth = 0.7,
            inherit.aes = FALSE
          ) +
          geom_label_repel(
            data = arrow_df %>% filter(kind == "numeric"),
            aes(x = xend * 1.1, y = yend * 1.1, label = variable),
            size = 2.8, colour = "#c0392b",
            fill = alpha("white", 0.7), label.size = 0,
            max.overlaps = 30, seed = 42,
            inherit.aes = FALSE
          )
      }
    }

    p
  })

  output$brush_table <- renderTable({
    brushed <- brushedPoints(site_scores, input$nmds_brush)
    if (nrow(brushed) == 0) {
      return(NULL)
    }
    brushed %>%
      arrange(plot_number, sample) %>%
      select(plot_number, sample, everything())
  })

  output$envfit_table <- renderTable({
    fit <- envfit_result()
    if (is.null(fit)) {
      return(NULL)
    }

    vec_tbl <- tibble(
      variable = names(fit$vectors$r),
      r2       = as.numeric(fit$vectors$r),
      p_value  = as.numeric(fit$vectors$pvals),
      kind     = "numeric"
    )

    if (!is.null(fit$factors) && length(fit$factors$pvals) > 0) {
      fac_tbl <- tibble(
        variable = names(fit$factors$r),
        r2       = as.numeric(fit$factors$r),
        p_value  = as.numeric(fit$factors$pvals),
        kind     = "factor"
      )
    } else {
      fac_tbl <- tibble(variable = character(), r2 = numeric(),
                        p_value = numeric(), kind = character())
    }

    bind_rows(vec_tbl, fac_tbl) %>%
      arrange(p_value)
  })
}

shinyApp(ui, server)

