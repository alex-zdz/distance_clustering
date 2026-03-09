# =============================================================================
# Interactive Shiny Dashboard for 1D Clustering Simulation Study – robust version
# =============================================================================

library(shiny)
library(shinyWidgets)
library(DT)
library(ggplot2)
library(viridis)
library(tidyverse)

# Source functions (adjust paths if necessary)
print(getwd())
source("src/search_algorithm.R")
source("src/mixture_utils.R")
source("src/distances.R")

# Data generation functions (unchanged)
generate_scenario1 <- function() {
  set.seed(1)
  N <- 100
  K <- 3
  mu_true    <- c(-2, 0, 1)
  sd_true    <- c(1, 2, 1) * 0.2
  weights_true <- c(0.3, 0.5, 0.2)
  cluster_true <- rep(1:K, round(weights_true * N))
  y <- numeric(N)
  for (k in 1:K) {
    idx <- cluster_true == k
    ni  <- sum(idx)
    y[idx] <- rnorm(ni, mean = mu_true[k], sd = sd_true[k])
  }
  list(y = y, cluster_true = cluster_true, K = K)
}

generate_scenario2 <- function() {
  set.seed(2)
  N <- 150
  K <- 4
  mu_true    <- c(-4, -1, 1, 4)
  sd_true    <- rep(0.5, 4)
  weights_true <- rep(0.25, 4)
  cluster_true <- rep(1:K, round(weights_true * N))
  y <- numeric(N)
  for (k in 1:K) {
    idx <- cluster_true == k
    ni  <- sum(idx)
    y[idx] <- rnorm(ni, mean = mu_true[k], sd = sd_true[k])
  }
  list(y = y, cluster_true = cluster_true, K = K)
}

scen1_data <- generate_scenario1()
scen2_data <- generate_scenario2()

# ── Load results with safety ─────────────────────────────────────────────────
out_dir <- "results/simstudy/1D"
files <- list.files(out_dir, pattern = "\\.rds$", full.names = TRUE)

if (length(files) == 0) {
  message("No .rds files found in ", out_dir, " — dashboard will show placeholder")
  all_results <- tibble()  # empty tibble
} else {
  all_results <- map_dfr(files, readRDS, .id = "file")
  message("Loaded ", nrow(all_results), " result rows from ", length(files), " files")
}

# Add labels
all_results <- all_results %>%
  mutate(
    scenario_label = case_when(
      scenario == "scen1" ~ "Scenario 1: Moderately overlapping (K=3, unequal σ)",
      scenario == "scen2" ~ "Scenario 2: Well-separated (K=4, equal σ)",
      TRUE ~ as.character(scenario)
    ),
    clusterwise_label = ifelse(clusterwise, "Yes", "No")
  )

# ── Scenario plot function (unchanged) ───────────────────────────────────────
plot_scenario_data <- function(dat, scen_label) {
  df <- tibble(y = dat$y, cluster = as.factor(dat$cluster_true))
  ggplot(df, aes(x = y)) +
    geom_density(fill = "grey85", colour = "grey60", alpha = 0.5) +
    geom_rug(aes(colour = cluster), sides = "b", length = unit(0.05, "npc")) +
    geom_point(aes(y = 0, colour = cluster), size = 1.8, alpha = 0.85, 
               position = position_jitter(height = 0.01)) +
    #scale_colour_viridis_d(option = "plasma", name = "True cluster") +
    labs(title = scen_label,
         subtitle = sprintf("N = %d | K = %d components", length(dat$y), dat$K),
         x = "Data values (y)", y = "Density / Rug plot") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
}

# ── UI ───────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  titlePanel("Clustering Simulation Dashboard – 1D Gaussian Mixtures"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      pickerInput("scenario_select", "Scenario(s):",
                  choices = unique(all_results$scenario_label),
                  selected = unique(all_results$scenario_label),
                  multiple = TRUE, options = list(`actions-box` = TRUE)),
      
      pickerInput("method_select", "Method(s):",
                  choices = unique(all_results$method),
                  selected = unique(all_results$method),
                  multiple = TRUE, options = list(`actions-box` = TRUE)),
      
      checkboxGroupInput("clusterwise_select", "Clusterwise:",
                         choices = c("Yes", "No"), selected = c("Yes", "No")),
      
      pickerInput("version_select", "Version:",
                  choices = unique(all_results$version),
                  selected = unique(all_results$version),
                  multiple = TRUE),
      
      uiOutput("nsweet_slider"),
      uiOutput("nms_slider"),
      uiOutput("a0_slider"),
      
      downloadButton("download_summary", "Download filtered summary CSV")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Scenario Data",
                 fluidRow(
                   column(6, plotOutput("scen1_plot")),
                   column(6, plotOutput("scen2_plot"))
                 )),
        tabPanel("Performance Summary Table", DTOutput("summary_table")),
        tabPanel("ARI by Method & Clusterwise", plotOutput("plot_method_clusterwise", height = "650px")),
        tabPanel("ARI by n_sweet & n_ms", plotOutput("plot_nsweet_nms", height = "650px")),
        tabPanel("ARI by alpha_0", plotOutput("plot_alpha0", height = "650px"))
      )
    )
  )
)

# ── Server ───────────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  # Reactive filtered data
  filtered_results <- reactive({
    req(nrow(all_results) > 0)  # prevent execution if no data
    
    all_results %>%
      filter(
        scenario_label %in% input$scenario_select,
        method %in% input$method_select,
        clusterwise_label %in% input$clusterwise_select,
        version %in% input$version_select,
        n_sweet >= input$nsweet_range[1],
        n_sweet <= input$nsweet_range[2],
        n_ms >= input$nms_range[1],
        n_ms <= input$nms_range[2],
        alpha_0 >= input$a0_range[1],
        alpha_0 <= input$a0_range[2]
      )
  })
  
  # Dynamic sliders – defensive version
  output$nsweet_slider <- renderUI({
    req(nrow(all_results) > 0)
    ns_min <- min(all_results$n_sweet, na.rm = TRUE)
    ns_max <- max(all_results$n_sweet, na.rm = TRUE)
    sliderInput("nsweet_range", "n_sweet:",
                min = ns_min, max = ns_max,
                value = c(ns_min, ns_max),
                step = 1)
  })
  
  output$nms_slider <- renderUI({
    req(nrow(all_results) > 0)
    nms_min <- min(all_results$n_ms, na.rm = TRUE)
    nms_max <- max(all_results$n_ms, na.rm = TRUE)
    sliderInput("nms_range", "n_ms:",
                min = nms_min, max = nms_max,
                value = c(nms_min, nms_max),
                step = 1)
  })
  
  output$a0_slider <- renderUI({
    req(nrow(all_results) > 0)
    a0_min <- min(all_results$alpha_0, na.rm = TRUE)
    a0_max <- max(all_results$alpha_0, na.rm = TRUE)
    sliderInput("a0_range", HTML("α₀:"), 
                min = a0_min, max = a0_max,
                value = c(a0_min, a0_max),
                step = 1)
  })
  # Scenario plots (always shown)
  output$scen1_plot <- renderPlot({
    plot_scenario_data(scen1_data, "Scenario 1: Moderately overlapping (K=3, unequal σ)")
  })
  
  output$scen2_plot <- renderPlot({
    plot_scenario_data(scen2_data, "Scenario 2: Well-separated (K=4, equal σ)")
  })
  
  # Summary table
  output$summary_table <- renderDT({
    req(filtered_results())
    filtered_results() %>%
      group_by(scenario_label, method, clusterwise_label, version, n_sweet, n_ms, alpha_0) %>%
      summarise(
        `Mean ARI`     = mean(ARI, na.rm = TRUE) %>% round(3),
        `SD ARI`       = sd(ARI, na.rm = TRUE) %>% round(3),
        `Mean D`       = mean(final_D, na.rm = TRUE) %>% round(4),
        `Mean sweets`  = mean(n_accepted_sweets, na.rm = TRUE) %>% round(1),
        `Mean merges`  = mean(n_accepted_merges, na.rm = TRUE) %>% round(1),
        `Mean splits`  = mean(n_accepted_splits, na.rm = TRUE) %>% round(1),
        .groups = "drop"
      ) %>%
      datatable(options = list(pageLength = 15, scrollX = TRUE),
                extensions = "Buttons",
                filter = "top")
  })
  
  # Plots (only render if data exists)
  output$plot_method_clusterwise <- renderPlot({
    req(filtered_results(), nrow(filtered_results()) > 0)
    filtered_results() %>%
      ggplot(aes(x = method, y = ARI, colour = as.factor(n_sweet))) +
      geom_boxplot(alpha = 0.3, outlier.shape = NA) +
      geom_point(position = position_jitterdodge(), size = 2.8, alpha = 0.9) +
      facet_grid(scenario_label ~ clusterwise_label) +
      #scale_colour_viridis_d(option = "C") +
      labs(title = "ARI: method × clusterwise × n_sweet") +
      theme_minimal() + theme(legend.position = "bottom")
  })
  
  output$plot_nsweet_nms <- renderPlot({
    req(filtered_results(), nrow(filtered_results()) > 0)
    filtered_results() %>%
      ggplot(aes(x = as.factor(n_sweet), y = ARI, colour = as.factor(n_ms))) +
      geom_boxplot(alpha = 0.3, outlier.shape = NA) +
      geom_point(position = position_jitterdodge(), size = 2.8, alpha = 0.9) +
      facet_grid(scenario_label ~ method) +
      #scale_colour_viridis_d(option = "plasma") +
      labs(title = "ARI: n_sweet × n_ms × method") +
      theme_minimal() + theme(legend.position = "bottom")
  })
  
  output$plot_alpha0 <- renderPlot({
    req(filtered_results(), nrow(filtered_results()) > 0)
    filtered_results() %>%
      ggplot(aes(x = as.factor(alpha_0), y = ARI, colour = method)) +
      geom_boxplot(alpha = 0.3, outlier.shape = NA) +
      geom_point(position = position_jitterdodge(), size = 2.8, alpha = 0.9) +
      facet_grid(scenario_label ~ clusterwise_label) +
      #scale_colour_viridis_d(option = "viridis") +
      labs(title = "ARI: alpha_0 × method × clusterwise") +
      theme_minimal() + theme(legend.position = "bottom")
  })
  
  # Download
  output$download_summary <- downloadHandler(
    filename = "simulation_summary_filtered.csv",
    content = function(file) write_csv(filtered_results(), file)
  )
}

# Launch
shinyApp(ui, server)