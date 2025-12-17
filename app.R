# MIT License
# 
# Copyright (c) 2025 Frank Simmen
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#   
#   The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

data_url <- "https://github.com/fsimghub/populism-shiny/releases/download/v1.0-data/ShinyPopDat.qs"

# Download if not already present (runs only once per deployment/session as needed)
if (!file.exists("data/ShinyPopDat.qs")) {
  dir.create("data", showWarnings = FALSE)
  download.file(data_url, "data/ShinyPopDat.qs", mode = "wb")
}

rm(list=ls())
gc()

# Check if the 'pacman' package is installed; if not, install it from CRAN
if (!require("pacman", character.only = TRUE)) {
  install.packages("pacman")
}

pacman::p_load(
  data.table,
  dplyr,
  ggnewscale,
  ggplot2,
  grid,
  cowplot,
  viridis,
  patchwork,
  parallel,
  qs,
  shiny,
  rlang
)


###############################
### FUNCTIONS
###############################

plot_effects <- function(dt,
                         in_casename,
                         in_covariates,
                         in_method,
                         in_preperiod,
                         in_variable,
                         in_database,
                         baseval = NULL,
                         set_backtracing_3yr = FALSE,
                         set_backtracing_5yr = FALSE,
                         set_case_permutation = FALSE,
                         set_loo_permutation = FALSE,
                         set_ci = FALSE,
                         series = "att") {
  
  series <- match.arg(series, c("att", "level", "origlvl"))
  
  # Force case permutation off if not in "att" mode
  if (series != "att" && set_case_permutation) {
    message("set_case_permutation is only available in series = 'att'. Forcing it to FALSE.")
    set_case_permutation <- FALSE
  }
  
  # Subset the data ---------------------------------------------------------
  selected <- dt[casename == in_casename &
                   covariates == in_covariates &
                   method == in_method &
                   preperiod == in_preperiod &
                   variable == in_variable &
                   database == in_database]
  
  if (nrow(selected) == 0) stop("No data found for the selected combination.")
  
  # Variable labels ---------------------------------------------------------
  var_labels <- c(
    "Pre-Tax Income Share: Bottom 50%",
    "Pre-Tax Income Share: Middle 40%",
    "Pre-Tax Income Share: Top 10%",
    "Pre-Tax Income Share: Top 1%",
    "Disposable Income Share: Bottom 50%",
    "Disposable Income Share: Middle 40%",
    "Disposable Income Share: Top 10%",
    "Disposable Income Share: Top 1%",
    "Pre-Tax National Income Share (WID)",
    "Disposable National Income Share (WID)",
    "Disposable Gini (SWIID)",
    "Market Gini (SWIID)",
    "Real GDP per Capita",
    "Inflation Rate",
    "Government Debt / GDP",
    "Institutions",
    "Civil Society"
  )
  names(var_labels) <- c(
    "sptinc992j_p0p50", "sptinc992j_p50p90", "sptinc992j_p90p100", "sptinc992j_p99p100",
    "sdiinc992j_p0p50", "sdiinc992j_p50p90", "sdiinc992j_p90p100", "sdiinc992j_p99p100",
    "gptinc992j", "gdiinc992j", "swiid_gini_disp", "swiid_gini_mkt",
    "mad_gdppc", "inflrate", "gmd_govdebt_GDP", "institution_pc1", "society_pc1"
  )
  plot_subtitle <- var_labels[in_variable]
  if (is.na(plot_subtitle)) plot_subtitle <- in_variable
  
  # Title -------------------------------------------------------------------
  treat_year <- as.integer(substr(in_casename, nchar(in_casename) - 3, nchar(in_casename)))
  plot_title <- in_casename
  if (!is.null(baseval) && is.data.table(baseval)) {
    base_row <- baseval[casename == in_casename]
    if (nrow(base_row) > 0) {
      ctry <- base_row$ctryname[1]
      leader <- base_row$leadershort[1]
      if (!is.na(ctry) && !is.na(leader) && !is.na(treat_year)) {
        plot_title <- paste0(ctry, " (", leader, ", ", treat_year, ")")
      }
    }
  }
  
  # Base value for origlvl --------------------------------------------------
  base_value <- 0
  if (series == "origlvl") {
    if (is.null(baseval) || !is.data.table(baseval)) {
      stop("For series = 'origlvl', 'baseval' must be a data.table with base values.")
    }
    base_row <- baseval[casename == in_casename]
    if (nrow(base_row) > 0 && in_variable %in% names(base_row)) {
      base_value <- as.numeric(base_row[[in_variable]])
    }
  }
  
  # Level transformations ---------------------------------------------------
  if (series %in% c("level", "origlvl")) {
    selected[, `:=`(
      treated_level       = observed - ATT,
      ci_low_level        = observed - ATT_high,
      ci_high_level       = observed - ATT_low,
      backtrace_3yr_level = ifelse(!is.na(backtrace_3yr_effect), observed - backtrace_3yr_effect, NA),
      backtrace_5yr_level = ifelse(!is.na(backtrace_5yr_effect), observed - backtrace_5yr_effect, NA)
    )]
    if (series == "origlvl") {
      selected[, `:=`(
        observed            = observed + base_value,
        treated_level       = treated_level + base_value,
        ci_low_level        = ci_low_level + base_value,
        ci_high_level       = ci_high_level + base_value,
        backtrace_3yr_level = backtrace_3yr_level + base_value,
        backtrace_5yr_level = backtrace_5yr_level + base_value
      )]
    }
  }
  
  # Grid --------------------------------------------------------------------
  time_range <- range(selected$time)
  grid_5 <- seq(5 * floor(time_range[1]/5), 5 * ceiling(time_range[2]/5), by = 5)
  
  # Start plot --------------------------------------------------------------
  p <- ggplot(selected, aes(x = time))
  
  # Main lines --------------------------------------------------------------
  if (series == "att") {
    p <- p + geom_line(aes(y = ATT, linetype = "ATT"), 
                       color = "black", linewidth = 0.8)
  } else {
    p <- p +
      geom_line(aes(y = observed, linetype = "Observed"), 
                color = "black", linewidth = 0.8) +
      geom_line(aes(y = treated_level, linetype = "Counterfactual"), 
                color = "black", linewidth = 0.8)
  }
  
  # Confidence ribbon -------------------------------------------------------
  if (set_ci) {
    if (series == "att") {
      p <- p + geom_ribbon(aes(ymin = ATT_low, ymax = ATT_high, fill = "95% CI"),
                           alpha = 0.3)
    } else {
      p <- p + geom_ribbon(aes(ymin = ci_low_level, ymax = ci_high_level, fill = "95% CI"),
                           alpha = 0.3)
    }
  }
  
  # Zero line (only in att) — fixed to avoid warning ------------------------
  if (series == "att") {
    p <- p + geom_hline(yintercept = 0, linetype = "solid", color = "grey60", linewidth = 0.5,
                        aes(linetype = "Zero"))
  }
  
  # Event line --------------------------------------------------------------
  p <- p + geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.7)
  
  # Backtracing -------------------------------------------------------------
  if (set_backtracing_3yr) {
    y_var <- if (series == "att") "backtrace_3yr_effect" else "backtrace_3yr_level"
    p <- p + geom_line(aes(y = .data[[y_var]], linetype = "Backtrace (3yr)"),
                       color = "#1b9e77", alpha = 0.6, linewidth = 0.8)
    p <- p + geom_vline(xintercept = -3, linetype = "longdash", color = "#1b9e77", linewidth = 0.8)
  }
  if (set_backtracing_5yr) {
    y_var <- if (series == "att") "backtrace_5yr_effect" else "backtrace_5yr_level"
    p <- p + geom_line(aes(y = .data[[y_var]], linetype = "Backtrace (5yr)"),
                       color = "#1b9e77", alpha = 0.6, linewidth = 0.8)
    p <- p + geom_vline(xintercept = -5, linetype = "dotdash", color = "#1b9e77", linewidth = 0.8)
  }
  
  # Case permutations (only allowed in att) ---------------------------------
  has_permutations <- FALSE
  if (set_case_permutation && length(selected$placebo_effects[[1]]) > 0) {
    K <- length(selected$placebo_effects[[1]])
    placebo_dt <- data.table(
      time  = rep(selected$time, each = K),
      value = unlist(selected$placebo_effects),
      group = factor(rep(1:K, times = nrow(selected)))
    )
    has_permutations <- TRUE
    p <- p + new_scale_color() +
      geom_line(data = placebo_dt, 
                aes(x = time, y = value, color = "Placebo Cases", group = group),
                alpha = 0.2, linewidth = 0.3)
  }
  
  # LOO permutations --------------------------------------------------------
  if (set_loo_permutation && length(selected$loo_effects[[1]]) > 0) {
    M <- length(selected$loo_effects[[1]])
    vals <- if (series == "att") unlist(selected$loo_effects) else 
      rep(selected$observed, each = M) - unlist(selected$loo_effects)
    loo_dt <- data.table(
      time  = rep(selected$time, each = M),
      value = vals,
      group = factor(rep(1:M, times = nrow(selected)))
    )
    if (!has_permutations) p <- p + new_scale_color()
    p <- p +
      geom_line(data = loo_dt, 
                aes(x = time, y = value, color = "Leave-One-Out", group = group),
                alpha = 0.2, linewidth = 0.3)
  }
  
  # Y-label -----------------------------------------------------------------
  y_lab <- switch(series,
                  att     = "Treatment Effect",
                  level   = "Normalized Series",
                  origlvl = "Original Level")
  
  # Legend scales -----------------------------------------------------------
  p <- p +
    scale_linetype_manual(
      values = c(
        "ATT"             = "solid",
        "Observed"        = "solid",
        "Counterfactual"  = "dashed",
        "Zero"            = "solid",
        "Backtrace (3yr)" = "longdash",
        "Backtrace (5yr)" = "dotdash"
      ),
      breaks = {
        brks <- character()
        if (series == "att") brks <- c("ATT", "Zero")
        else brks <- c("Observed", "Counterfactual")
        if (set_backtracing_3yr) brks <- c(brks, "Backtrace (3yr)")
        if (set_backtracing_5yr) brks <- c(brks, "Backtrace (5yr)")
        brks
      }
    ) +
    scale_fill_manual(
      values = c("95% CI" = "grey40"),
      guide = guide_legend(override.aes = list(alpha = 0.3))
    ) +
    scale_color_manual(
      values = c("Placebo Cases" = "steelblue", "Leave-One-Out" = "darkred"),
      guide = guide_legend(override.aes = list(alpha = 0.5, linewidth = 0.8))
    )
  
  # Final plot --------------------------------------------------------------
  p <- p +
    geom_vline(xintercept = grid_5[grid_5 != 0], linetype = "dotted", color = "grey80", linewidth = 0.5) +
    labs(title = plot_title,
         subtitle = plot_subtitle,
         x = "Event time (t)",
         y = y_lab,
         linetype = "Series",
         fill     = "Confidence",
         color    = "Permutations") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = rel(1.1)),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = margin(t = 10),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  if (grepl("\\%", plot_subtitle)) {
    p <- p + labs(subtitle = parse(text = plot_subtitle))
  }
  
  return(p)
}


###############################
### DATA AND SHINY APP
###############################

ShinyPopDat <- qread(file = "data/ShinyPopDat.qs")

# Variable labels (human-readable, for display)
var_labels <- c(
  "Pre-Tax Income Share: Bottom 50%",
  "Pre-Tax Income Share: Middle 40%",
  "Pre-Tax Income Share: Top 10%",
  "Pre-Tax Income Share: Top 1%",
  "Disposable Income Share: Bottom 50%",
  "Disposable Income Share: Middle 40%",
  "Disposable Income Share: Top 10%",
  "Disposable Income Share: Top 1%",
  "Pre-Tax National Income Share (WID)",
  "Disposable National Income Share (WID)",
  "Disposable Gini (SWIID)",
  "Market Gini (SWIID)",
  "Real GDP per Capita",
  "Inflation Rate",
  "Government Debt / GDP",
  "Institutions",
  "Civil Society"
)

# Corresponding internal variable codes
var_names <- c(
  "sptinc992j_p0p50", "sptinc992j_p50p90", "sptinc992j_p90p100", "sptinc992j_p99p100",
  "sdiinc992j_p0p50", "sdiinc992j_p50p90", "sdiinc992j_p90p100", "sdiinc992j_p99p100",
  "gptinc992j", "gdiinc992j", "swiid_gini_disp", "swiid_gini_mkt",
  "mad_gdppc", "inflrate", "gmd_govdebt_GDP", "institution_pc1", "society_pc1"
)

# Mapping from displayed label to internal code
label_to_code <- setNames(var_names, var_labels)

ui <- fluidPage(
  
  fluidRow(
    column(12,
           uiOutput("plots_ui")
    )
  ),
  
  hr(),
  
  fluidRow(
    column(6,
           radioButtons("view_mode", "View Mode",
                        choices = c("Dual Plot" = "dual", "Single Plot" = "single"),
                        selected = "dual", inline = TRUE)
    ),
    column(6,
           conditionalPanel(
             condition = "input.view_mode == 'dual'",
             radioButtons("control_mode", "Control Mode (Dual View)",
                          choices = c("Shared", "Separate"),
                          selected = "Shared", inline = TRUE)
           )
    )
  ),
  
  # Single Plot mode
  conditionalPanel(
    condition = "input.view_mode == 'single'",
    fluidRow(
      column(6,
             h4("Parameters"),
             selectInput("casename_single", "Case Name",
                         choices = sort(unique(ShinyPopDat$caselevel_dt$casename)),
                         selected = "BRA2019"),
             selectInput("variable_single", "Outcome Variable",
                         choices = character(0)),
             selectInput("covariates_single", "Covariates",
                         choices = sort(unique(ShinyPopDat$caselevel_dt$covariates)),
                         selected = "none"),
             selectInput("method_single", "Method",
                         choices = character(0)),
             selectInput("preperiod_single", "Preperiod",
                         choices = character(0)),
             selectInput("database_single", "Database",
                         choices = character(0))
      ),
      column(6,
             h4("Display Options"),
             radioButtons("series_single", "Series Mode",
                          choices = c("Treatment Effect" = "att",
                                      "Counterfactual"  = "level"),
                          selected = "att", inline = TRUE),
             checkboxInput("set_ci_single", "Show 95% Confidence Interval", value = FALSE),
             checkboxInput("set_backtracing_3yr_single", "Backtrace (3 years)", value = FALSE),
             checkboxInput("set_backtracing_5yr_single", "Backtrace (5 years)", value = FALSE),
             checkboxInput("set_loo_permutation_single", "Leave-One-Donor-Out Permutations", value = FALSE),
             conditionalPanel(
               condition = "input.series_single == 'att'",
               checkboxInput("set_case_permutation_single", "Donor-Placebo Permutations", value = FALSE)
             )
      )
    )
  ),
  
  # Dual Shared mode
  conditionalPanel(
    condition = "input.view_mode == 'dual' && input.control_mode == 'Shared'",
    fluidRow(
      column(6,
             selectInput("casename_shared", "Case Name (Left)",
                         choices = sort(unique(ShinyPopDat$caselevel_dt$casename)),
                         selected = "BRA2019"),
             selectInput("variable_shared", "Outcome Variable (Left)",
                         choices = character(0)),
             selectInput("covariates_shared", "Covariates",
                         choices = sort(unique(ShinyPopDat$caselevel_dt$covariates)),
                         selected = "none"),
             selectInput("method_shared", "Method",
                         choices = character(0)),
             selectInput("preperiod_shared", "Preperiod",
                         choices = character(0)),
             selectInput("database_shared", "Database",
                         choices = character(0))
      ),
      column(6,
             selectInput("casename_right_shared", "Case Name (Right)",
                         choices = sort(unique(ShinyPopDat$caselevel_dt$casename)),
                         selected = "BGR2009"),
             selectInput("variable_right_shared", "Outcome Variable (Right)",
                         choices = character(0)),
             radioButtons("series_shared", "Series Mode",
                          choices = c("Treatment Effect" = "att",
                                      "Counterfactual"  = "level"),
                          selected = "att", inline = TRUE),
             checkboxInput("set_ci_shared", "Show 95% Confidence Interval", value = FALSE),
             checkboxInput("set_backtracing_3yr_shared", "Backtrace (3 years)", value = FALSE),
             checkboxInput("set_backtracing_5yr_shared", "Backtrace (5 years)", value = FALSE),
             checkboxInput("set_loo_permutation_shared", "Leave-One-Donor-Out Permutations", value = FALSE),
             conditionalPanel(
               condition = "input.series_shared == 'att'",
               checkboxInput("set_case_permutation_shared", "Donor-Placebo Permutations", value = FALSE)
             )
      )
    )
  ),
  
  # Dual Separate mode
  conditionalPanel(
    condition = "input.view_mode == 'dual' && input.control_mode == 'Separate'",
    fluidRow(
      column(6,
             selectInput("casename_left", "Case Name (Left)",
                         choices = sort(unique(ShinyPopDat$caselevel_dt$casename)),
                         selected = "BRA2019"),
             selectInput("variable_left", "Outcome Variable (Left)",
                         choices = character(0)),
             selectInput("covariates_left", "Covariates",
                         choices = sort(unique(ShinyPopDat$caselevel_dt$covariates)),
                         selected = "none"),
             selectInput("method_left", "Method",
                         choices = character(0)),
             selectInput("preperiod_left", "Preperiod",
                         choices = character(0)),
             selectInput("database_left", "Database",
                         choices = character(0)),
             radioButtons("series_left", "Series Mode",
                          choices = c("Treatment Effect" = "att",
                                      "Counterfactual"  = "level"),
                          selected = "att", inline = TRUE),
             checkboxInput("set_ci_left", "Show 95% Confidence Interval", value = FALSE),
             checkboxInput("set_backtracing_3yr_left", "Backtrace (3 years)", value = FALSE),
             checkboxInput("set_backtracing_5yr_left", "Backtrace (5 years)", value = FALSE),
             checkboxInput("set_loo_permutation_left", "Leave-One-Donor-Out Permutations", value = FALSE),
             conditionalPanel(
               condition = "input.series_left == 'att'",
               checkboxInput("set_case_permutation_left", "Donor-Placebo Permutations", value = FALSE)
             )
      ),
      column(6,
             selectInput("casename_right", "Case Name (Right)",
                         choices = sort(unique(ShinyPopDat$caselevel_dt$casename)),
                         selected = "BGR2009"),
             selectInput("variable_right", "Outcome Variable (Right)",
                         choices = character(0)),
             selectInput("covariates_right", "Covariates",
                         choices = sort(unique(ShinyPopDat$caselevel_dt$covariates)),
                         selected = "none"),
             selectInput("method_right", "Method",
                         choices = character(0)),
             selectInput("preperiod_right", "Preperiod",
                         choices = character(0)),
             selectInput("database_right", "Database",
                         choices = character(0)),
             radioButtons("series_right", "Series Mode",
                          choices = c("Treatment Effect" = "att",
                                      "Counterfactual"  = "level"),
                          selected = "att", inline = TRUE),
             checkboxInput("set_ci_right", "Show 95% Confidence Interval", value = FALSE),
             checkboxInput("set_backtracing_3yr_right", "Backtrace (3 years)", value = FALSE),
             checkboxInput("set_backtracing_5yr_right", "Backtrace (5 years)", value = FALSE),
             checkboxInput("set_loo_permutation_right", "Leave-One-Donor-Out Permutations", value = FALSE),
             conditionalPanel(
               condition = "input.series_right == 'att'",
               checkboxInput("set_case_permutation_right", "Donor-Placebo Permutations", value = FALSE)
             )
      )
    )
  ),
  
  hr(),
  
  fluidRow(
    column(12,
           h4("Notes"),
           HTML("<p>This Shiny App visualizes causal estimates of the effects of populist leaders on economic growth and distributional outcomes (e.g., inequality) using the traditional synthetic control method (SCM) and the augmented synthetic control method (ASCM), where the outcome model employs matrix completion.</p>
            <p>Users can select specific populist episodes, outcomes, and donor pools from three databases (PLE, GDP, or combined). For each case, the app displays either (i) observed versus counterfactual outcome trajectories or (ii) the estimated treatment effect series.</p>
            <p>Available options include:</p>
            <ul>
              <li>95% jackknife+ confidence bands (<a href='https://arxiv.org/abs/1905.02928' target='_blank'>Barber et al., 2021, Ann. Statist.</a>).</li>
              <li>Backtrace (3 years) and Backtrace (5 years): in-time placebo exercises that counterfactually advance the treatment date by 3 or 5 years and estimate the corresponding treatment effect and counterfactual.</li>
              <li>Leave-One-Donor-Out permutations: robustness check omitting one donor at a time.</li>
              <li>Donor-Placebo permutations: placebo treatments assigned counterfactually to donor units (not available for counterfactual series).</li>
            </ul>
            <p>The app supports single-panel view or Dual Plot View Mode (side-by-side panels). In Dual Plot View, users choose Shared Mode (identical options across panels) or Separate Mode (independent options per panel) for direct method or specification comparisons.</p>
            <p>Cases and options are restricted to those meeting pre-estimation requirements; availability varies by database and data coverage.</p>
            <p>Case names follow the scheme ISO3 country code and treatment year, indicating when the populist comes to power. Detailed information on the individual populist cases is in the Appendix of the accompanying paper.</p>")
    )
  ),
 
  tags$footer(
    style = "
    text-align: center;
    font-size: 12px;
    color: #666;
    padding: 20px 10px;
    margin-top: 40px;
    border-top: 1px solid #eee;
    background-color: #f9f9f9;
  ",
    HTML("© 2025 Frank Simmen. Released under the MIT License. Free to use, modify, and distribute with attribution. No warranty provided."),
    tags$br(),
    tags$p(
      style = "margin-top: 8px;",
      "If you use this app or its outputs in academic work, please cite: ",
      tags$em("Frank Simmen (2026). Does Populism Cause Economic Outcomes?, Working Paper.")
    )
  )
)


server <- function(input, output, session) {
  
  # Helper functions
  label_method <- function(meth) {
    ifelse(meth == "asyn", "ASCM MCPanel",
           ifelse(meth == "syn", "SCM", meth))
  }
  
  label_preperiod <- function(pre) {
    ifelse(pre == "class1", "15 years",
           ifelse(pre == "class2", "10 years", pre))
  }
  
  label_database <- function(db) {
    ifelse(db == "both", "Combined Database",
           ifelse(db == "ple", "PLE (Funke et al.)",
                  ifelse(db == "gpd", "GPD (Hawkins et al.)", db)))
  }
  
  # Reverse helpers
  get_method_code <- function(label) ifelse(label == "ASCM MCPanel", "asyn",
                                            ifelse(label == "SCM", "syn", label))
  get_preperiod_code <- function(label) ifelse(label == "15 years", "class1",
                                               ifelse(label == "10 years", "class2", label))
  get_database_code <- function(label) ifelse(label == "Combined Database", "both",
                                              ifelse(label == "PLE (Funke et al.)", "ple",
                                                     ifelse(label == "GPD (Hawkins et al.)", "gpd", label)))
  
  # Dynamic UI for plots
  output$plots_ui <- renderUI({
    if (input$view_mode == "single") {
      fluidRow(column(12, plotOutput("plot_single", height = "750px")))
    } else {
      fluidRow(
        column(6, plotOutput("plot_left", height = "650px")),
        column(6, plotOutput("plot_right", height = "650px"))
      )
    }
  })
  
  # === Single Plot mode ===
  observe({
    dt <- ShinyPopDat$caselevel_dt
    updateSelectInput(session, "casename_single",
                      choices = sort(unique(dt$casename)))
  })
  
  observe({
    req(input$casename_single)
    dt <- ShinyPopDat$caselevel_dt[casename == input$casename_single]
    available_codes <- sort(unique(dt$variable))
    available_labels <- var_labels[var_names %in% available_codes]
    selected_label <- if ("Real GDP per Capita" %in% available_labels) "Real GDP per Capita" else available_labels[1]
    updateSelectInput(session, "variable_single", choices = available_labels, selected = selected_label)
  })
  
  observe({
    req(input$casename_single, input$variable_single)
    code <- label_to_code[input$variable_single]
    dt <- ShinyPopDat$caselevel_dt[casename == input$casename_single & variable == code]
    updateSelectInput(session, "covariates_single", choices = sort(unique(dt$covariates)))
  })
  
  observe({
    req(input$casename_single, input$variable_single, input$covariates_single)
    code <- label_to_code[input$variable_single]
    dt <- ShinyPopDat$caselevel_dt[casename == input$casename_single & variable == code & covariates == input$covariates_single]
    meths <- sort(unique(dt$method))
    updateSelectInput(session, "method_single", choices = label_method(meths))
  })
  
  observe({
    req(input$casename_single, input$variable_single, input$covariates_single, input$method_single)
    code <- label_to_code[input$variable_single]
    raw_method <- get_method_code(input$method_single)
    dt <- ShinyPopDat$caselevel_dt[casename == input$casename_single & variable == code & covariates == input$covariates_single & method == raw_method]
    pres <- sort(unique(dt$preperiod))
    updateSelectInput(session, "preperiod_single", choices = label_preperiod(pres))
  })
  
  observe({
    req(input$casename_single, input$variable_single, input$covariates_single, input$method_single, input$preperiod_single)
    code <- label_to_code[input$variable_single]
    raw_method <- get_method_code(input$method_single)
    raw_pre <- get_preperiod_code(input$preperiod_single)
    dt <- ShinyPopDat$caselevel_dt[casename == input$casename_single & variable == code & covariates == input$covariates_single & method == raw_method & preperiod == raw_pre]
    dbs <- sort(unique(dt$database))
    updateSelectInput(session, "database_single", choices = label_database(dbs))
  })
  
  # === Dual Shared mode ===
  observe({
    dt <- ShinyPopDat$caselevel_dt
    updateSelectInput(session, "casename_shared", choices = sort(unique(dt$casename)))
  })
  
  observe({
    req(input$casename_shared)
    dt <- ShinyPopDat$caselevel_dt[casename == input$casename_shared]
    available_codes <- sort(unique(dt$variable))
    available_labels <- var_labels[var_names %in% available_codes]
    selected_label <- if ("Real GDP per Capita" %in% available_labels) "Real GDP per Capita" else available_labels[1]
    updateSelectInput(session, "variable_shared", choices = available_labels, selected = selected_label)
  })
  
  observe({
    req(input$casename_shared, input$variable_shared, input$casename_right_shared)
    left_code <- label_to_code[input$variable_shared]
    dt <- ShinyPopDat$caselevel_dt
    left_sub <- unique(dt[casename == input$casename_shared & variable == left_code, .(covariates, method, preperiod, database)])
    right_case <- input$casename_right_shared
    right_vars_raw <- unique(dt[casename == right_case]$variable)
    possible_labels <- var_labels[var_names %in% right_vars_raw]
    current_label <- input$variable_right_shared %||% input$variable_shared
    selected <- if (current_label %in% possible_labels) current_label else input$variable_shared
    if (!(selected %in% possible_labels)) selected <- possible_labels[1]
    updateSelectInput(session, "variable_right_shared", choices = possible_labels, selected = selected)
  })
  
  observe({
    req(input$casename_shared, input$variable_shared, input$casename_right_shared, input$variable_right_shared)
    left_code <- label_to_code[input$variable_shared]
    right_code <- label_to_code[input$variable_right_shared]
    dt <- ShinyPopDat$caselevel_dt
    left_sub <- unique(dt[casename == input$casename_shared & variable == left_code, .(covariates, method, preperiod, database)])
    right_sub <- unique(dt[casename == input$casename_right_shared & variable == right_code, .(covariates, method, preperiod, database)])
    merged <- merge(left_sub, right_sub, by = c("covariates", "method", "preperiod", "database"))
    common_cov <- sort(unique(merged$covariates))
    selected <- if (length(common_cov) > 0) common_cov[1] else NULL
    updateSelectInput(session, "covariates_shared", choices = common_cov, selected = selected)
  })
  
  # Remaining shared observers (method, preperiod, database) unchanged except using label_to_code
  observe({
    req(input$casename_shared, input$variable_shared, input$casename_right_shared, input$variable_right_shared, input$covariates_shared)
    left_code <- label_to_code[input$variable_shared]
    right_code <- label_to_code[input$variable_right_shared]
    dt <- ShinyPopDat$caselevel_dt
    left_sub <- unique(dt[casename == input$casename_shared & variable == left_code & covariates == input$covariates_shared, .(method, preperiod, database)])
    right_sub <- unique(dt[casename == input$casename_right_shared & variable == right_code & covariates == input$covariates_shared, .(method, preperiod, database)])
    merged <- merge(left_sub, right_sub, by = c("method", "preperiod", "database"))
    common_meth <- sort(unique(merged$method))
    choices <- label_method(common_meth)
    selected <- choices[1]
    updateSelectInput(session, "method_shared", choices = choices, selected = selected)
  })
  
  observe({
    req(input$casename_shared, input$variable_shared, input$casename_right_shared, input$variable_right_shared, input$covariates_shared, input$method_shared)
    left_code <- label_to_code[input$variable_shared]
    right_code <- label_to_code[input$variable_right_shared]
    raw_method <- get_method_code(input$method_shared)
    dt <- ShinyPopDat$caselevel_dt
    left_sub <- unique(dt[casename == input$casename_shared & variable == left_code & covariates == input$covariates_shared & method == raw_method, .(preperiod, database)])
    right_sub <- unique(dt[casename == input$casename_right_shared & variable == right_code & covariates == input$covariates_shared & method == raw_method, .(preperiod, database)])
    merged <- merge(left_sub, right_sub, by = c("preperiod", "database"))
    common_pre <- sort(unique(merged$preperiod))
    choices <- label_preperiod(common_pre)
    selected <- choices[1]
    updateSelectInput(session, "preperiod_shared", choices = choices, selected = selected)
  })
  
  observe({
    req(input$casename_shared, input$variable_shared, input$casename_right_shared, input$variable_right_shared, input$covariates_shared, input$method_shared, input$preperiod_shared)
    left_code <- label_to_code[input$variable_shared]
    right_code <- label_to_code[input$variable_right_shared]
    raw_method <- get_method_code(input$method_shared)
    raw_pre <- get_preperiod_code(input$preperiod_shared)
    dt <- ShinyPopDat$caselevel_dt
    left_sub <- unique(dt[casename == input$casename_shared & variable == left_code & covariates == input$covariates_shared & method == raw_method & preperiod == raw_pre, .(database)])
    right_sub <- unique(dt[casename == input$casename_right_shared & variable == right_code & covariates == input$covariates_shared & method == raw_method & preperiod == raw_pre, .(database)])
    merged <- merge(left_sub, right_sub, by = "database")
    common_db <- sort(unique(merged$database))
    choices <- label_database(common_db)
    selected <- choices[1]
    updateSelectInput(session, "database_shared", choices = choices, selected = selected)
  })
  
  # === Dual Separate mode (Left) ===
  observe({
    dt <- ShinyPopDat$caselevel_dt
    updateSelectInput(session, "casename_left", choices = sort(unique(dt$casename)))
  })
  
  observe({
    req(input$casename_left)
    dt <- ShinyPopDat$caselevel_dt[casename == input$casename_left]
    available_codes <- sort(unique(dt$variable))
    available_labels <- var_labels[var_names %in% available_codes]
    selected_label <- if ("Real GDP per Capita" %in% available_labels) "Real GDP per Capita" else available_labels[1]
    updateSelectInput(session, "variable_left", choices = available_labels, selected = selected_label)
  })
  
  observe({
    req(input$casename_left, input$variable_left)
    code <- label_to_code[input$variable_left]
    dt <- ShinyPopDat$caselevel_dt[casename == input$casename_left & variable == code]
    updateSelectInput(session, "covariates_left", choices = sort(unique(dt$covariates)))
  })
  
  observe({
    req(input$casename_left, input$variable_left, input$covariates_left)
    code <- label_to_code[input$variable_left]
    dt <- ShinyPopDat$caselevel_dt[casename == input$casename_left & variable == code & covariates == input$covariates_left]
    meths <- sort(unique(dt$method))
    updateSelectInput(session, "method_left", choices = label_method(meths))
  })
  
  observe({
    req(input$casename_left, input$variable_left, input$covariates_left, input$method_left)
    code <- label_to_code[input$variable_left]
    raw_method <- get_method_code(input$method_left)
    dt <- ShinyPopDat$caselevel_dt[casename == input$casename_left & variable == code & covariates == input$covariates_left & method == raw_method]
    pres <- sort(unique(dt$preperiod))
    updateSelectInput(session, "preperiod_left", choices = label_preperiod(pres))
  })
  
  observe({
    req(input$casename_left, input$variable_left, input$covariates_left, input$method_left, input$preperiod_left)
    code <- label_to_code[input$variable_left]
    raw_method <- get_method_code(input$method_left)
    raw_pre <- get_preperiod_code(input$preperiod_left)
    dt <- ShinyPopDat$caselevel_dt[casename == input$casename_left & variable == code & covariates == input$covariates_left & method == raw_method & preperiod == raw_pre]
    dbs <- sort(unique(dt$database))
    updateSelectInput(session, "database_left", choices = label_database(dbs))
  })
  
  # === Dual Separate mode (Right) ===
  observe({
    dt <- ShinyPopDat$caselevel_dt
    updateSelectInput(session, "casename_right", choices = sort(unique(dt$casename)))
  })
  
  observe({
    req(input$casename_right)
    dt <- ShinyPopDat$caselevel_dt[casename == input$casename_right]
    available_codes <- sort(unique(dt$variable))
    available_labels <- var_labels[var_names %in% available_codes]
    selected_label <- if ("Disposable Gini (SWIID)" %in% available_labels) "Disposable Gini (SWIID)" else available_labels[1]
    updateSelectInput(session, "variable_right", choices = available_labels, selected = selected_label)
  })
  
  observe({
    req(input$casename_right, input$variable_right)
    code <- label_to_code[input$variable_right]
    dt <- ShinyPopDat$caselevel_dt[casename == input$casename_right & variable == code]
    updateSelectInput(session, "covariates_right", choices = sort(unique(dt$covariates)))
  })
  
  observe({
    req(input$casename_right, input$variable_right, input$covariates_right)
    code <- label_to_code[input$variable_right]
    dt <- ShinyPopDat$caselevel_dt[casename == input$casename_right & variable == code & covariates == input$covariates_right]
    meths <- sort(unique(dt$method))
    updateSelectInput(session, "method_right", choices = label_method(meths))
  })
  
  observe({
    req(input$casename_right, input$variable_right, input$covariates_right, input$method_right)
    code <- label_to_code[input$variable_right]
    raw_method <- get_method_code(input$method_right)
    dt <- ShinyPopDat$caselevel_dt[casename == input$casename_right & variable == code & covariates == input$covariates_right & method == raw_method]
    pres <- sort(unique(dt$preperiod))
    updateSelectInput(session, "preperiod_right", choices = label_preperiod(pres))
  })
  
  observe({
    req(input$casename_right, input$variable_right, input$covariates_right, input$method_right, input$preperiod_right)
    code <- label_to_code[input$variable_right]
    raw_method <- get_method_code(input$method_right)
    raw_pre <- get_preperiod_code(input$preperiod_right)
    dt <- ShinyPopDat$caselevel_dt[casename == input$casename_right & variable == code & covariates == input$covariates_right & method == raw_method & preperiod == raw_pre]
    dbs <- sort(unique(dt$database))
    updateSelectInput(session, "database_right", choices = label_database(dbs))
  })
  
  # === Plot rendering ===
  output$plot_single <- renderPlot({
    req(input$view_mode == "single")
    plot_effects(
      dt = ShinyPopDat$caselevel_dt,
      in_casename = input$casename_single,
      in_variable = label_to_code[input$variable_single],
      in_covariates = input$covariates_single,
      in_method = get_method_code(input$method_single),
      in_preperiod = get_preperiod_code(input$preperiod_single),
      in_database = get_database_code(input$database_single),
      baseval = ShinyPopDat$basevalues_dt,
      set_backtracing_3yr = input$set_backtracing_3yr_single,
      set_backtracing_5yr = input$set_backtracing_5yr_single,
      set_case_permutation = input$set_case_permutation_single,
      set_loo_permutation = input$set_loo_permutation_single,
      set_ci = input$set_ci_single,
      series = input$series_single
    )
  })
  
  output$plot_left <- renderPlot({
    req(input$view_mode == "dual")
    if (input$control_mode == "Shared") {
      plot_effects(
        dt = ShinyPopDat$caselevel_dt,
        in_casename = input$casename_shared,
        in_variable = label_to_code[input$variable_shared],
        in_covariates = input$covariates_shared,
        in_method = get_method_code(input$method_shared),
        in_preperiod = get_preperiod_code(input$preperiod_shared),
        in_database = get_database_code(input$database_shared),
        baseval = ShinyPopDat$basevalues_dt,
        set_backtracing_3yr = input$set_backtracing_3yr_shared,
        set_backtracing_5yr = input$set_backtracing_5yr_shared,
        set_case_permutation = input$set_case_permutation_shared,
        set_loo_permutation = input$set_loo_permutation_shared,
        set_ci = input$set_ci_shared,
        series = input$series_shared
      )
    } else {
      plot_effects(
        dt = ShinyPopDat$caselevel_dt,
        in_casename = input$casename_left,
        in_variable = label_to_code[input$variable_left],
        in_covariates = input$covariates_left,
        in_method = get_method_code(input$method_left),
        in_preperiod = get_preperiod_code(input$preperiod_left),
        in_database = get_database_code(input$database_left),
        baseval = ShinyPopDat$basevalues_dt,
        set_backtracing_3yr = input$set_backtracing_3yr_left,
        set_backtracing_5yr = input$set_backtracing_5yr_left,
        set_case_permutation = input$set_case_permutation_left,
        set_loo_permutation = input$set_loo_permutation_left,
        set_ci = input$set_ci_left,
        series = input$series_left
      )
    }
  })
  
  output$plot_right <- renderPlot({
    req(input$view_mode == "dual")
    if (input$control_mode == "Shared") {
      plot_effects(
        dt = ShinyPopDat$caselevel_dt,
        in_casename = input$casename_right_shared,
        in_variable = label_to_code[input$variable_right_shared],
        in_covariates = input$covariates_shared,
        in_method = get_method_code(input$method_shared),
        in_preperiod = get_preperiod_code(input$preperiod_shared),
        in_database = get_database_code(input$database_shared),
        baseval = ShinyPopDat$basevalues_dt,
        set_backtracing_3yr = input$set_backtracing_3yr_shared,
        set_backtracing_5yr = input$set_backtracing_5yr_shared,
        set_case_permutation = input$set_case_permutation_shared,
        set_loo_permutation = input$set_loo_permutation_shared,
        set_ci = input$set_ci_shared,
        series = input$series_shared
      )
    } else {
      plot_effects(
        dt = ShinyPopDat$caselevel_dt,
        in_casename = input$casename_right,
        in_variable = label_to_code[input$variable_right],
        in_covariates = input$covariates_right,
        in_method = get_method_code(input$method_right),
        in_preperiod = get_preperiod_code(input$preperiod_right),
        in_database = get_database_code(input$database_right),
        baseval = ShinyPopDat$basevalues_dt,
        set_backtracing_3yr = input$set_backtracing_3yr_right,
        set_backtracing_5yr = input$set_backtracing_5yr_right,
        set_case_permutation = input$set_case_permutation_right,
        set_loo_permutation = input$set_loo_permutation_right,
        set_ci = input$set_ci_right,
        series = input$series_right
      )
    }
  })
}

shinyApp(ui = ui, server = server)
