# =============================================================================
# Staggered DiD Simulation: AI Rollout Panel
# Discrete DGP + Continuous Extension (CGBS)
# =============================================================================
library(tidyverse)
library(fixest)
library(bacondecomp)
library(gridExtra)
library(cowplot)

COL_SA   <- "#1B6CA8"   
COL_CGBS <- "#E84E0F" 
COL_HS   <- "#E84E0F"
COL_MS   <- "#1B6CA8"

# =============================================================================
# PART I — DISCRETE DGP
# =============================================================================
set.seed(0001)
N_REGIONS <- 1000
N_YEARS   <- 20
BASE_POP  <- 10000

region_id  <- 1:N_REGIONS
inf_capcty <- runif(N_REGIONS, 0, 1)

share_H_raw <- pmax(rnorm(N_REGIONS, mean = 20 + 15 * inf_capcty, sd = 5), 1)
share_M_raw <- rnorm(N_REGIONS, mean = 30, sd = 5)
share_L_raw <- pmax(rnorm(N_REGIONS, mean = 50 - 15 * inf_capcty, sd = 5), 1)
total_raw   <- share_H_raw + share_M_raw + share_L_raw

L_H0 <- round(BASE_POP * share_H_raw / total_raw)
L_M0 <- round(BASE_POP * share_M_raw / total_raw)
L_L0 <- BASE_POP - L_H0 - L_M0
K_H0 <- 1.0 * L_H0
K_M0 <- 1.0 * L_M0
K_L0 <- 1.0 * L_L0

# Treatment year: depends on infrastructural capacity; below 0.5 = control
treat_year <- ifelse(inf_capcty >= 0.9, 5,
              ifelse(inf_capcty >= 0.8, 7,
              ifelse(inf_capcty >= 0.7, 9,
              ifelse(inf_capcty >= 0.6, 11,
              ifelse(inf_capcty >= 0.5, 13, 0)))))

cohort_label <- ifelse(treat_year == 5,  "Cohort_05",
                ifelse(treat_year == 7,  "Cohort_07",
                ifelse(treat_year == 9,  "Cohort_09",
                ifelse(treat_year == 11, "Cohort_11",
                ifelse(treat_year == 13, "Cohort_13", "Control")))))

regions <- data.frame(region_id, inf_capcty, treat_year, cohort_label,
                      BASE_POP, L_H0, L_M0, L_L0, K_H0, K_M0, K_L0)

panel <- expand.grid(year = 1:N_YEARS, region_id = 1:N_REGIONS) %>%
  merge(regions, by = "region_id") %>%
  arrange(region_id, year)

# Treatment indicators
panel$treated_post <- as.integer(panel$treat_year > 0 & panel$year >= panel$treat_year)
panel$years_since  <- ifelse(panel$treat_year > 0, panel$year - panel$treat_year, -100)
panel$phase        <- ifelse(panel$years_since < 0, 0,
                      ifelse(panel$years_since >= 10, 1, panel$years_since / 10))

# Production parameters
tfp_mult          <- (1.02)^(panel$year - 1)
panel$alpha_H_eff <- 0.60 - 0.25 * panel$phase * exp(rnorm(nrow(panel), 0, 0.25))     # labor share: 0.60 -> 0.35 with AI
panel$A_H_eff     <- (2.0 + 2.0 * panel$inf_capcty) * tfp_mult
panel$A_M_eff     <- 2.0 * tfp_mult * (1 + 0.80 * panel$phase)
panel$A_L_eff     <- 1.0 * tfp_mult

# Labor flows: AI displaces 25% of HS over 10 years
HS_displaced <- panel$phase * 0.25 * panel$L_H0 * exp(rnorm(nrow(panel), 0, 0.05))
panel$L_H    <- panel$L_H0 - HS_displaced
panel$L_M    <- panel$L_M0 + 0.2 * HS_displaced
panel$L_L    <- panel$L_L0 + 0.6 * HS_displaced
panel$L_total  <- panel$L_H + panel$L_M + panel$L_L
 
# Cobb-Douglas output and wages
panel$Y_H     <- panel$A_H_eff * panel$L_H^panel$alpha_H_eff * panel$K_H0^(1 - panel$alpha_H_eff)
panel$Y_M     <- panel$A_M_eff * panel$L_M^0.50              * panel$K_M0^0.50
panel$Y_L     <- panel$A_L_eff * panel$L_L^0.40              * panel$K_L0^0.60
panel$Y_total <- panel$Y_H + panel$Y_M + panel$Y_L

panel$Wage_H <- panel$alpha_H_eff * (panel$Y_H / panel$L_H) * exp(rnorm(nrow(panel), 0, 0.025))
panel$Wage_M <- 0.50              * (panel$Y_M / panel$L_M) * exp(rnorm(nrow(panel), 0, 0.025))
panel$Wage_L <- 0.40              * (panel$Y_L / panel$L_L) * exp(rnorm(nrow(panel), 0, 0.025))
panel$Wage_avg <- (panel$Wage_H*panel$L_H + panel$Wage_M*panel$L_M + panel$Wage_L*panel$L_L)/
                  (panel$L_H + panel$L_M + panel$L_L)

panel$Labor_Share   <- (panel$Wage_H*panel$L_H + panel$Wage_M*panel$L_M + panel$Wage_L*panel$L_L) / panel$Y_total
panel$Capital_Share <- 1 - panel$Labor_Share

panel <- panel %>% mutate(
  ln_Wage_H = log(Wage_H), ln_Wage_M = log(Wage_M), ln_Wage_L = log(Wage_L),ln_Wage_avg = log(Wage_avg),
  ln_L_H    = log(L_H),    ln_L_M    = log(L_M),    ln_L_L    = log(L_L),   ln_L_total = log(L_total),
  ln_Y_H    = log(Y_H),    ln_Y_M    = log(Y_M),    ln_Y_L    = log(Y_L),   ln_Y_total = log(Y_total)
)

# =============================================================================
# PART II — DISCRETE ESTIMATION
# =============================================================================
vars_disc <- c("ln_Wage_H", "ln_Wage_M", "ln_Wage_L", "ln_Wage_avg",
               "ln_L_H",    "ln_L_M",    "ln_L_L",  "ln_L_total",
               "ln_Y_H",    "ln_Y_M",    "ln_Y_L","ln_Y_total",
               "Labor_Share")

mod_ols_disc  <- feols(.[vars_disc] ~ treated_post,                               data = panel)
mod_twfe_disc <- feols(.[vars_disc] ~ treated_post | region_id + year,            data = panel)
mod_sa_disc   <- feols(.[vars_disc] ~ sunab(treat_year, year) | region_id + year, data = panel)

bacon_disc    <- bacon(ln_L_H ~ treated_post, data = panel,
                       id_var = "region_id", time_var = "year")

# =============================================================================
# PART III — DISCRETE PLOTS and TABLE
# =============================================================================

# Extract SA dynamic coefficients
extract_sa <- function(sa_model, label) {
  df <- as.data.frame(coeftable(sa_model))
  names(df)[1:2] <- c("Estimate", "StdError")
  df$term     <- rownames(df)
  df          <- df[grepl("year", df$term) & !grepl("Intercept", df$term), ]
  df$rel_time <- as.numeric(gsub(".*[^0-9-](-?[0-9]+).*", "\\1", df$term))
  ref         <- df[1, ]; ref$Estimate <- 0; ref$StdError <- 0; ref$rel_time <- -1
  df          <- rbind(df, ref)[order(rbind(df, ref)$rel_time), ]
  df$ci_lower  <- df$Estimate - 1.96 * df$StdError
  df$ci_upper  <- df$Estimate + 1.96 * df$StdError
  df$estimator <- label
  df
}

plot_sa <- function(sa_mod, ols_mod, twfe_mod, title, subtitle, ylab) {
  df        <- extract_sa(sa_mod, "Sun & Abraham")
  ols_coef  <- coef(ols_mod)["treated_post"]
  twfe_coef <- coef(twfe_mod)["treated_post"]
  ggplot(df, aes(x = rel_time, y = Estimate)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = COL_SA, alpha = 0.15) +
    geom_hline(yintercept = 0,         color = "black",   linewidth = 0.4) +
    geom_vline(xintercept = -0.5,      linetype = "dashed",  color = "black", linewidth = 0.5) +
    geom_hline(yintercept = ols_coef,  linetype = "dashed",  color = "grey50", linewidth = 0.7) +
    geom_hline(yintercept = twfe_coef, linetype = "dotted",  color = "grey50", linewidth = 0.7) +
    annotate("text", x = max(df$rel_time) - 1, y = ols_coef,
             label = paste("OLS =", round(ols_coef, 3)),
             color = "grey30", vjust = -0.5, fontface = "bold", size = 5) +   # was 3
    annotate("text", x = max(df$rel_time) - 1, y = twfe_coef,
             label = paste("TWFE =", round(twfe_coef, 3)),
             color = "grey30", vjust = 1.5, fontface = "bold", size = 5) +    # was 3
    geom_line(color = COL_SA, linewidth = 0.9) +
    geom_point(shape = 21, size = 2.5, fill = COL_SA, color = "black") +
    theme_classic(base_size = 20) +                                            # was default 11
    labs(title = title, subtitle = subtitle, x = "Years Relative to AI Adoption", y = ylab) +
    theme(plot.title    = element_text(face = "bold", size = 18),              # was 13
          plot.subtitle = element_text(size = 14),
          legend.position = "none")
}

# Plot 1: HS Employment displacement
print(plot_sa(mod_sa_disc[[5]], mod_ols_disc[[5]], mod_twfe_disc[[5]],
              "AI Impact on High-Skill Employment",
              "Sun & Abraham ATT with OLS/TWFE references",
              "Log High-Skill Employment"))

# Plot 2: Low-Skill Wages (supply spillover)
print(plot_sa(mod_sa_disc[[3]], mod_ols_disc[[3]], mod_twfe_disc[[3]],
              "AI Impact on Low-Skill Wages",
              "Sun & Abraham ATT with OLS/TWFE references",
              "Log Low-Skill Wage"))

# Plot 3: Observed wage paths by cohort (TWFE forbidden comparisons)
plot_cohorts <- panel %>%
  filter(cohort_label %in% c("Cohort_05", "Cohort_11", "Control")) %>%
  group_by(year, cohort_label) %>%
  summarise(Wage = mean(ln_Wage_H), .groups = "drop")

print(
  ggplot(plot_cohorts, aes(x = year, y = Wage, color = cohort_label)) +
    geom_line(linewidth = 1.4) +
    geom_vline(xintercept = 5,  linetype = "solid", color = "black", linewidth = 0.4) +
    geom_vline(xintercept = 11, linetype = "solid", color = "black", linewidth = 0.4) +
    annotate("text", x = 3,  y = max(plot_cohorts$Wage)*0.99, label = "PRE",       fontface = "bold", size = 3.5) +
    annotate("text", x = 8,  y = max(plot_cohorts$Wage)*0.99, label = "MID(5,11)", fontface = "bold", size = 3.5) +
    annotate("text", x = 16, y = max(plot_cohorts$Wage)*0.99, label = "POST",      fontface = "bold", size = 3.5) +
    scale_color_manual(values = c("Cohort_05" = COL_HS, "Cohort_11" = COL_MS, "Control" = "black"),
                       labels = c("Early Adopters (Yr 5)", "Late Adopters (Yr 11)", "Never Treated")) +
    theme_classic(base_size = 20) +
    labs(title    = "The TWFE Trap: Forbidden Comparisons in Staggered Adoption",
         subtitle = "Early adopters serve as controls for late adopters mid-window",
         x = "Simulation Year", y = "Log High-Skill Wage", color = NULL) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 13))
)


#-----------------------Final Discrete Table------------------------------------
sa_att_se <- function(mod) {
  s <- summary(mod, agg = "att")$coeftable
  c(att = s[1, 1], se = s[1, 2])
}

fmt <- function(est, se) sprintf("%.4f\n(%.4f)", est, se)

disc_tbl_vars   <- c("ln_L_H", "ln_L_total", "ln_Wage_H", "ln_Wage_M","ln_Wage_L",
                     "ln_Y_H","ln_Y_total","Labor_Share")
disc_tbl_labels <- c("High-Skill Employment", "Total Employment", "High-Skill Wage", "Mid-Skill Wage", "Low-Skill Wage",
                      "Output - High Skill","Output Total","Labor Share")

gb_disc_att <- weighted.mean(bacon_disc$estimate, bacon_disc$weight)

disc_table <- lapply(seq_along(disc_tbl_vars), function(i) {
  v   <- disc_tbl_vars[i]
  idx <- which(vars_disc == v)
  sa  <- sa_att_se(mod_sa_disc[[idx]])
  data.frame(
    Variable   = disc_tbl_labels[i],
    OLS        = fmt(coef(mod_ols_disc[[idx]])["treated_post"],  se(mod_ols_disc[[idx]])["treated_post"]),
    TWFE       = fmt(coef(mod_twfe_disc[[idx]])["treated_post"], se(mod_twfe_disc[[idx]])["treated_post"]),
    SA_ATT     = fmt(sa["att"], sa["se"])
  )
}) %>% bind_rows()

cat("\n", strrep("=", 72), "\n")
cat("  SECTION A — DISCRETE AI IMPLEMENTATION\n")
cat("  Standard errors in parentheses | GB ATT = Goodman-Bacon (ln_L_H only)\n")
cat(strrep("=", 72), "\n")
print(disc_table, row.names = FALSE)


# =============================================================================
# PART IV — CONTINUOUS DGP
# =============================================================================
set.seed(0002)

investment_annual  <- ifelse(regions$treat_year > 0,
                              100000 * (regions$inf_capcty * 2) + rnorm(N_REGIONS, 0, 5000), 0)
regions$investment <- investment_annual * 10
panel$investment   <- regions$investment[match(panel$region_id, regions$region_id)]
panel$invest_post  <- panel$investment * panel$treated_post
panel$years_active <- ifelse(panel$treat_year == 0, 0, pmin(pmax(panel$years_since, 0), 10))

# Capital: accumulates only in HS sector from AI investment
panel$K_H_cont <- panel$K_H0 + (investment_annual[match(panel$region_id, regions$region_id)] / 1000) * panel$years_active
panel$K_M_cont <- panel$K_M0
panel$K_L_cont <- panel$K_L0

# Displacement dependent on Investment (mean investment = -25pp in 10 years) + noise
mean_investment <- mean(100000 * (regions$inf_capcty[regions$treat_year > 0] * 2) * 10)
panel$disp_frac_cont <- (panel$investment / mean_investment) * (panel$years_active / 10) * 0.25
panel$L_H_cont       <- panel$L_H0 - panel$disp_frac_cont * panel$L_H0  * exp(rnorm(nrow(panel), 0, 0.1))
panel$L_M_cont       <- panel$L_M0 + 0.2 * panel$disp_frac_cont * panel$L_H0
panel$L_L_cont       <- panel$L_L0 + 0.6 * panel$disp_frac_cont * panel$L_H0
panel$L_total_cont   <- panel$L_H_cont + panel$L_M_cont + panel$L_L_cont

# Mid-skill TFP (continuous)
panel$A_M_eff_cont   <- 2.0 * tfp_mult * (1 + 0.80 * (panel$investment / mean_investment) * (panel$years_active / 10))

# Cobb-Douglas output and wages (continuous)
panel$Y_H_cont     <- panel$A_H_eff      * panel$L_H_cont^panel$alpha_H_eff * panel$K_H_cont^(1 - panel$alpha_H_eff)
panel$Y_M_cont     <- panel$A_M_eff_cont * panel$L_M_cont^0.50              * panel$K_M_cont^0.50
panel$Y_L_cont     <- panel$A_L_eff      * panel$L_L_cont^0.40              * panel$K_L_cont^0.60
panel$Y_total_cont <- panel$Y_H_cont + panel$Y_M_cont + panel$Y_L_cont

# Wages (continuous)
panel$Wage_H_cont <- panel$alpha_H_eff * (panel$Y_H_cont / panel$L_H_cont) * exp(rnorm(nrow(panel), 0, 0.025))
panel$Wage_M_cont <- 0.50              * (panel$Y_M_cont / panel$L_M_cont) * exp(rnorm(nrow(panel), 0, 0.025))
panel$Wage_L_cont <- 0.40              * (panel$Y_L_cont / panel$L_L_cont) * exp(rnorm(nrow(panel), 0, 0.025))

panel$Wage_avg_cont     <- (panel$Wage_H_cont*panel$L_H_cont + panel$Wage_M_cont*panel$L_M_cont +
                             panel$Wage_L_cont*panel$L_L_cont) /
                           (panel$L_H_cont + panel$L_M_cont + panel$L_L_cont)
# Labor and Capital share (continuous)
panel$Labor_Share_cont  <- (panel$Wage_H_cont*panel$L_H_cont + panel$Wage_M_cont*panel$L_M_cont +
                             panel$Wage_L_cont*panel$L_L_cont) / panel$Y_total_cont
panel$Capital_Share_cont <- 1 - panel$Labor_Share_cont

# Log-transformation of continuous variables
panel <- panel %>% mutate(
  ln_L_H_cont      = log(L_H_cont),    ln_L_M_cont      = log(L_M_cont),    ln_L_L_cont      = log(L_L_cont),
  ln_L_total_cont  = log(L_total_cont),
  ln_Wage_H_cont   = log(Wage_H_cont), ln_Wage_M_cont   = log(Wage_M_cont), ln_Wage_L_cont   = log(Wage_L_cont),
  ln_Wage_avg_cont = log(Wage_avg_cont),
  ln_Y_H_cont      = log(Y_H_cont),    ln_Y_M_cont      = log(Y_M_cont),    ln_Y_L_cont      = log(Y_L_cont),
  ln_Y_total_cont  = log(Y_total_cont),
)


# =============================================================================
# PART V — CONTINUOUS ESTIMATION
# =============================================================================
vars_cont <- c("ln_L_H_cont",  "ln_L_M_cont",    "ln_L_L_cont", "ln_L_total_cont",
              "ln_Wage_H_cont",  "ln_Wage_M_cont", "ln_Wage_L_cont",
               "ln_Y_H_cont",  "ln_Y_M_cont",    "ln_Y_L_cont", "ln_Y_total_cont",
               "Capital_Share_cont", "Labor_Share_cont")

mod_ols_cont  <- feols(.[vars_cont] ~ invest_post,                               data = panel)
mod_twfe_cont <- feols(.[vars_cont] ~ invest_post | region_id + year,            data = panel)
mod_sa_cont   <- feols(.[vars_cont] ~ sunab(treat_year, year) | region_id + year, data = panel)

mean_dose  <- mean(panel$investment[panel$treat_year > 0 & panel$treated_post == 1])
bacon_cont <- bacon(ln_L_H_cont ~ treated_post, data = panel,
                   id_var = "region_id", time_var = "year")

# Function to run CGBS for one variable, return ATT, ATT_se, and dynamic series

run_cgbs <- function(varname) {
  # Part i. extracting baseline year
  base_yr <- panel %>% filter(year == 1) %>% select(region_id, yb = all_of(varname))
  
  #Part ii. estimating pre-treatment trends
  pre_sl  <- panel %>%
    filter(treat_year == 0 | year < treat_year) %>%
    group_by(region_id) %>%
    summarise(pre_slope = coef(lm(get(varname) ~ year))[["year"]], .groups = "drop")
  
  # Part iii. computing the y change
  dat <- panel %>%
    left_join(pre_sl,  by = "region_id") %>%
    left_join(base_yr, by = "region_id") %>%
    mutate(y_change = get(varname) - yb)
  dyn <- lapply(sort(unique(dat$year)), function(t) {
    if (t == 1) return(NULL)
    df_t <- dat %>% filter(year == t, !is.na(y_change), !is.na(investment))
    if (nrow(df_t) < 10) return(NULL)
    
    # Part iv. estimating y_change on the investment dose and the pre-slope captured trends
    m <- tryCatch(lm(y_change ~ investment + pre_slope, data = df_t), error = function(e) NULL)
    if (is.null(m)) return(NULL)
    data.frame(event_time = t - 5,
               ATT_slope  = coef(m)[["investment"]],
               SE         = sqrt(diag(vcov(m)))[["investment"]])
  }) %>% bind_rows()
  
  # Compute ATT and ATT_se for the final table
  post <- dyn$event_time >= 0
  list(dyn    = dyn,
       att    = mean(dyn$ATT_slope[post]) * mean_dose,
       att_se = sqrt(mean(dyn$SE[post]^2)) * mean_dose)
}

# Run CGBS once per variable and results
cgbs_results <- setNames(lapply(vars_cont, run_cgbs), vars_cont)

# =============================================================================
# PART VI — CONTINUOUS PLOTS
# =============================================================================

# Plot 5: High-skill factor share evolution (stacked area)
area_data <- panel %>%
  filter(treat_year > 0) %>%
  group_by(years_since) %>%
  summarise(Labor_Share_H   = mean((Wage_H_cont * L_H_cont) / Y_H_cont),
            Capital_Share_H = mean(1 - (Wage_H_cont * L_H_cont) / Y_H_cont),
            .groups = "drop") %>%
  filter(years_since >= -5 & years_since <= 10) %>%
  pivot_longer(c(Labor_Share_H, Capital_Share_H), names_to = "Share_Type", values_to = "Pct")

print(
  ggplot(area_data, aes(x = years_since, y = Pct, fill = Share_Type)) +
    geom_area(alpha = 0.85, color = "white", linewidth = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    scale_fill_manual(values = c("Capital_Share_H" = COL_CGBS, "Labor_Share_H" = COL_SA),
                      labels = c("Capital Share", "Labor Share")) +
    scale_y_continuous(labels = scales::percent_format()) +
    theme_classic(base_size = 20) +
    labs(title    = "Evolution of Factor Shares in High-Skill Production",
         subtitle = "Continuous DGP: $100k/yr investment accumulates capital over 10 years",
         x = "Years Relative to AI Adoption", y = "Share of Output") +
    theme(legend.position = "bottom", legend.title = element_blank(),
          plot.title = element_text(face = "bold", size = 13))
)

# Plot 6: SA vs CGBS event-study overlay (ln_L_H_cont)
cgbs_lh    <- cgbs_results[["ln_L_H_cont"]]
df_sa_cont <- extract_sa(mod_sa_cont[[1]], "Sun & Abraham")
cgbs_es <- data.frame(
  rel_time  = cgbs_lh$dyn$event_time,
  Estimate  = cgbs_lh$dyn$ATT_slope * mean_dose,
  ci_lower  = (cgbs_lh$dyn$ATT_slope - 1.96 * cgbs_lh$dyn$SE) * mean_dose,
  ci_upper  = (cgbs_lh$dyn$ATT_slope + 1.96 * cgbs_lh$dyn$SE) * mean_dose,
  estimator = "CGBS"
)
df_overlay <- rbind(df_sa_cont[, c("rel_time","Estimate","ci_lower","ci_upper","estimator")], cgbs_es)
ols_ref    <- coef(mod_ols_cont[[1]])["invest_post"]  * mean_dose
twfe_ref   <- coef(mod_twfe_cont[[1]])["invest_post"] * mean_dose

print(
  ggplot(df_overlay, aes(x = rel_time, y = Estimate, color = estimator, fill = estimator)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.15, color = NA) +
    geom_hline(yintercept = 0,        color = "black",  linewidth = 0.4) +
    geom_vline(xintercept = -0.5,     linetype = "dashed",  color = "black", linewidth = 0.5) +
    geom_hline(yintercept = ols_ref,  linetype = "dashed",  color = "grey50", linewidth = 0.7) +
    geom_hline(yintercept = twfe_ref, linetype = "dotted",  color = "grey50", linewidth = 0.7) +
    annotate("text", x = max(df_overlay$rel_time) - 1, y = ols_ref,
             label = paste("OLS =",  round(ols_ref,  3)), color = "grey30", vjust = -0.5, fontface = "bold", size = 3) +
    annotate("text", x = max(df_overlay$rel_time) - 1, y = twfe_ref,
             label = paste("TWFE =", round(twfe_ref, 3)), color = "grey30", vjust = 1.5,  fontface = "bold", size = 3) +
    geom_line(linewidth = 0.9) +
    geom_point(shape = 21, size = 2.5, color = "black") +
    scale_color_manual(values = c("Sun & Abraham" = COL_SA, "CGBS" = COL_CGBS)) +
    scale_fill_manual( values = c("Sun & Abraham" = COL_SA, "CGBS" = COL_CGBS)) +
    theme_classic(base_size = 20) +
    labs(title    = "Continuous Treatment Event Study: SA vs CGBS",
         subtitle = "Dose rescaled to mean cumulative AI investment (~$1M)",
         x = "Years Relative to AI Adoption",
         y = "ATT: Log High-Skill Employment (at mean dose)",
         color = NULL, fill = NULL) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 13))
)

# =============================================================================
# PART VII — COMPREHENSIVE GRID (4 rows x 4 cols)
# Rows: Employment, Wage, Capital Share, Output
# Cols: High Skill, Mid Skill, Low Skill, Total
# Each cell: SA (blue) + CGBS (orange) + OLS/TWFE reference lines
# =============================================================================
grid_vars <- list(
  Employment    = c(HS = "ln_L_H_cont",       MS = "ln_L_M_cont",       LS = "ln_L_L_cont",       Total = "ln_L_total_cont"),
  Wage          = c(HS = "ln_Wage_H_cont",     MS = "ln_Wage_M_cont",     LS = "ln_Wage_L_cont",     Total = "ln_Wage_avg_cont"),
  Output        = c(HS = "ln_Y_H_cont",        MS = "ln_Y_M_cont",        LS = "ln_Y_L_cont",        Total = "ln_Y_total_cont")
)

vars_grid_all <- unique(c(vars_cont,
                          "ln_Wage_H_cont", "ln_Wage_L_cont", "ln_Wage_avg_cont",
                          "ln_L_M_cont",    "ln_L_L_cont"))

mod_ols_grid  <- feols(.[vars_grid_all] ~ invest_post,                               data = panel)
mod_twfe_grid <- feols(.[vars_grid_all] ~ invest_post | region_id + year,            data = panel)
mod_sa_grid   <- feols(.[vars_grid_all] ~ sunab(treat_year, year) | region_id + year, data = panel)
cgbs_grid     <- setNames(lapply(vars_grid_all, run_cgbs), vars_grid_all)

make_cell <- function(varname, row_label, col_label) {
  idx <- which(vars_grid_all == varname)
  if (length(idx) == 0) return(ggplot() + theme_void())
  cg      <- cgbs_grid[[varname]]
  df_sa   <- extract_sa(mod_sa_grid[[idx]], "Sun & Abraham")
  cgbs_df <- data.frame(
    rel_time  = cg$dyn$event_time,
    Estimate  = cg$dyn$ATT_slope * mean_dose,
    ci_lower  = (cg$dyn$ATT_slope - 1.96 * cg$dyn$SE) * mean_dose,
    ci_upper  = (cg$dyn$ATT_slope + 1.96 * cg$dyn$SE) * mean_dose,
    estimator = "CGBS"
  )
  df_combined <- rbind(df_sa[, c("rel_time","Estimate","ci_lower","ci_upper","estimator")], cgbs_df)
  ols_c  <- coef(mod_ols_grid[[idx]])["invest_post"]  * mean_dose
  twfe_c <- coef(mod_twfe_grid[[idx]])["invest_post"] * mean_dose

  ggplot(df_combined, aes(x = rel_time, y = Estimate, color = estimator, fill = estimator)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.12, color = NA) +
    geom_hline(yintercept = 0,      color = "black",  linewidth = 0.3) +
    geom_vline(xintercept = -0.5,   linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_hline(yintercept = ols_c,  linetype = "dashed", color = "grey60", linewidth = 0.5) +
    geom_hline(yintercept = twfe_c, linetype = "dotted", color = "grey60", linewidth = 0.5) +
    geom_line(linewidth = 0.7) +
    geom_point(shape = 21, size = 1.8, color = "black") +
    scale_color_manual(values = c("Sun & Abraham" = COL_SA, "CGBS" = COL_CGBS)) +
    scale_fill_manual( values = c("Sun & Abraham" = COL_SA, "CGBS" = COL_CGBS)) +
    theme_classic(base_size = 20) +
    labs(title = paste(col_label, "|", row_label), x = NULL, y = NULL, color = NULL, fill = NULL) +
    theme(legend.position = "none",
          plot.title      = element_text(face = "bold", size = 7.5),
          axis.text       = element_text(size = 6))
}

cell_plots <- list()
for (r in c("Employment","Wage","Output"))
  for (cc in c("HS","MS","LS","Total"))
    cell_plots[[length(cell_plots) + 1]] <- make_cell(grid_vars[[r]][[cc]], r, cc)

legend_strip <- ggplot(data.frame(x=1, y=1, g=c("Sun & Abraham","CGBS")), aes(x, y, color=g)) +
  geom_line() +
  scale_color_manual(values = c("Sun & Abraham" = COL_SA, "CGBS" = COL_CGBS)) +
  theme_void() +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 9))

gridExtra::grid.arrange(
  grobs  = cell_plots, ncol = 4,
  top    = grid::textGrob("Comprehensive ATT Grid: SA vs CGBS across Skill Groups and Outcomes",
                          gp = grid::gpar(fontface = "bold", fontsize = 11)),
  bottom = cowplot::get_legend(legend_strip)
)

# =============================================================================
# PART VIII — SUMMARY TABLE
# Section A: Discrete AI implementation
# Section B: Continuous AI investment
# Columns: OLS | TWFE | SA ATT | CGBS ATT (cont. only) | Goodman-Bacon ATT
# Standard errors in parentheses below each estimate
# =============================================================================


# --- Section A: Continuous -----------------------------------------------------
cont_tbl_vars   <- c("ln_L_H_cont","ln_L_total_cont", "ln_Wage_H_cont", "ln_Wage_M_cont", "ln_Wage_L_cont",
                     "ln_Y_H_cont","ln_Y_total_cont", "Labor_Share_cont")
cont_tbl_labels <- c("High-Skill Employment", "Total Employment", "High-Skill Wage", "Mid-Skill Wage", "Low-Skill Wage",
                     "Output - High Skill","Output Total","Labor Share")

cont_table <- lapply(seq_along(cont_tbl_vars), function(i) {
  v   <- cont_tbl_vars[i]
  idx <- which(vars_cont == v)
  sa  <- sa_att_se(mod_sa_cont[[idx]])
  cg  <- cgbs_results[[v]]
  data.frame(
    Variable   = cont_tbl_labels[i],
    OLS        = fmt(coef(mod_ols_cont[[idx]])["invest_post"]  * mean_dose, se(mod_ols_cont[[idx]])["invest_post"]  * mean_dose),
    TWFE       = fmt(coef(mod_twfe_cont[[idx]])["invest_post"] * mean_dose, se(mod_twfe_cont[[idx]])["invest_post"] * mean_dose),
    SA_ATT     = fmt(sa["att"],  sa["se"]),
    CGBS_ATT   = fmt(cg$att,     cg$att_se)
  )
}) %>% bind_rows()


cat("\n", strrep("=", 72), "\n")
cat("  SECTION B — CONTINUOUS AI INVESTMENT\n")
cat("  All estimates at mean dose (~$1M) | SEs in parentheses\n")
cat(strrep("=", 72), "\n")
print(cont_table, row.names = FALSE)


