source("00_criteria.R")
suppressMessages({library(dplyr); library(tidyr); library(ggplot2)})

INTERIM_TPS <- c(12, 18, 24, 30)
FINE_GRID <- seq(0.30, 0.995, by = 0.005)
NO_BENEFIT <- c("sc_null", "sc_snull")
MEDIATED   <- c("sc_modl", "sc_med", "sc_mixed")

files <- list.files(SIM_RAW, pattern = "\\.rds$", full.names = TRUE)
sim <- bind_rows(lapply(files, readRDS)) %>%
  filter(rep <= 500, tp %in% INTERIM_TPS,
         converged, !is.na(p_nie_mv), !is.na(p_te_mv)) %>%
  mutate(gold = logrank_full < 0.05)

matlab <- sim %>% group_by(tp) %>%
  summarise(m = round(100*mean(maturity)), .groups = "drop") %>%
  mutate(lab = sprintf("Month %d (%d%% mature)", tp, m))
MATLAB <- setNames(matlab$lab, matlab$tp)

RULES <- c("Proposed criteria (NIE with TE gate)", "NIE alone")

VLINES <- tibble::tibble(tp = as.integer(names(DELTA_STAR)),
                         dstar = unname(DELTA_STAR))

sweep1 <- function(De) {
  bind_rows(lapply(RULES, function(rn) {
    te_gate <- if (rn == "NIE alone") -Inf else TE_GATE
    sim %>% mutate(
        med = mediation_criteria(p_nie_mv, p_te_mv, p_nde_mv, De, te_gate),
        ok  = vindicated_correct(med, gold)) %>%
      group_by(tp, scenario) %>%
      summarise(N = n(), grant = mean(med == "grant"),
                grant_neg = mean(med[!gold] == "grant"),
                grant_pos = mean(med[gold] == "grant"),
                n_neg = sum(!gold), n_pos = sum(gold),
                ok = mean(ok), .groups = "drop") %>%
      mutate(Delta = De, Rule = rn)
  }))
}
grid <- bind_rows(lapply(FINE_GRID, sweep1))

band <- function(r, n) 1.96 * sqrt(pmax(r*(1 - r), 0)/n)
ROW1 <- "Overall correct decision"
ROW2 <- "Power and type I error"
M_OV <- "Overall correct decision rate — mediation criteria"
M_TR <- "Overall correct decision rate — traditional ORR criteria"
M_PW <- "Mean power (mediated benefit)"
M_T1 <- "Worst-case type I error"

f5m <- grid %>% filter(Rule == RULES[1]) %>% group_by(tp, Delta) %>% summarise(
    t1 = max(grant_neg[scenario %in% NO_BENEFIT]),
    n_t1 = n_neg[scenario %in% NO_BENEFIT][which.max(grant_neg[scenario %in% NO_BENEFIT])],
    pw = mean(grant_pos[scenario %in% MEDIATED]), n_pw = sum(n_pos[scenario %in% MEDIATED]),
    ov = mean(ok), n_ov = sum(N), .groups = "drop")
f5 <- bind_rows(
  f5m %>% transmute(tp, Delta, metric = M_OV, rate = ov, n = n_ov, Row = ROW1),
  f5m %>% transmute(tp, Delta, metric = M_PW, rate = pw, n = n_pw, Row = ROW2),
  f5m %>% transmute(tp, Delta, metric = M_T1, rate = t1, n = n_t1, Row = ROW2)) %>%
  mutate(lo = pmax(0, rate - band(rate, n)), hi = pmin(1, rate + band(rate, n)))

f5_trad <- sim %>% mutate(trad = traditional_criteria(chisq_p),
                          ok = vindicated_correct(trad, gold)) %>%
  group_by(tp) %>% summarise(rate = mean(ok), .groups = "drop") %>%
  crossing(Delta = range(FINE_GRID)) %>%
  mutate(metric = M_TR, Row = ROW1, lo = NA_real_, hi = NA_real_)

f5 <- bind_rows(f5, f5_trad) %>%
  mutate(maturity = factor(MATLAB[as.character(tp)], levels = MATLAB),
         Row = factor(Row, levels = c(ROW1, ROW2)),
         metric = factor(metric, levels = c(M_OV, M_TR, M_PW, M_T1)))

ref_t1 <- expand.grid(y = c(5, 10), maturity = factor(MATLAB, levels = MATLAB),
                      Row = factor(ROW2, levels = c(ROW1, ROW2)))

p5 <- ggplot(f5, aes(Delta, 100*rate, colour = metric, fill = metric,
                     linetype = metric)) +
  geom_hline(data = ref_t1, aes(yintercept = y), inherit.aes = FALSE,
             linetype = "dotted", colour = "grey55") +
  geom_ribbon(aes(ymin = 100*lo, ymax = 100*hi), alpha = .18, colour = NA,
              na.rm = TRUE) +
  geom_line(linewidth = .6, na.rm = TRUE) +
  geom_vline(data = VLINES %>%
               mutate(maturity = factor(MATLAB[as.character(tp)], levels = MATLAB)),
             aes(xintercept = dstar), inherit.aes = FALSE,
             linetype = "dashed", colour = "grey35") +
  facet_grid(Row ~ maturity) +
  scale_colour_manual(values = setNames(c("#6a7f2a", "#c9962e", "#2f7d8f", "#b3402f"),
                                        c(M_OV, M_TR, M_PW, M_T1)),
                      aesthetics = c("colour", "fill")) +
  scale_linetype_manual(values = setNames(c("solid", "22", "solid", "solid"),
                                          c(M_OV, M_TR, M_PW, M_T1))) +
  scale_x_continuous(breaks = seq(0.3, 1.0, 0.1)) +
  labs(x = expression("Posterior probability threshold " * Delta),
       y = "Rate (%)", colour = NULL, fill = NULL, linetype = NULL,
       caption = paste("Bands: 95% Monte-Carlo interval. Dotted lines: 5% and 10%",
                       "type I references. Dashed vertical: the calibrated \u0394* of each",
                       "interim (0.825, 0.835, 0.83, 0.775). Type I and power condition",
                       "on each trial's own final OS log-rank result.")) +
  theme_bw(base_size = 10) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey93"),
        plot.caption = element_text(size = 7.5, colour = "grey30")) +
  guides(colour = guide_legend(nrow = 2), fill = guide_legend(nrow = 2),
         linetype = guide_legend(nrow = 2))
ggsave(file.path(OUT_FIG, "fig5_delta_calibration_by_maturity.pdf"), p5,
       width = 10.5, height = 5.6, device = cairo_pdf)
ggsave(file.path(OUT_FIG, "fig5_delta_calibration_by_maturity.png"), p5,
       width = 10.5, height = 5.6, dpi = 300)

PANELS <- c("Mediation criteria (proposed)", "Traditional ORR criteria")
err_scen <- c("sc_null","sc_snull","sc_modl","sc_med","sc_mixed")

f6_med <- grid %>% filter(Rule == RULES[1], scenario %in% err_scen) %>%
  transmute(tp, scenario, Delta,
            err = ifelse(scenario %in% NO_BENEFIT, grant_neg, 1 - grant_pos),
            Panel = PANELS[1])

f6_trad <- sim %>% mutate(trad = traditional_criteria(chisq_p)) %>%
  group_by(tp, scenario) %>%
  summarise(grant_neg = mean(trad[!gold] == "grant"),
            grant_pos = mean(trad[gold] == "grant"), .groups = "drop") %>%
  filter(scenario %in% err_scen) %>%
  mutate(err = ifelse(scenario %in% NO_BENEFIT, grant_neg, 1 - grant_pos),
         Panel = PANELS[2]) %>%
  crossing(Delta = range(FINE_GRID)) %>%
  select(tp, scenario, Delta, err, Panel)

f6 <- bind_rows(f6_med, f6_trad) %>%
  mutate(type = ifelse(scenario %in% NO_BENEFIT, "False approval", "Missed approval"),
         Scenario = factor(SCEN_LAB[scenario], levels = SCEN_LAB[err_scen]),
         maturity = factor(MATLAB[as.character(tp)], levels = MATLAB),
         Panel = factor(Panel, levels = PANELS))
p6 <- ggplot(f6, aes(Delta, 100*err, colour = Scenario, linetype = type)) +
  geom_line(linewidth = .55) +
  geom_vline(data = VLINES %>%
               mutate(maturity = factor(MATLAB[as.character(tp)], levels = MATLAB)),
             aes(xintercept = dstar), inherit.aes = FALSE,
             linetype = "dotted", colour = "grey25") +
  facet_grid(Panel ~ maturity) +
  scale_colour_manual(values = c("Null" = "#e07a5f", "Survival-Null" = "#a68a00",
                                 "Moderate-Local" = "#2a9d8f",
                                 "Mediated-Only" = "#33a1fd", "Mixed" = "#e26ee5")) +
  scale_linetype_manual(values = c("False approval" = "solid",
                                   "Missed approval" = "22")) +
  scale_x_continuous(breaks = seq(0.3, 0.9, 0.2)) +
  labs(x = expression("Posterior probability threshold " * Delta),
       y = "Error rate (%)", colour = NULL, linetype = NULL,
       caption = paste("Solid: false approval. Dashed: missed approval. Bottom row:",
                       "the traditional criterion has no threshold, so its error",
                       "rates are constant in \u0394. Dotted vertical: the calibrated",
                       "\u0394* of each interim.")) +
  theme_bw(base_size = 10) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey93"),
        plot.caption = element_text(size = 7.5, colour = "grey30"))
ggsave(file.path(OUT_FIG, "fig6_error_tradeoff_by_scenario.pdf"), p6,
       width = 10.5, height = 5.6, device = cairo_pdf)
ggsave(file.path(OUT_FIG, "fig6_error_tradeoff_by_scenario.png"), p6,
       width = 10.5, height = 5.6, dpi = 300)

PANEL_LAB <- c(sc_null  = "Null (false approval)",
               sc_snull = "Survival-Null (false approval)",
               sc_modl  = "Moderate-Local (missed approval)",
               sc_med   = "Mediated-Only (missed approval)",
               sc_mixed = "Mixed (missed approval)")
MCOL <- setNames(c("#9ecad6", "#5ba3b6", "#2f7d8f", "#153f4a"), c(12, 18, 24, 30))

f6b_med <- grid %>% filter(Rule == RULES[1], scenario %in% err_scen) %>%
  transmute(tp, scenario, Delta,
            err = ifelse(scenario %in% NO_BENEFIT, grant_neg, 1 - grant_pos),
            Criteria = "Mediation (proposed)")
f6b_trad <- sim %>% mutate(trad = traditional_criteria(chisq_p)) %>%
  group_by(tp, scenario) %>%
  summarise(grant_neg = mean(trad[!gold] == "grant"),
            grant_pos = mean(trad[gold] == "grant"), .groups = "drop") %>%
  filter(scenario %in% err_scen) %>%
  mutate(err = ifelse(scenario %in% NO_BENEFIT, grant_neg, 1 - grant_pos),
         Criteria = "Traditional (ORR)") %>%
  crossing(Delta = range(FINE_GRID)) %>%
  select(tp, scenario, Delta, err, Criteria)

f6b <- bind_rows(f6b_med, f6b_trad) %>%
  mutate(Panel = factor(PANEL_LAB[scenario], levels = PANEL_LAB),
         Maturity = factor(MATLAB[as.character(tp)], levels = MATLAB),
         Criteria = factor(Criteria,
                           levels = c("Mediation (proposed)", "Traditional (ORR)")))
p6b <- ggplot(f6b, aes(Delta, 100*err, colour = Maturity, linetype = Criteria)) +
  geom_line(linewidth = .6) +
  geom_vline(data = VLINES %>%
               mutate(Maturity = factor(MATLAB[as.character(tp)], levels = MATLAB)),
             aes(xintercept = dstar, colour = Maturity), inherit.aes = FALSE,
             linetype = "dotted", linewidth = .45, show.legend = FALSE) +
  facet_wrap(~ Panel, ncol = 3) +
  scale_colour_manual(values = setNames(MCOL, MATLAB)) +
  scale_linetype_manual(values = c("Mediation (proposed)" = "solid",
                                   "Traditional (ORR)" = "22")) +
  scale_x_continuous(breaks = seq(0.3, 0.9, 0.2)) +
  labs(x = expression("Posterior probability threshold " * Delta),
       y = "Error rate (%)", colour = NULL, linetype = NULL,
       caption = paste("Solid: mediation criteria as a function of \u0394. Dashed:",
                       "traditional ORR criteria, constant in \u0394. Dotted",
                       "verticals: the calibrated \u0394* of each interim.")) +
  theme_bw(base_size = 10) +
  theme(legend.position = "top", legend.box = "horizontal",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey93"),
        plot.caption = element_text(size = 7.5, colour = "grey30")) +
  guides(colour = guide_legend(order = 1, nrow = 2),
         linetype = guide_legend(order = 2, nrow = 2,
                                 override.aes = list(colour = "grey20")))
ggsave(file.path(OUT_FIG, "fig6alt_error_tradeoff_by_scenario.pdf"), p6b,
       width = 9.6, height = 6.4, device = cairo_pdf)
ggsave(file.path(OUT_FIG, "fig6alt_error_tradeoff_by_scenario.png"), p6b,
       width = 9.6, height = 6.4, dpi = 300)

cur <- grid %>% filter(Rule == RULES[1]) %>% group_by(tp, Delta) %>%
  summarise(t1 = 100*max(grant_neg[scenario %in% NO_BENEFIT]),
            pw = 100*mean(grant_pos[scenario %in% MEDIATED]),
            acc = 100*mean(ok), .groups = "drop")
tabS7 <- cur %>% group_by(tp) %>% summarise(
  Maturity_pct = round(100*mean(sim$maturity[sim$tp == tp[1]])),
  Delta_star_5pct = min(Delta[t1 <= 100*ALPHA_BUDGET]),
  TypeI_at_star = t1[Delta == Delta_star_5pct][1],
  Power_at_star = pw[Delta == Delta_star_5pct][1],
  Accuracy_at_star = acc[Delta == Delta_star_5pct][1],
  Delta_star_10pct = min(Delta[t1 <= 10]),
  Accuracy_plateau = sprintf("%.2f-%.2f", min(Delta[acc >= max(acc) - 1]),
                             max(Delta[acc >= max(acc) - 1])), .groups = "drop") %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))
write.csv(tabS7, file.path(OUT_TAB, "supp_tableS7_delta_calibration.csv"),
          row.names = FALSE)
cat("[05] Table S7 (error-budget calibration):\n")
print(as.data.frame(tabS7), row.names = FALSE)
stopifnot(all(tabS7$Delta_star_5pct == unname(DELTA_STAR[as.character(tabS7$tp)])))

tabD <- grid %>%
  inner_join(VLINES, by = "tp") %>% filter(abs(Delta - dstar) < 1e-9) %>%
  group_by(Rule, tp) %>%
  summarise(worst_typeI = round(100*max(grant_neg[scenario %in% NO_BENEFIT]), 1),
            mean_power  = round(100*mean(grant_pos[scenario %in% MEDIATED]), 1),
            overall_correct = round(100*mean(ok), 1),
            missed_ModerateLocal = round(100*(1 - grant_pos[scenario == "sc_modl"]), 1),
            .groups = "drop")
tabD_trad <- sim %>% mutate(trad = traditional_criteria(chisq_p),
                            ok = vindicated_correct(trad, gold)) %>%
  group_by(tp, scenario) %>%
  summarise(N = n(), grant_neg = mean(trad[!gold] == "grant"),
            grant_pos = mean(trad[gold] == "grant"), ok = mean(ok), .groups = "drop") %>%
  group_by(tp) %>%
  summarise(Rule = "Traditional ORR (Delta-independent)",
            worst_typeI = round(100*max(grant_neg[scenario %in% NO_BENEFIT]), 1),
            mean_power  = round(100*mean(grant_pos[scenario %in% MEDIATED]), 1),
            overall_correct = round(100*mean(ok), 1),
            missed_ModerateLocal = round(100*(1 - grant_pos[scenario == "sc_modl"]), 1),
            .groups = "drop") %>% select(Rule, everything())
tabD <- bind_rows(tabD, tabD_trad)
write.csv(tabD, file.path(OUT_TAB, "tableD_calibration_at_delta090.csv"),
          row.names = FALSE)
cat("[05] At the calibrated delta_star of each interim:\n")
print(as.data.frame(tabD), row.names = FALSE)

cur24 <- tabD %>% filter(Rule == RULES[1], tp == 24)
stopifnot(nrow(cur24) == 1, cur24$worst_typeI <= 5)
trad24 <- tabD %>% filter(grepl("Traditional", Rule), tp == 24)
stopifnot(nrow(trad24) == 1, trad24$worst_typeI == 100,
          trad24$overall_correct == 60.6)
cat("[05] verification PASSED; figures 5-6 and table D written\n")
