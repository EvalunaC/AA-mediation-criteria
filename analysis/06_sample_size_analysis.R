source("00_criteria.R")
suppressMessages({library(dplyr); library(tidyr); library(ggplot2)})

SS_RAW <- file.path(BASE, "data", "sample_size")
NO_B <- c("sc_null","sc_snull"); MED <- c("sc_modl","sc_med","sc_mixed")
DSTAR24 <- delta_star(24)

fs <- list.files(SS_RAW, pattern = "\\.rds$", full.names = TRUE)
if (length(fs) == 0)
  stop("no sample-size fits found — run generate_sample_size.sh first")
if (length(fs) < 9000)
  cat(sprintf("[06] WARNING: %d of 9000 sample-size files present\n", length(fs)))
ss <- bind_rows(lapply(fs, readRDS))

f500 <- list.files(SIM_RAW, pattern = "\\.rds$", full.names = TRUE)
d500 <- bind_rows(lapply(f500, readRDS)) %>%
  filter(rep <= 500, tp == 24) %>% mutate(n_target = 500L) %>%
  select(any_of(names(ss)))

d <- bind_rows(ss, d500)
fit_rates <- d %>% group_by(n_target, scenario) %>%
  summarise(N_attempted = n(), fit_rate = mean(converged), .groups = "drop")
d <- d %>% filter(converged, !is.na(p_nie_mv), !is.na(p_te_mv)) %>%
  mutate(gold = logrank_full < 0.05,
         med  = mediation_criteria(p_nie_mv, p_te_mv, p_nde_mv, DSTAR24),
         trad = traditional_criteria(chisq_p),
         ok_med  = vindicated_correct(med,  gold),
         ok_trad = vindicated_correct(trad, gold))

cell <- d %>% group_by(n_target, scenario) %>% summarise(
    N = n(),
    P_AA = 100*mean(med == "grant"), P_Continue = 100*mean(med == "continue"),
    P_Stop = 100*mean(med == "stop"),
    Traditional_AA = 100*mean(trad == "grant"),
    Gold_positive = 100*mean(gold),
    N_gold_neg = sum(!gold), N_gold_pos = sum(gold),
    TypeI_med  = 100*mean(med[!gold] == "grant"),
    TypeI_trad = 100*mean(trad[!gold] == "grant"),
    Power_med  = 100*mean(med[gold] == "grant"),
    Power_trad = 100*mean(trad[gold] == "grant"),
    Correct_mediation = 100*mean(ok_med),
    Correct_traditional = 100*mean(ok_trad), .groups = "drop") %>%
  left_join(fit_rates, by = c("n_target", "scenario")) %>%
  mutate(Scenario = SCEN_LAB[scenario], Truth = SCEN_TRUTH[scenario],
         Fit_rate_pct = round(100*fit_rate, 1), .after = n_target) %>%
  select(-scenario, -fit_rate, -N_attempted) %>%
  mutate(across(where(is.numeric), ~round(.x, 1)))
ov <- cell %>% group_by(n_target) %>%
  summarise(Scenario = "OVERALL (mean of 6)", Truth = "", Fit_rate_pct = NA,
            N = sum(N), across(c(P_AA:Gold_positive, Correct_mediation,
                                 Correct_traditional), ~round(mean(.x), 1)),
            .groups = "drop")
tabS8 <- bind_rows(cell, ov) %>% arrange(n_target)

recal <- bind_rows(lapply(sort(unique(d$n_target)), function(nn) {
  x <- d %>% filter(n_target == nn)
  sw <- bind_rows(lapply(seq(0.30, 0.995, 0.005), function(De) {
    x %>% mutate(m = mediation_criteria(p_nie_mv, p_te_mv, p_nde_mv, De)) %>%
      group_by(scenario) %>%
      summarise(g_neg = 100*mean(m[!gold] == "grant"),
                g_pos = 100*mean(m[gold] == "grant"), .groups = "drop") %>%
      summarise(Delta = De, t1 = max(g_neg[scenario %in% NO_B]),
                pw = mean(g_pos[scenario %in% MED]))
  }))
  ok <- sw %>% filter(t1 <= 5)
  tibble(n_target = nn,
         Delta_star_recalibrated = if (nrow(ok)) min(ok$Delta) else NA,
         TypeI_at_recal = if (nrow(ok)) ok$t1[which.min(ok$Delta)] else NA,
         Power_at_recal = if (nrow(ok)) round(ok$pw[which.min(ok$Delta)], 1) else NA)
}))
write.csv(tabS8, file.path(OUT_TAB, "supp_tableS8_sample_size_sensitivity.csv"),
          row.names = FALSE)
write.csv(recal, file.path(OUT_TAB, "supp_tableS8b_recalibrated_delta_by_n.csv"),
          row.names = FALSE)

cat("[06] Table S8 overall rows:\n")
print(as.data.frame(tabS8 %>% filter(Scenario == "OVERALL (mean of 6)")), row.names = FALSE)
cat("[06] Recalibrated Delta* by n (5% budget):\n")
print(as.data.frame(recal), row.names = FALSE)

long <- bind_rows(
  cell %>% transmute(n_target, scenario = Scenario, Metric = "Correct decision rate",
                     Mediation = Correct_mediation, Traditional = Correct_traditional),
  cell %>%
    transmute(n_target, scenario = Scenario,
              Metric = "Type I error: grant | final OS log-rank negative",
              Mediation = TypeI_med, Traditional = TypeI_trad),
  cell %>%
    transmute(n_target, scenario = Scenario,
              Metric = "Power: grant | final OS log-rank positive",
              Mediation = Power_med, Traditional = Power_trad)) %>%
  pivot_longer(c(Mediation, Traditional), names_to = "Criteria", values_to = "rate") %>%
  mutate(Scenario = factor(scenario, levels = SCEN_LAB[SCEN_ORD]),
         Metric = factor(Metric, levels = c("Correct decision rate",
                 "Type I error: grant | final OS log-rank negative",
                 "Power: grant | final OS log-rank positive")))

NO_BENEFIT_LAB <- unname(SCEN_LAB[names(SCEN_TRUTH)[SCEN_TRUTH == "No benefit"]])
BENEFIT_LAB    <- unname(SCEN_LAB[names(SCEN_TRUTH)[SCEN_TRUTH != "No benefit"]])
M_T1  <- "Type I error: grant | final OS log-rank negative"
M_PWR <- "Power: grant | final OS log-rank positive"
long <- long %>%
  filter(as.character(Metric) == "Correct decision rate" |
         (as.character(Metric) == M_T1  & scenario %in% NO_BENEFIT_LAB) |
         (as.character(Metric) == M_PWR & scenario %in% BENEFIT_LAB))
cat(sprintf("[06] figure S1: type I panel restricted to %s; power panel to %s\n",
            paste(NO_BENEFIT_LAB, collapse = "/"), paste(BENEFIT_LAB, collapse = "/")))
SCOL <- c("Null"="#e07a5f","Survival-Null"="#a68a00","Moderate-Local"="#2a9d8f",
          "Mediated-Only"="#33a1fd","Direct-Only"="#6b705c","Mixed"="#e26ee5")
pS1 <- ggplot(long, aes(n_target, rate, colour = Scenario,
                        linetype = Criteria, shape = Criteria)) +
  geom_hline(data = data.frame(
             Metric = factor("Type I error: grant | final OS log-rank negative",
             levels = levels(long$Metric)), y = 5), aes(yintercept = y),
             linetype = "dotted", colour = "grey40") +
  geom_line(linewidth = .6) + geom_point(size = 2, fill = "white") +
  facet_wrap(~Metric, ncol = 3) +
  scale_colour_manual(values = SCOL, drop = FALSE) +
  scale_linetype_manual(values = c(Mediation = "solid", Traditional = "22")) +
  scale_shape_manual(values = c(Mediation = 16, Traditional = 21)) +
  scale_x_continuous(breaks = c(100, 200, 300, 500)) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Total sample size", y = "Rate (%)",
       colour = NULL, linetype = NULL, shape = NULL,
       caption = paste0("Interim at 24 months; decisions at the design-calibrated ",
                sprintf("Δ* = %.2f. ", DSTAR24), "A trial counts as having a benefit ",
                "only if its own final OS log-rank is significant.\n",
                "Type I error is shown for the two no-benefit scenarios and power for ",
                "the four benefit scenarios, the only settings in which each is ",
                "interpretable; the correct decision rate is shown for all six. ",
                "Dotted line: the 5% false-approval budget.")) +
  theme_bw(base_size = 10) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey93"),
        plot.caption = element_text(size = 7.5, colour = "grey30"))
ggsave(file.path(OUT_FIG, "supp_figS1_sample_size.pdf"), pS1,
       width = 10.2, height = 4.6, device = cairo_pdf)
ggsave(file.path(OUT_FIG, "supp_figS1_sample_size.png"), pS1,
       width = 10.2, height = 4.6, dpi = 300)

v500 <- tabS8 %>% filter(n_target == 500, Scenario == "OVERALL (mean of 6)")
stopifnot(v500$Correct_traditional == 60.6, v500$Correct_mediation > 90)
cat("[06] verification PASSED; table S8 and figure S1 written\n")
