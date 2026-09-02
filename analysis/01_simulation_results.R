source("00_criteria.R")
suppressMessages({library(dplyr); library(tidyr); library(ggplot2)})

INTERIM_TPS <- c(12, 18, 24, 30)
REF_TP <- 24

files <- list.files(SIM_RAW, pattern = "\\.rds$", full.names = TRUE)
stopifnot(length(files) > 0)
raw <- bind_rows(lapply(files, readRDS)) %>%
  filter(rep <= 500, tp %in% INTERIM_TPS)

sim <- raw %>% filter(converged, !is.na(p_nie_mv), !is.na(p_te_mv))
cat(sprintf("[01] %d of %d replication-cuts usable\n", nrow(sim), nrow(raw)))

decide <- function(d) {
  d %>% mutate(
    gold = logrank_full < 0.05,
    med  = mediation_criteria(p_nie_mv, p_te_mv, p_nde_mv, delta_star(tp)),
    trad = traditional_criteria(chisq_p),
    ok_med  = vindicated_correct(med,  gold),
    ok_trad = vindicated_correct(trad, gold))
}

d24 <- decide(sim %>% filter(tp == REF_TP))
tabS4 <- d24 %>% group_by(scenario) %>% summarise(
    N = n(),
    P_AA = 100*mean(med == "grant"),
    P_Continue = 100*mean(med == "continue"),
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
  mutate(Scenario = SCEN_LAB[scenario], Truth = SCEN_TRUTH[scenario],
         .before = 1) %>% select(-scenario) %>%
  mutate(across(where(is.numeric), ~round(.x, 1)))
tabS4 <- bind_rows(tabS4,
  tabS4 %>% summarise(Scenario = "OVERALL (mean of 6)", Truth = "", N = sum(N),
                      across(c(P_AA:Gold_positive, Correct_mediation,
                               Correct_traditional), ~round(mean(.x), 1))))
write.csv(tabS4, file.path(OUT_TAB, "supp_tableS4_decision_distribution.csv"),
          row.names = FALSE)

tabC <- decide(sim) %>% group_by(tp, scenario) %>%
  summarise(N = n(), maturity_pct = round(100*mean(maturity), 0),
            Correct_mediation = round(100*mean(ok_med), 1),
            Correct_traditional = round(100*mean(ok_trad), 1), .groups = "drop") %>%
  mutate(Scenario = SCEN_LAB[scenario], .after = tp) %>% select(-scenario)
tabC <- bind_rows(tabC,
  tabC %>% group_by(tp) %>%
    summarise(Scenario = "OVERALL (mean of 6)", N = sum(N),
              maturity_pct = round(mean(maturity_pct), 0),
              Correct_mediation = round(mean(Correct_mediation), 1),
              Correct_traditional = round(mean(Correct_traditional), 1),
              .groups = "drop")) %>% arrange(tp)
write.csv(tabC, file.path(OUT_TAB, "tableC_maturity_sensitivity.csv"),
          row.names = FALSE)

tab4 <- decide(sim) %>% group_by(tp, scenario) %>% summarise(
    N = n(),
    Gold_pos_pct = round(100*mean(gold), 1),
    N_gold_neg = sum(!gold), N_gold_pos = sum(gold),
    Acc_med  = round(100*mean(ok_med), 1),  Acc_trad = round(100*mean(ok_trad), 1),
    TypeI_med  = round(100*mean(med[!gold] == "grant"), 1),
    TypeI_trad = round(100*mean(trad[!gold] == "grant"), 1),
    Power_med  = round(100*mean(med[gold] == "grant"), 1),
    Power_trad = round(100*mean(trad[gold] == "grant"), 1), .groups = "drop") %>%
  mutate(Scenario = SCEN_LAB[scenario],
         Maturity = paste0("Month ", tp),
         Delta_star = delta_star(tp), .before = 1) %>%
  select(Maturity, Delta_star, Scenario, N, Gold_pos_pct, N_gold_neg, N_gold_pos,
         Acc_med, Acc_trad, TypeI_med, TypeI_trad, Power_med, Power_trad) %>%
  arrange(match(Maturity, paste0("Month ", c(12,18,24,30))),
          match(Scenario, SCEN_LAB[SCEN_ORD]))
write.csv(tab4, file.path(OUT_TAB, "table4_oc_by_maturity_scenario.csv"),
          row.names = FALSE)
cat(sprintf("[01] table 4 written (%d rows)\n", nrow(tab4)))

f4 <- d24 %>% group_by(scenario) %>%
  summarise(`Mediation (proposed)` = 100*mean(ok_med),
            `Traditional (ORR)` = 100*mean(ok_trad), .groups = "drop") %>%
  pivot_longer(-scenario, names_to = "Criteria", values_to = "rate") %>%
  mutate(Scenario = factor(SCEN_LAB[scenario], levels = SCEN_LAB[SCEN_ORD]),
         Criteria = factor(Criteria, levels = names(COL2)))
p4 <- ggplot(f4, aes(Scenario, rate, fill = Criteria)) +
  geom_col(position = position_dodge(width = .78), width = .7) +
  geom_text(aes(label = sprintf("%.1f", rate)),
            position = position_dodge(width = .78), vjust = -0.45, size = 3.0) +
  scale_fill_manual(values = COL2) +
  scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(x = NULL, y = "Correct decision rate (%)", fill = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(OUT_FIG, "fig4_correct_decision_rate_v4.pdf"), p4,
       width = 7.2, height = 4.2, device = cairo_pdf)
ggsave(file.path(OUT_FIG, "fig4_correct_decision_rate_v4.png"), p4,
       width = 7.2, height = 4.2, dpi = 300)

f6 <- tabC %>% filter(Scenario == "OVERALL (mean of 6)") %>%
  transmute(month = tp, maturity_pct,
            `Mediation (proposed)` = Correct_mediation,
            `Traditional (ORR)` = Correct_traditional) %>%
  pivot_longer(c(`Mediation (proposed)`, `Traditional (ORR)`),
               names_to = "Criteria", values_to = "rate") %>%
  mutate(Criteria = factor(Criteria, levels = names(COL2)))
p6 <- ggplot(f6, aes(month, rate, colour = Criteria)) +
  geom_line(linewidth = .9) + geom_point(size = 2.4) +
  geom_text(data = distinct(f6, month, maturity_pct),
            aes(month, 101.5, label = sprintf("%d%%", maturity_pct)),
            inherit.aes = FALSE, size = 2.9, colour = "grey35") +
  annotate("text", x = 12.2, y = 106, label = "OS data maturity",
           hjust = 0, size = 2.9, colour = "grey35") +
  scale_colour_manual(values = COL2) +
  scale_x_continuous(breaks = INTERIM_TPS) +
  scale_y_continuous(limits = c(45, 107), breaks = seq(50, 100, 10)) +
  labs(x = "Interim analysis time (months from trial start)",
       y = "Overall correct decision rate (%)", colour = NULL) +
  theme_classic(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(OUT_FIG, "fig7_maturity_sensitivity.pdf"), p6,
       width = 6.4, height = 4.2, device = cairo_pdf)
ggsave(file.path(OUT_FIG, "fig7_maturity_sensitivity.png"), p6,
       width = 6.4, height = 4.2, dpi = 300)

ov <- tabS4[tabS4$Scenario == "OVERALL (mean of 6)", ]
cat(sprintf("[01] month-24 overall at delta_star=%.2f: mediation %.1f | traditional %.1f\n", delta_star(24),
            ov$Correct_mediation, ov$Correct_traditional))
cat("[01] month-24 per-scenario at delta_star:\n")
print(as.data.frame(tabS4), row.names = FALSE)
cat("[01] maturity overall rows:\n")
print(as.data.frame(tabC[tabC$Scenario == "OVERALL (mean of 6)", ]), row.names = FALSE)
stopifnot(ov$Correct_traditional == 60.6)
cat("[01] figures 4 and 7 and tables S4/C written\n")
