source("00_criteria.R")
suppressMessages({library(dplyr); library(tidyr)})

files <- list.files(CONF_RAW, pattern = "\\.rds$", full.names = TRUE)
stopifnot(length(files) > 0)
conf <- bind_rows(lapply(files, readRDS)) %>%
  filter(!is.na(p_nie_omit), !is.na(p_te_omit), !is.na(p_nie_adj))
cat(sprintf("[03] %d usable confounding replications\n", nrow(conf)))

conf <- conf %>% mutate(
  med_omit = mediation_criteria(p_nie_omit, p_te_omit, p_nde_omit, delta_star(24)),
  med_adj  = mediation_criteria(p_nie_adj,  p_te_adj,  p_nde_adj,  delta_star(24)))

cell <- conf %>% mutate(gold = logrank_full < 0.05) %>%
  group_by(scenario, u) %>% summarise(
    N = n(), N_gold_neg = sum(!gold),
    AA_pct_omitting_U  = round(100*mean(med_omit == "grant"), 1),
    AA_pct_adjusting_U = round(100*mean(med_adj  == "grant"), 1),
    TypeI_omitting_U   = round(100*mean(med_omit[!gold] == "grant"), 1),
    TypeI_adjusting_U  = round(100*mean(med_adj[!gold]  == "grant"), 1),
    mean_Pr_NIE_omit = round(mean(p_nie_omit), 3),
    mean_Pr_NIE_adj  = round(mean(p_nie_adj),  3),
    gap    = mean(p_nie_omit - p_nie_adj),
    gap_se = sd(p_nie_omit - p_nie_adj)/sqrt(n()), .groups = "drop")

u0 <- cell %>% filter(u == 0) %>% select(scenario, gap0 = gap, gap0_se = gap_se)
tabS6 <- cell %>% left_join(u0, by = "scenario") %>%
  mutate(U_hazard_ratio = round(exp(u), 2),
         excess_over_u0 = round(gap - gap0, 4),
         excess_se = sqrt(gap_se^2 + gap0_se^2),
         excess_CI95 = sprintf("(%.4f, %.4f)",
                               (gap - gap0) - 1.96*excess_se,
                               (gap - gap0) + 1.96*excess_se),
         confounding_bias_detected =
           ifelse(u == 0, "-- (control)",
                  ifelse((gap - gap0) - 1.96*excess_se > 0, "yes", "no"))) %>%
  transmute(Scenario = SCEN_LAB[scenario], u, U_hazard_ratio, N, N_gold_neg,
            AA_pct_omitting_U, AA_pct_adjusting_U,
            TypeI_omitting_U, TypeI_adjusting_U,
            mean_Pr_NIE_omit, mean_Pr_NIE_adj,
            excess_over_u0, excess_CI95, confounding_bias_detected) %>%
  arrange(factor(Scenario, levels = SCEN_LAB[SCEN_ORD]), u)
write.csv(tabS6, file.path(OUT_TAB, "supp_tableS6_confounding_sensitivity.csv"),
          row.names = FALSE)

sn <- tabS6 %>% filter(Scenario == "Survival-Null")
ex <- sn$excess_over_u0[sn$u == 0.5]
cat(sprintf("[03] Survival-Null: AA omitting U %.0f%% (u=0) -> %.0f%% (u=0.5); excess bias at u=0.5 = %.4f %s\n",
            sn$AA_pct_omitting_U[sn$u == 0], sn$AA_pct_omitting_U[sn$u == 0.5],
            ex, sn$excess_CI95[sn$u == 0.5]))
stopifnot(round(ex, 2) == 0.06)
cat(sprintf("[03] SNull at delta_star(24)=%.2f: AA omit by u: %s | TypeI (grant|gold-) omit by u: %s\n",
            delta_star(24), paste(sn$AA_pct_omitting_U, collapse = " / "),
            paste(sn$TypeI_omitting_U, collapse = " / ")))
cat("[03] verification PASSED; table S6 written\n")
