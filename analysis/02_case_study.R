source("00_criteria.R")
suppressMessages({library(dplyr); library(tidyr)})

COHORTS <- c("309 all", "309 KRAS wt")
CUT_ORD <- c("6mo", "9mo", "12mo", "18mo", "24mo", "final")
CUT_LAB <- c(`6mo` = "6 months", `9mo` = "9 months", `12mo` = "12 months",
             `18mo` = "18 months", `24mo` = "24 months", `final` = "Final analysis")

files <- list.files(CASE_RAW, pattern = "\\.rds$", full.names = TRUE)
stopifnot(length(files) > 0)
cells <- lapply(files, readRDS)
grid <- bind_rows(lapply(cells, function(x) {
  if (!is.data.frame(x$eff)) return(NULL)
  pri <- x$eff[x$eff$is_primary, , drop = FALSE]
  if (!nrow(pri)) return(NULL)
  cbind(x$meta, pri[1, c("tau", "RR_TE", "TE_lo", "TE_hi", "P_TE",
                         "RR_NDE", "NDE_lo", "NDE_hi", "P_NDE",
                         "RR_NIE", "NIE_lo", "NIE_hi", "P_NIE")])
})) %>% filter(cohort %in% COHORTS) %>%
  mutate(cut = factor(cut, levels = CUT_ORD)) %>% arrange(cohort, cut)
cat(sprintf("[02] %d cells loaded for %s\n", nrow(grid),
            paste(COHORTS, collapse = " & ")))

CASE_DELTA <- 0.835
grid <- grid %>% mutate(
  med  = mediation_criteria(P_NIE, P_TE, P_NDE, CASE_DELTA),
  trad = traditional_criteria(p_ORR_int),
  med_display  = c(grant = "Grant", continue = "Continue", stop = "Decline")[med],
  trad_display = ifelse(trad == "grant", "Grant", "Decline"),

  gold = p_OS_final < 0.05,
  med_correct  = vindicated_correct(med, gold),
  trad_correct = vindicated_correct(trad, gold))

fmt_ci <- function(m, lo, hi) sprintf("%.3f (%.3f–%.3f)", m, lo, hi)

t3 <- grid %>% filter(cut == "12mo") %>% transmute(
  Cohort = ifelse(cohort == "309 all",
                  "Panitumumab + FOLFOX, all randomized",
                  "Panitumumab + FOLFOX, KRAS wild-type"),
  n,
  `ORR, control -> treated` = sprintf("%.1f%% -> %.1f%% (p = %.3f)",
                                      100*ORR0, 100*ORR1, p_ORR_int),
  `OS events, interim / final (maturity)` =
    sprintf("%d / %d (%.1f%%)", os_ev_int, os_ev_final, 100*maturity),
  `RR_NIE (95% CrI); Pr(>1)` = sprintf("%s; %.2f", fmt_ci(RR_NIE, NIE_lo, NIE_hi), P_NIE),
  `RR_TE (95% CrI); Pr(>1)`  = sprintf("%s; %.3f", fmt_ci(RR_TE, TE_lo, TE_hi), P_TE),
  `Mediation decision` = med_display,
  `Traditional decision` = trad_display,
  `Final OS log-rank p` = sprintf("%.2f", p_OS_final))
write.csv(t3, file.path(OUT_TAB, "table3_case_study.csv"), row.names = FALSE)

tS5 <- grid %>% transmute(
  Cohort = ifelse(cohort == "309 all",
                  "Panitumumab + FOLFOX, all randomized",
                  "Panitumumab + FOLFOX, KRAS wild-type"),
  Interim = CUT_LAB[as.character(cut)],
  Information_fraction_pct = round(100*maturity, 1),
  Pr_NIE_gt1 = round(P_NIE, 3),
  Pr_TE_gt1 = round(P_TE, 3),
  Mediation_decision = sprintf("%s (%s)", med_display,
                               ifelse(med_correct, "correct", "incorrect")),
  Traditional_decision = sprintf("%s (%s)", trad_display,
                                 ifelse(trad_correct, "correct", "incorrect")))
write.csv(tS5, file.path(OUT_TAB, "supp_tableS5_case_timing_sensitivity.csv"),
          row.names = FALSE)

chk <- grid %>% filter(cut == "12mo")
a <- chk[chk$cohort == "309 all", ]
stopifnot(round(a$RR_TE, 3) == 0.922, round(a$P_TE, 3) == 0.023,
          round(a$RR_NIE, 3) == 1.041, a$med_display == "Decline",
          a$trad_display == "Grant")
k <- chk[chk$cohort == "309 KRAS wt", ]
stopifnot(k$med_display == "Decline", k$trad_display == "Decline")
cat("[02] verification PASSED; table 3 and table S5 written\n")
