source("00_criteria.R")
suppressMessages({library(dplyr); library(tidyr); library(ggplot2)})

TP <- 24
SCENARIOS <- c("sc_null", "sc_snull", "sc_modl", "sc_med", "sc_direct", "sc_mixed")
SCEN_LAB  <- c(sc_null = "Null", sc_snull = "Survival-Null",
               sc_modl = "Moderate-Local", sc_med = "Mediated-Only",
               sc_direct = "Direct-Only", sc_mixed = "Mixed")
TEAL <- "#2f7d8f"

DEC_MED   <- c("Grant AA", "Continue", "Stop")
DEC_TRAD  <- c("Approved", "Not approved")
DEC_LEVELS  <- c(DEC_MED, DEC_TRAD)
CRIT_LEVELS <- c("Traditional (ORR)", "Mediation (proposed)")
GOLD_LEVELS <- c("OS success", "OS insignificant")

files <- list.files(SIM_RAW, pattern = "\\.rds$", full.names = TRUE)
d <- bind_rows(lapply(files, readRDS)) %>%
  filter(tp == TP, converged, !is.na(p_nie_mv), !is.na(p_te_mv)) %>%
  mutate(gold = logrank_full < 0.05,
         med  = mediation_criteria(p_nie_mv, p_te_mv, p_nde_mv, delta_star(TP)),
         trad = traditional_criteria(chisq_p))

long <- bind_rows(
  d %>% transmute(scenario, gold, Criteria = "Mediation (proposed)",
                  Decision = recode(as.character(med), grant = "Grant AA",
                                    continue = "Continue", stop = "Stop")),
  d %>% transmute(scenario, gold, Criteria = "Traditional (ORR)",
                  Decision = recode(as.character(trad), grant = "Approved",
                                    decline = "Not approved"))) %>%
  mutate(Scenario = factor(SCEN_LAB[scenario], levels = SCEN_LAB[SCENARIOS]),
         Criteria = factor(Criteria, levels = CRIT_LEVELS),
         Decision = factor(Decision, levels = DEC_LEVELS),
         Gold     = factor(ifelse(gold, GOLD_LEVELS[1], GOLD_LEVELS[2]),
                           levels = GOLD_LEVELS))

grid <- bind_rows(
  expand.grid(Criteria = "Mediation (proposed)", Decision = DEC_MED,
              Gold = GOLD_LEVELS, stringsAsFactors = FALSE),
  expand.grid(Criteria = "Traditional (ORR)", Decision = DEC_TRAD,
              Gold = GOLD_LEVELS, stringsAsFactors = FALSE))
grid <- merge(data.frame(Scenario = SCEN_LAB[SCENARIOS], stringsAsFactors = FALSE), grid)

cm <- long %>% count(Scenario, Criteria, Decision, Gold, name = "n") %>%
  right_join(grid, by = c("Scenario", "Criteria", "Decision", "Gold")) %>%
  mutate(n = ifelse(is.na(n), 0L, n),
         Scenario = factor(Scenario, levels = SCEN_LAB[SCENARIOS]),
         Criteria = factor(Criteria, levels = CRIT_LEVELS),
         Decision = factor(Decision, levels = DEC_LEVELS),
         Gold     = factor(Gold, levels = GOLD_LEVELS)) %>%
  group_by(Scenario, Criteria) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>%
  mutate(correct = (Decision %in% c("Grant AA", "Continue", "Approved") &
                      Gold == GOLD_LEVELS[1]) |
                   (Decision %in% c("Stop", "Not approved") & Gold == GOLD_LEVELS[2]))

out <- cm %>%
  group_by(Scenario, Criteria) %>% mutate(N_panel = sum(n)) %>% ungroup() %>%
  mutate(scenario_order = match(Scenario, SCEN_LAB[SCENARIOS]),
         criteria_order = match(Criteria, CRIT_LEVELS),
         decision_order = match(Decision, DEC_LEVELS),
         gold_order     = match(Gold, GOLD_LEVELS),
         pct = round(pct, 1), interim_month = TP, Delta = delta_star(TP)) %>%
  arrange(criteria_order, scenario_order, decision_order, gold_order) %>%
  select(Scenario, scenario_order, Criteria, criteria_order, Decision, decision_order,
         Gold, gold_order, n, pct, N_panel, correct, interim_month, Delta)

write.csv(out, file.path(OUT_TAB, "supp_tableS9_confusion_matrix.csv"), row.names = FALSE)
cat(sprintf("[08] wrote supp_tableS9 (%d rows)\n", nrow(out)))

NMAX <- max(cm$n)
p <- ggplot(cm, aes(Decision, Gold)) +
  geom_tile(aes(fill = n), colour = "white", linewidth = 0.7) +

  geom_tile(data = subset(cm, correct), fill = NA, colour = "#16181c",
            linewidth = 0.7) +
  geom_text(aes(label = n, colour = n > 0.55 * NMAX), size = 2.9) +
  facet_grid(Scenario ~ Criteria, scales = "free_x", space = "free_x", switch = "y") +
  scale_fill_gradient(low = "#ffffff", high = TEAL, limits = c(0, NMAX),
                      name = "Trials (n)") +
  scale_colour_manual(values = c(`FALSE` = "#16181c", `TRUE` = "#ffffff"),
                      guide = "none") +
  scale_y_discrete(limits = rev(GOLD_LEVELS)) +
  labs(x = "Predicted outcome", y = "Actual outcome",
       title = sprintf("Confusion matrix of interim decisions against the final OS log-rank (month %d, %s* = %.3g)",
                       TP, "Δ", delta_star(TP)),
       subtitle = paste0(
         "Actual outcome is the trial's own final OS log-rank; OS success = significant at complete follow-up.\n",
         "Cells are counts of trials; denominators in the row headings. Outlined cells are vindicated.\n",
         "The traditional criterion is binary (two columns); the mediation criteria add a continuation arm.")) +
  theme_classic(base_size = 11) +
  theme(legend.position = "top",
        legend.key.height = unit(9, "pt"), legend.key.width = unit(34, "pt"),
        legend.title = element_text(size = 8.5), legend.text = element_text(size = 8),
        plot.title = element_text(face = "bold", size = 11.5),
        plot.subtitle = element_text(size = 8.2, colour = "grey30", lineheight = 1.35,
                                     margin = margin(b = 8)),
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold", size = 9.5),
        strip.text.y.left = element_text(face = "bold", size = 8.6, angle = 0,
                                         hjust = 1, lineheight = 1.1),
        strip.placement = "outside",
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 22, hjust = 1),
        axis.title.x = element_text(size = 9.5, margin = margin(t = 7)),
        axis.title.y = element_text(size = 9.5, margin = margin(r = 7)),
        axis.line = element_blank(), axis.ticks = element_blank(),
        panel.spacing.x = unit(7, "pt"), panel.spacing.y = unit(9, "pt"))

ggsave(file.path(OUT_FIG, "supp_figS2_confusion_matrix.pdf"), p,
       width = 7.3, height = 9.6, device = cairo_pdf)
ggsave(file.path(OUT_FIG, "supp_figS2_confusion_matrix.png"), p,
       width = 7.3, height = 9.6, dpi = 300)

acc <- cm %>% group_by(Scenario, Criteria) %>%
  summarise(correct = 100 * sum(n[correct]) / sum(n), .groups = "drop") %>%
  pivot_wider(names_from = Criteria, values_from = correct)
MED <- "Mediation (proposed)"; TRAD <- "Traditional (ORR)"
g <- function(s, col) round(acc[[col]][acc$Scenario == s], 1)

stopifnot(
  g("Null", MED) == 92.6,          g("Null", TRAD) == 93.2,
  g("Survival-Null", MED) == 89.4, g("Survival-Null", TRAD) == 5.8,
  g("Moderate-Local", MED) == 92.2,g("Moderate-Local", TRAD) == 70.2,
  g("Mediated-Only", MED) == 92.6, g("Mediated-Only", TRAD) == 92.6,
  g("Direct-Only", MED) == 96.4,   g("Direct-Only", TRAD) == 3.8,
  g("Mixed", MED) == 97.8,         g("Mixed", TRAD) == 98.0)
cat("[08] verification PASSED against Table 4; confusion matrix written\n")
print(as.data.frame(acc), row.names = FALSE)
