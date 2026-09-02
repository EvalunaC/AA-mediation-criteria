source("00_criteria.R")
suppressMessages({library(dplyr); library(tidyr); library(ggplot2)})

t4 <- read.csv(file.path(OUT_TAB, "table4_oc_by_maturity_scenario.csv"))
SC_ORD3 <- c("Null", "Survival-Null", "Moderate-Local",
             "Mediated-Only", "Direct-Only", "Mixed", "Mean of 6 scenarios")

w <- t4 %>%
  mutate(month = as.integer(sub("Month ", "", Maturity))) %>%
  select(month, Scenario, med = Acc_med, trad = Acc_trad)
w <- bind_rows(w, w %>% group_by(month) %>%
                 summarise(med = round(mean(med), 1),
                           trad = round(mean(trad), 1), .groups = "drop") %>%
                 mutate(Scenario = "Mean of 6 scenarios")) %>%
  mutate(Scenario = factor(Scenario, levels = SC_ORD3))

d <- w %>%
  pivot_longer(c(med, trad), names_to = "crit", values_to = "rate") %>%
  mutate(Criteria = factor(ifelse(crit == "med", "Mediation (proposed)",
                                  "Traditional (ORR)"), levels = names(COL2)))

lab <- d %>%
  left_join(w %>% transmute(month, Scenario, gap = abs(med - trad)),
            by = c("month", "Scenario")) %>%
  filter(gap >= 5) %>%
  mutate(y = rate + ifelse(crit == "med", 5.6, -6.0))

OS_EVT <- c(`12` = "9%", `18` = "30%", `24` = "54%", `30` = "80%")
xlabs <- function(m) sprintf("%d\n%s\n%s", m, OS_EVT[as.character(m)],
                             formatC(DELTA_STAR[as.character(m)], format = "g"))
LEG <- c("Traditional ORR criterion", "Mediation criterion")

p <- ggplot(d, aes(month, rate, colour = Criteria)) +
  geom_segment(data = w, aes(x = month, xend = month, y = trad, yend = med),
               inherit.aes = FALSE, colour = "grey85",
               linewidth = 2, lineend = "round") +
  geom_point(aes(shape = Criteria, size = Criteria)) +
  geom_text(data = lab, aes(y = y, label = sprintf("%.1f", rate)),
            size = 2.3, show.legend = FALSE) +
  facet_wrap(~Scenario, nrow = 2) +
  scale_colour_manual(values = COL2, breaks = names(COL2), labels = LEG) +
  scale_shape_manual(values = setNames(c(16, 18), names(COL2)),
                     breaks = names(COL2), labels = LEG) +
  scale_size_manual(values = setNames(c(2.2, 3), names(COL2)),
                    breaks = names(COL2), labels = LEG) +
  scale_x_continuous(breaks = c(12, 18, 24, 30), labels = xlabs,
                     expand = expansion(add = 2.4)) +
  scale_y_continuous(breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(0.07, 0.09))) +
  coord_cartesian(ylim = c(0, 100), clip = "off") +
  labs(x = paste("Interim analysis:  month",
                 "·  % of final OS events observed",
                 "·  calibrated Δ*", sep = "  "),
       y = "Correct decision rate (%)",
       colour = NULL, shape = NULL, size = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.35),
        strip.text = element_text(face = "bold", hjust = 0, size = 8.5,
                                  colour = "grey15"),
        axis.title.y = element_text(size = 9, colour = "grey20"),
        axis.title.x = element_text(size = 8.5, colour = "grey30",
                                    margin = margin(t = 7)),
        axis.text = element_text(size = 7.3, colour = "grey30"),
        axis.text.x = element_text(lineheight = 1.15),
        legend.position = "bottom",
        legend.text = element_text(size = 8.5),
        panel.spacing.x = unit(16, "pt"), panel.spacing.y = unit(10, "pt"),
        plot.margin = margin(6, 10, 4, 6))

ggsave(file.path(OUT_FIG, "Figure3_correct_decision_by_scenario_maturity.pdf"),
       p, width = 7.6, height = 4.452, device = cairo_pdf)
ggsave(file.path(OUT_FIG, "Figure3_correct_decision_by_scenario_maturity.png"),
       p, width = 7.6, height = 4.452, dpi = 300)

chk <- w %>% filter(month == 24, Scenario == "Survival-Null")
stopifnot(chk$med == 89.4, chk$trad == 5.8)
m30 <- w %>% filter(month == 30, Scenario == "Mean of 6 scenarios")
stopifnot(m30$med == 94.0, m30$trad == 61.2)
cat("[11] main-text figure 3 (dumbbell) written\n")
