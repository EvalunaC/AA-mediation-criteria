source("00_criteria.R")
suppressMessages({library(dplyr); library(tidyr); library(ggplot2)})

t4 <- read.csv(file.path(OUT_TAB, "table4_oc_by_maturity_scenario.csv"))
SC_ORD2 <- c("Null", "Survival-Null", "Moderate-Local",
             "Mediated-Only", "Direct-Only", "Mixed")
d <- t4 %>%
  mutate(month = as.integer(sub("Month ", "", Maturity)),
         Scenario = factor(Scenario, levels = SC_ORD2)) %>%
  select(month, Scenario, `Mediation (proposed)` = Acc_med,
         `Traditional (ORR)` = Acc_trad) %>%
  pivot_longer(-c(month, Scenario), names_to = "Criteria", values_to = "rate") %>%
  mutate(Criteria = factor(Criteria, levels = names(COL2)))

p <- ggplot(d, aes(month, rate, colour = Criteria)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.9) +
  facet_wrap(~Scenario, nrow = 2) +
  scale_colour_manual(values = COL2) +
  scale_x_continuous(breaks = c(12, 18, 24, 30)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(x = "Interim analysis time (months from trial start)",
       y = "Correct decision rate (%)", colour = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey93"))

ggsave(file.path(OUT_FIG, "fig3_scenario_maturity.pdf"), p,
       width = 7.6, height = 4.452, device = cairo_pdf)
ggsave(file.path(OUT_FIG, "fig3_scenario_maturity.png"), p,
       width = 7.6, height = 4.452, dpi = 300)

chk <- d %>% filter(month == 24, Scenario == "Survival-Null")
stopifnot(chk$rate[chk$Criteria == "Mediation (proposed)"] == 89.4,
          chk$rate[chk$Criteria == "Traditional (ORR)"] == 5.8)
cat("[09] figure 3 (scenario x maturity) written\n")
