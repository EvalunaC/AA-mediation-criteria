source("00_criteria.R")
suppressMessages({library(dplyr); library(tidyr); library(ggplot2)})

COHORTS <- c("309 all", "309 KRAS wt")
CLAB <- c(`309 all` = "All randomized (n = 935)",
          `309 KRAS wt` = "KRAS wild-type (n = 514)")
CUT_ORD <- c("6mo","9mo","12mo","18mo","24mo","final")
CUT_LAB <- c("6","9","12","18","24","Final")

files <- list.files(CASE_RAW, pattern = "\\.rds$", full.names = TRUE)
cells <- lapply(files, readRDS)
grid <- bind_rows(lapply(cells, function(x) {
  if (!is.data.frame(x$eff)) return(NULL)
  pri <- x$eff[x$eff$is_primary, , drop = FALSE]
  if (!nrow(pri)) return(NULL)
  cbind(x$meta, pri[1, c("RR_TE","TE_lo","TE_hi","P_TE",
                         "RR_NDE","NDE_lo","NDE_hi","P_NDE",
                         "RR_NIE","NIE_lo","NIE_hi","P_NIE")])
})) %>% filter(cohort %in% COHORTS) %>%
  mutate(Cohort = factor(CLAB[cohort], levels = CLAB),
         cutf = factor(cut, levels = CUT_ORD, labels = CUT_LAB))

CCOL <- c("#2f7d8f", "#c9962e"); names(CCOL) <- CLAB

g12 <- grid %>% filter(cut == "12mo")
fa <- bind_rows(
  g12 %>% transmute(Cohort, Effect = "Indirect (NIE)", RR = RR_NIE, lo = NIE_lo, hi = NIE_hi),
  g12 %>% transmute(Cohort, Effect = "Direct (NDE)",   RR = RR_NDE, lo = NDE_lo, hi = NDE_hi),
  g12 %>% transmute(Cohort, Effect = "Total (TE)",     RR = RR_TE,  lo = TE_lo,  hi = TE_hi)) %>%
  mutate(Effect = factor(Effect, levels = c("Indirect (NIE)", "Direct (NDE)", "Total (TE)")))
pa <- ggplot(fa, aes(RR, Effect, colour = Cohort)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_pointrange(aes(xmin = lo, xmax = hi),
                  position = position_dodge(width = .5), size = .45) +
  scale_colour_manual(values = CCOL) +
  labs(x = "Posterior risk ratio for survival at the interim horizon (95% CrI)",
       y = NULL, colour = NULL,
       title = "(a) Effect decomposition at the specified 12-month interim") +
  theme_classic(base_size = 10) +
  theme(legend.position = "top", plot.title = element_text(size = 10, face = "bold"))

pb <- ggplot(grid, aes(cutf, P_TE, colour = Cohort, group = Cohort)) +
  geom_hline(yintercept = TE_GATE, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = .6) + geom_point(size = 2.2) +
  annotate("text", x = 0.7, y = TE_GATE + 0.045, label = "TE gate (0.5)",
           hjust = 0, size = 3, colour = "grey30") +
  scale_colour_manual(values = CCOL) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Interim analysis time (months; Final = complete follow-up)",
       y = "Pr(RR_TE > 1 | data)",
       colour = NULL,
       title = "(b) Total-effect coherence across interim timings") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none", plot.title = element_text(size = 10, face = "bold"))

p <- ggpubr_fallback <- NULL
suppressMessages({
  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork); p <- pa / pb
  } else {
    library(gridExtra); p <- arrangeGrob(pa, pb, ncol = 1)
  }
})
ggsave(file.path(OUT_FIG, "fig6_case_study.pdf"), p,
       width = 6.8, height = 6.6, device = cairo_pdf)
ggsave(file.path(OUT_FIG, "fig6_case_study.png"), p,
       width = 6.8, height = 6.6, dpi = 300)

a <- grid %>% filter(cohort == "309 all", cut == "12mo")
stopifnot(round(a$RR_TE, 3) == 0.922, round(a$P_TE, 3) == 0.023,
          round(a$RR_NIE, 3) == 1.041)
cat("[07] verification PASSED; case-study figure written\n")
