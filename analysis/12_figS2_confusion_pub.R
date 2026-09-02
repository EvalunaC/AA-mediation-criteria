source("00_criteria.R")
suppressMessages(library(ggplot2))

TEAL <- "#2f7d8f"
INK  <- "#16181c"

d <- read.csv(file.path(OUT_TAB, "supp_tableS9_confusion_matrix.csv"),
              stringsAsFactors = FALSE, check.names = FALSE)
ord <- function(lab, key) unique(d[[lab]][order(d[[key]])])
d$Decision <- factor(d$Decision, levels = ord("Decision", "decision_order"))
d$Gold     <- factor(d$Gold,     levels = ord("Gold",     "gold_order"))
d$Criteria <- factor(d$Criteria, levels = ord("Criteria", "criteria_order"))
d$Scenario <- factor(d$Scenario, levels = ord("Scenario", "scenario_order"))

NMAX <- max(d$n)
p <- ggplot(d, aes(Decision, Gold)) +
  geom_tile(aes(fill = n), colour = "white", linewidth = 0.7) +
  geom_tile(data = d[d$correct, ], fill = NA, colour = INK, linewidth = 0.7) +
  geom_text(aes(label = n, colour = n > 0.55 * NMAX), size = 2.9) +
  facet_grid(Scenario ~ Criteria, scales = "free_x", space = "free_x",
             switch = "y") +
  scale_fill_gradient(low = "#ffffff", high = TEAL, limits = c(0, NMAX),
                      name = "Trials (n)") +
  scale_colour_manual(values = c(`FALSE` = INK, `TRUE` = "#ffffff"),
                      guide = "none") +
  scale_y_discrete(limits = rev(levels(d$Gold))) +
  labs(x = "Predicted outcome", y = "Actual outcome") +
  theme_classic(base_size = 11) +
  theme(legend.position = "top",
        legend.key.height = unit(9, "pt"), legend.key.width = unit(34, "pt"),
        legend.title = element_text(size = 8.5), legend.text = element_text(size = 8),
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

ggsave(file.path(OUT_FIG, "FigureS2_confusion_matrix.pdf"), p,
       width = 7.3, height = 8.85, device = cairo_pdf)
ggsave(file.path(OUT_FIG, "FigureS2_confusion_matrix.png"), p,
       width = 7.3, height = 8.85, dpi = 300)

num <- aggregate(n ~ Criteria + Scenario, data = d[d$correct, ], FUN = sum)
den <- aggregate(n ~ Criteria + Scenario, data = d, FUN = sum)
a <- merge(num, den, by = c("Criteria", "Scenario"), suffixes = c("_ok", "_all"))
a$acc <- round(100 * a$n_ok / a$n_all, 1)
g <- function(s, cr) a$acc[a$Scenario == s & a$Criteria == cr]
MED <- "Mediation (proposed)"; TRAD <- "Traditional (ORR)"
stopifnot(
  g("Null", MED) == 92.6,          g("Null", TRAD) == 93.2,
  g("Survival-Null", MED) == 89.4, g("Survival-Null", TRAD) == 5.8,
  g("Moderate-Local", MED) == 92.2,g("Moderate-Local", TRAD) == 70.2,
  g("Mediated-Only", MED) == 92.6, g("Mediated-Only", TRAD) == 92.6,
  g("Direct-Only", MED) == 96.4,   g("Direct-Only", TRAD) == 3.8,
  g("Mixed", MED) == 97.8,         g("Mixed", TRAD) == 98.0)
cat("[12] supplement figure S2 (confusion matrix, publication) written\n")
