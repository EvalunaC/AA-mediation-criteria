source("00_criteria.R")
suppressMessages({library(ggplot2)})

PUB <- file.path(OUT_FIG, "pub"); dir.create(PUB, showWarnings = FALSE)
strip <- function(p) p + labs(title = NULL, subtitle = NULL, caption = NULL)

grab <- function(script, obj, out, w, h) {
  e <- new.env(parent = globalenv())
  suppressWarnings(suppressMessages(sys.source(script, envir = e)))
  p <- strip(get(obj, envir = e))
  ggsave(file.path(PUB, paste0(out, ".png")), p, width = w, height = h, dpi = 300)
  ggsave(file.path(PUB, paste0(out, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  cat(sprintf("[10] %-38s %.1f x %.1f in\n", paste0(out, ".png"), w, h))
}

grab("05_calibration_figures.R", "p5",  "fig_calibration",     10.5, 5.6)
grab("06_sample_size_analysis.R", "pS1","fig_sample_size",     10.2, 4.6)
grab("07_case_figure.R",           "p",  "fig_case_study",       6.8,  6.6)
grab("08_confusion_matrix.R",      "p",  "fig_confusion_matrix", 7.3,  9.0)
cat("[10] publication figures written to", PUB, "\n")
