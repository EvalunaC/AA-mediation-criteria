t0 <- Sys.time()
for (s in c("01_simulation_results.R", "02_case_study.R",
            "03_confounding.R", "04_design_tables.R",
            "05_calibration_figures.R", "06_sample_size_analysis.R",
            "07_case_figure.R", "08_confusion_matrix.R",
            "09_fig3_scenario_maturity.R", "10_pub_figures.R",
            "11_fig3_dumbbell.R", "12_figS2_confusion_pub.R")) {
  cat("\n========", s, "========\n")
  source(s, echo = FALSE)
  rm(list = setdiff(ls(), c("s", "t0")))
}
cat(sprintf("\nAll outputs reproduced and verified in %.1f min.\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
