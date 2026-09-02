source("00_criteria.R")

tabS1 <- data.frame(
  Scenario   = c("Null", "Survival-Null", "Moderate-Local",
                 "Mediated-Only", "Direct-Only", "Mixed"),
  gamma_a    = c(0, 1.1156, 0.6360, 1.1156, 0, 1.1156),
  ORR_control_to_treated = c("0.15 -> 0.15", "0.15 -> 0.35", "0.15 -> 0.25",
                             "0.15 -> 0.35", "0.15 -> 0.15", "0.15 -> 0.35"),
  gamma_12 = c(0, 0, -0.28,  0,   -0.4,  -0.2),
  eta_12   = c(0, 0, -0.36, -2,    0,    -1),
  gamma_13 = c(0, 0, -0.175, 0,   -0.25, -0.125),
  eta_13   = c(0, 0, -0.36, -2,    0,    -1),
  gamma_23 = c(0, 0, -0.07,  0,   -0.1,  -0.05),
  eta_23   = c(0, 0, -0.09, -0.5,  0,    -0.25))
write.csv(tabS1, file.path(OUT_TAB, "supp_tableS1_scenario_parameters.csv"),
          row.names = FALSE)

tabS2 <- data.frame(
  Component = c("Transition 1 -> 2 (progression)",
                "Transition 1 -> 3 (death without progression)",
                "Transition 2 -> 3 (death after progression)",
                "Time to response"),
  Weibull_shape = c(2, 1.33, 3, 5),
  Weibull_rate_or_scale = c("rate 1/55", "rate 1/600", "rate 1/780", "scale 1.5"))
write.csv(tabS2, file.path(OUT_TAB, "supp_tableS2_baseline_weibull.csv"),
          row.names = FALSE)

cat("[04] tables S1 and S2 written\n")
