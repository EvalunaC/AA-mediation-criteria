ALPHA_BUDGET <- 0.05
DELTA_STAR <- c(`12` = 0.825, `18` = 0.835, `24` = 0.830, `30` = 0.775)
delta_star <- function(tp) unname(DELTA_STAR[as.character(tp)])
TE_GATE   <- 0.50
ORR_ALPHA <- 0.05

mediation_criteria <- function(p_nie, p_te, p_nde,
                               delta, te_gate = TE_GATE) {
  grant    <- !is.na(p_nie) & p_nie > delta & !is.na(p_te) & p_te > te_gate
  continue <- !grant & !is.na(p_nde) & p_nde > delta
  factor(ifelse(grant, "grant", ifelse(continue, "continue", "stop")),
         levels = c("grant", "continue", "stop"))
}

traditional_criteria <- function(p_orr, alpha = ORR_ALPHA) {
  factor(ifelse(!is.na(p_orr) & p_orr < alpha, "grant", "decline"),
         levels = c("grant", "decline"))
}

vindicated_correct <- function(decision, gold_positive) {
  ifelse(decision %in% c("grant", "continue"), gold_positive, !gold_positive)
}

SCEN_LAB <- c(sc_null = "Null", sc_snull = "Survival-Null",
              sc_modl = "Moderate-Local", sc_med = "Mediated-Only",
              sc_direct = "Direct-Only", sc_mixed = "Mixed")
SCEN_ORD <- c("sc_null", "sc_snull", "sc_modl", "sc_med", "sc_direct", "sc_mixed")
SCEN_TRUTH <- c(sc_null = "No benefit", sc_snull = "No benefit",
                sc_modl = "Mediated benefit", sc_med = "Mediated benefit",
                sc_direct = "Direct-only benefit", sc_mixed = "Mediated benefit")
COL2 <- c("Traditional (ORR)" = "#c9962e", "Mediation (proposed)" = "#2f7d8f")

BASE     <- ".."
SIM_RAW  <- file.path(BASE, "data", "maturity")
CASE_RAW <- file.path(BASE, "data", "case_interim")
CONF_RAW <- file.path(BASE, "data", "confounding")
OUT_FIG  <- file.path(BASE, "output", "figures")
OUT_TAB  <- file.path(BASE, "output", "tables")
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_TAB, showWarnings = FALSE, recursive = TRUE)
