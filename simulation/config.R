MSM_DIR    <- Sys.getenv("MSM_DIR", "../MultiStateModels")
MODEL_FILE <- file.path(MSM_DIR, "R", "Models", "modelWeibMstateLik2.txt")
PATCHES    <- c("../patches/patch_fitjags.R", "../patches/patch_logrrcalc_fast.R")

MAT_RAW  <- "../data/maturity"
SS_RAW   <- "../data/sample_size"
CONF_RAW <- "../data/confounding"
CASE_RAW <- "../data/case_interim"

SCENARIOS  <- c("sc_null", "sc_snull", "sc_modl", "sc_med", "sc_direct", "sc_mixed")
N_PATIENTS <- 500
N_SIMS     <- 500
PRTRT      <- 0.5
ACC_TIME   <- 24
ADD_FU     <- 24

INTERIM_TPS <- c(12, 18, 24, 30)
TAU_FIXED   <- 24

MCMC <- list(niter = 3000, nburn = 1000, nchain = 2, nthin = 1)
INCLUDE_TRT <- list(include12 = c("trt", "resp"),
                    include13 = c("trt", "resp"),
                    include23 = c("trt", "resp"))
BASE_SEED     <- 1
FIT_MAX_TRIES <- 4

mat_seed <- function(scen_idx, rep) BASE_SEED + scen_idx * 1000000L + rep * 100L
ss_seed  <- function(n_target, si, rep) 7000000L + si * 100000L + n_target * 10L + rep

U_GRID      <- c(0, 0.20, 0.35, 0.50)
SCENARIOS_C <- c("sc_snull", "sc_med", "sc_null", "sc_modl", "sc_direct", "sc_mixed")
INTERIM_TP_C <- 24
TAU_C        <- 24
N_SIMS_C     <- 500
FIT_TRIES_C  <- 3
RHS_OMIT <- "trt+resp"
RHS_ADJ  <- "trt+resp+x1"

confound_parm <- function(p, u) {
  p$gamma[3] <- u
  for (nm in c("12", "13", "23"))
    p$beta[[paste0("beta", nm)]][3] <- -u
  p
}
conf_seed <- function(si, ui, rep) 7000000L + si * 1000000L + ui * 10000L + rep * 10L

scen_parm <- function(sc) {
  e <- new.env()
  utils::data(list = sc, package = "MultiStateModels", envir = e)
  e[[sc]]
}

msm_load <- function() {
  suppressMessages({
    library(MultiStateModels)
    library(survival)
    library(rjags)
    library(coda)
    library(msm)
    library(HDInterval)
    library(EnvStats)
  })
  ns <- asNamespace("MultiStateModels")
  for (nm in setdiff(ls(ns), getNamespaceExports("MultiStateModels")))
    assign(nm, get(nm, envir = ns), envir = globalenv())
  for (p in PATCHES) capture.output(source(p))
  invisible(TRUE)
}
