args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("usage: Rscript worker_sample_size.R <n_target> <scenario> <rep_from> <rep_to>")
N_TARGET <- as.integer(args[1]); SC <- args[2]
REP_FROM <- as.integer(args[3]); REP_TO <- as.integer(args[4])
stopifnot(N_TARGET %in% c(100L, 200L, 300L))

source("config.R")
stopifnot(SC %in% SCENARIOS)
SI <- match(SC, SCENARIOS)
dir.create(SS_RAW, showWarnings = FALSE, recursive = TRUE)
msm_load()

TP <- 24
parm <- scen_parm(SC)
f12 <- Surv(TTP, TTPevent) ~ trt + resp
f13 <- Surv(OS, OSevent) ~ trt + resp
f23 <- ~ trt + resp
stopifnot(file.exists(MODEL_FILE))

fit_cut <- function(idat, tp) {
  for (att in seq_len(FIT_MAX_TRIES)) {
    lrr <- tryCatch({
      fit <- FitJagsMstate(f12, f13, f23, idat, modelFile = MODEL_FILE,
                           include = INCLUDE_TRT, mcmc.par = MCMC)
      logRRCalc(fit, tau = tp, arm.name = "trt", med.name = "resp")
    }, error = function(e) NULL)
    if (!is.null(lrr)) {
      m <- as.numeric(lrr$Mediation[1, , 3])
      d <- as.numeric(lrr$Mediation[1, , 2])
      t <- as.numeric(lrr$Mediation[1, , 1])
      if (length(m) > 0 && !all(is.na(m)))
        return(list(ok = TRUE, tries = att,
          p_nie_mv = mean(m > 0, na.rm = TRUE),
          p_nde_mv = mean(d > 0, na.rm = TRUE),
          p_te_mv  = mean(t > 0, na.rm = TRUE),
          med_lRR_nie = median(m, na.rm = TRUE),
          med_lRR_nde = median(d, na.rm = TRUE),
          n_draws = length(m)))
    }
  }
  list(ok = FALSE, tries = FIT_MAX_TRIES, p_nie_mv = NA, p_nde_mv = NA,
       p_te_mv = NA, med_lRR_nie = NA, med_lRR_nde = NA, n_draws = 0L)
}

for (rep in REP_FROM:REP_TO) {
  outf <- file.path(SS_RAW, sprintf("n%03d_%s_rep%04d.rds", N_TARGET, SC, rep))
  if (file.exists(outf)) next
  t_rep <- Sys.time()

  set.seed(ss_seed(N_TARGET, SI, rep))
  dat <- gen_mstate(parm = parm, n = N_TARGET, nb = 1, nc = 1, prtrt = PRTRT,
                    accTime = ACC_TIME, accExN = NULL, addFUt = ADD_FU, CenUpLim = NULL)
  logrank_full <- tryCatch(logrank_z(dat$OS, dat$OSevent, dat$trt, "two.sided"),
                           error = function(e) NA_real_)
  os_events_fin <- sum(dat$OSevent)

  idat <- interim_mstate(dat = dat, tp = TP, delt = 0); idat <- idat[idat$enrolled, ]
  lr_int <- tryCatch(logrank_z(idat$PFS, idat$PFSevent, idat$trt, "two.sided"),
                     error = function(e) NA_real_)
  chi <- tryCatch(suppressWarnings(chisq.test(idat$resp, idat$trt)$p.value),
                  error = function(e) NA_real_)
  fr <- fit_cut(idat, TP)

  res <- data.frame(
    n_target = N_TARGET, scenario = SC, rep = rep, tp = TP,
    logrank_full = logrank_full, os_events_fin = os_events_fin,
    n_enrolled = nrow(idat), os_events_int = sum(idat$OSevent),
    pfs_events_int = sum(idat$PFSevent),
    maturity = sum(idat$OSevent) / os_events_fin,
    orr_int = mean(idat$resp),
    logrank_interim = lr_int, chisq_p = chi,
    p_nie_mv = fr$p_nie_mv, p_nde_mv = fr$p_nde_mv, p_te_mv = fr$p_te_mv,
    med_lRR_nie = fr$med_lRR_nie, med_lRR_nde = fr$med_lRR_nde,
    converged = fr$ok, n_tries = fr$tries, n_draws = fr$n_draws,
    stringsAsFactors = FALSE)
  attr(res, "elapsed_sec") <- as.numeric(difftime(Sys.time(), t_rep, units = "secs"))
  saveRDS(res, outf)
  cat(sprintf("[n=%d %s rep %d] %.0f s (converged=%s)\n",
              N_TARGET, SC, rep, attr(res, "elapsed_sec"), fr$ok))
  flush(stdout())
}
