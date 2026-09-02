args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("usage: Rscript worker_maturity.R <scenario> <rep_from> <rep_to>")
SC <- args[1]; REP_FROM <- as.integer(args[2]); REP_TO <- as.integer(args[3])

source("config.R")
stopifnot(SC %in% SCENARIOS)
SI <- match(SC, SCENARIOS)
dir.create(MAT_RAW, showWarnings = FALSE, recursive = TRUE)
msm_load()

parm <- scen_parm(SC)
f12 <- Surv(TTP, TTPevent) ~ trt + resp
f13 <- Surv(OS, OSevent) ~ trt + resp
f23 <- ~ trt + resp
stopifnot(file.exists(MODEL_FILE))

fit_cut <- function(idat, tp) {
  taus <- unique(c(tp, TAU_FIXED))
  for (att in seq_len(FIT_MAX_TRIES)) {
    lrr <- tryCatch({
      fit <- FitJagsMstate(f12, f13, f23, idat, modelFile = MODEL_FILE,
                           include = INCLUDE_TRT, mcmc.par = MCMC)
      logRRCalc(fit, tau = taus, arm.name = "trt", med.name = "resp")
    }, error = function(e) NULL)
    if (!is.null(lrr)) {
      i_mv <- match(tp, taus); i_fx <- match(TAU_FIXED, taus)
      m_mv <- as.numeric(lrr$Mediation[i_mv, , 3]); d_mv <- as.numeric(lrr$Mediation[i_mv, , 2])
      m_fx <- as.numeric(lrr$Mediation[i_fx, , 3]); d_fx <- as.numeric(lrr$Mediation[i_fx, , 2])
      t_mv <- as.numeric(lrr$Mediation[i_mv, , 1]); t_fx <- as.numeric(lrr$Mediation[i_fx, , 1])
      if (length(m_mv) > 0 && !all(is.na(m_mv)))
        return(list(ok = TRUE, tries = att,
          p_nie_mv = mean(m_mv > 0, na.rm = TRUE), p_nde_mv = mean(d_mv > 0, na.rm = TRUE),
          p_nie_fx = mean(m_fx > 0, na.rm = TRUE), p_nde_fx = mean(d_fx > 0, na.rm = TRUE),
          p_te_mv  = mean(t_mv > 0, na.rm = TRUE), p_te_fx  = mean(t_fx > 0, na.rm = TRUE),
          med_lRR_nie = median(m_mv, na.rm = TRUE), med_lRR_nde = median(d_mv, na.rm = TRUE),
          n_draws = length(m_mv)))
    }
  }
  list(ok = FALSE, tries = FIT_MAX_TRIES, p_nie_mv = NA, p_nde_mv = NA,
       p_nie_fx = NA, p_nde_fx = NA, p_te_mv = NA, p_te_fx = NA,
       med_lRR_nie = NA, med_lRR_nde = NA, n_draws = 0L)
}

for (rep in REP_FROM:REP_TO) {
  outf <- file.path(MAT_RAW, sprintf("%s_rep%04d.rds", SC, rep))
  if (file.exists(outf)) next
  t_rep <- Sys.time()

  set.seed(mat_seed(SI, rep))
  dat <- gen_mstate(parm = parm, n = N_PATIENTS, nb = 1, nc = 1, prtrt = PRTRT,
                    accTime = ACC_TIME, accExN = NULL, addFUt = ADD_FU, CenUpLim = NULL)
  logrank_full <- tryCatch(logrank_z(dat$OS, dat$OSevent, dat$trt, "two.sided"),
                           error = function(e) NA_real_)
  os_events_fin <- sum(dat$OSevent)

  out <- list()
  for (tp in INTERIM_TPS) {
    idat <- interim_mstate(dat = dat, tp = tp, delt = 0); idat <- idat[idat$enrolled, ]
    lr_int <- tryCatch(logrank_z(idat$PFS, idat$PFSevent, idat$trt, "two.sided"),
                       error = function(e) NA_real_)
    chi <- tryCatch(suppressWarnings(chisq.test(idat$resp, idat$trt)$p.value),
                    error = function(e) NA_real_)
    fr <- fit_cut(idat, tp)
    out[[length(out) + 1]] <- data.frame(
      scenario = SC, rep = rep, tp = tp,
      maturity_level = sprintf("Month %d", tp),
      logrank_full = logrank_full, os_events_fin = os_events_fin,
      n_enrolled = nrow(idat), os_events_int = sum(idat$OSevent),
      pfs_events_int = sum(idat$PFSevent),
      maturity = sum(idat$OSevent) / os_events_fin,
      orr_int = mean(idat$resp),
      logrank_interim = lr_int, chisq_p = chi,
      p_nie_mv = fr$p_nie_mv, p_nde_mv = fr$p_nde_mv,
      p_nie_fx = fr$p_nie_fx, p_nde_fx = fr$p_nde_fx,
      p_te_mv = fr$p_te_mv, p_te_fx = fr$p_te_fx,
      med_lRR_nie = fr$med_lRR_nie, med_lRR_nde = fr$med_lRR_nde,
      converged = fr$ok, n_tries = fr$tries, n_draws = fr$n_draws,
      stringsAsFactors = FALSE)
  }
  res <- do.call(rbind, out)
  attr(res, "elapsed_sec") <- as.numeric(difftime(Sys.time(), t_rep, units = "secs"))
  saveRDS(res, outf)
  cat(sprintf("[%s rep %d] %.0f s | converged %d/%d\n", SC, rep,
      attr(res, "elapsed_sec"), sum(res$converged), nrow(res)))
  flush(stdout())
}
cat(sprintf("WORKER DONE %s %d-%d\n", SC, REP_FROM, REP_TO))
