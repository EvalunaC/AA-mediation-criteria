args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("usage: Rscript worker_confound.R <scenario> <u_index> <rep_from> <rep_to>")
SC <- args[1]; UI <- as.integer(args[2]); RF <- as.integer(args[3]); RT <- as.integer(args[4])

source("config.R")
stopifnot(SC %in% SCENARIOS_C, UI >= 1, UI <= length(U_GRID))
U  <- U_GRID[UI]
SI <- match(SC, SCENARIOS_C)
dir.create(CONF_RAW, showWarnings = FALSE, recursive = TRUE)
msm_load()

PARM <- confound_parm(scen_parm(SC), U)
stopifnot(file.exists(MODEL_FILE))

inc_for <- function(rhs) {
  terms <- trimws(strsplit(rhs, "\\+")[[1]])
  list(include12 = terms, include13 = terms, include23 = terms)
}

fit_one <- function(idat, rhs) {
  inc <- inc_for(rhs)
  for (att in seq_len(FIT_TRIES_C)) {
    z <- tryCatch({
      f <- FitJagsMstate(as.formula(paste("Surv(TTP,TTPevent)~", rhs)),
                         as.formula(paste("Surv(OS,OSevent)~", rhs)),
                         as.formula(paste("~", rhs)), idat,
                         modelFile = MODEL_FILE, include = inc, mcmc.par = MCMC)
      logRRCalc(f, tau = TAU_C, arm.name = "trt", med.name = "resp")
    }, error = function(e) NULL)
    if (!is.null(z)) {
      m <- as.numeric(z$Mediation[, , 3]); d <- as.numeric(z$Mediation[, , 2])
      te <- as.numeric(z$Mediation[, , 1])
      if (length(m) > 0 && !all(is.na(m)))
        return(list(ok = TRUE, p_nie = mean(m > 0), p_nde = mean(d > 0),
                    p_te = mean(te > 0, na.rm = TRUE),
                    lrr_nie = median(m), lrr_nde = median(d), lrr_te = median(te),
                    lo = quantile(m, .025), hi = quantile(m, .975)))
    }
  }
  list(ok = FALSE, p_nie = NA, p_nde = NA, p_te = NA,
       lrr_nie = NA, lrr_nde = NA, lrr_te = NA, lo = NA, hi = NA)
}

for (rep in RF:RT) {
  outf <- file.path(CONF_RAW, sprintf("%s_u%d_rep%04d.rds", SC, UI, rep))
  if (file.exists(outf)) next
  t0 <- Sys.time()
  set.seed(conf_seed(SI, UI, rep))
  dat <- gen_mstate(parm = PARM, n = N_PATIENTS, nb = 1, nc = 1, prtrt = PRTRT,
                    accTime = ACC_TIME, accExN = NULL, addFUt = ADD_FU, CenUpLim = NULL)
  lr_full <- tryCatch(logrank_z(dat$OS, dat$OSevent, dat$trt, "two.sided"),
                      error = function(e) NA_real_)
  idat <- interim_mstate(dat = dat, tp = INTERIM_TP_C, delt = 0); idat <- idat[idat$enrolled, ]
  lr_int <- tryCatch(logrank_z(idat$PFS, idat$PFSevent, idat$trt, "two.sided"),
                     error = function(e) NA_real_)
  chi <- tryCatch(suppressWarnings(chisq.test(idat$resp, idat$trt)$p.value),
                  error = function(e) NA_real_)

  om <- fit_one(idat, RHS_OMIT)
  ad <- fit_one(idat, RHS_ADJ)

  res <- data.frame(
    scenario = SC, u = U, u_index = UI, rep = rep,
    logrank_full = lr_full, logrank_interim = lr_int, chisq_p = chi,
    n_enrolled = nrow(idat), os_events_int = sum(idat$OSevent),
    orr_int = mean(idat$resp),
    resp_x1_0 = mean(idat$resp[idat$x1 == 0]), resp_x1_1 = mean(idat$resp[idat$x1 == 1]),
    osev_x1_0 = mean(idat$OSevent[idat$x1 == 0]), osev_x1_1 = mean(idat$OSevent[idat$x1 == 1]),
    p_nie_omit = om$p_nie, p_nde_omit = om$p_nde, p_te_omit = om$p_te,
    lrr_nie_omit = om$lrr_nie, lrr_nde_omit = om$lrr_nde, lrr_te_omit = om$lrr_te,
    p_nie_adj = ad$p_nie, p_nde_adj = ad$p_nde, p_te_adj = ad$p_te,
    lrr_nie_adj = ad$lrr_nie, lrr_nde_adj = ad$lrr_nde, lrr_te_adj = ad$lrr_te,
    nie_lo_omit = om$lo, nie_hi_omit = om$hi,
    nie_lo_adj = ad$lo, nie_hi_adj = ad$hi,
    ok_omit = om$ok, ok_adj = ad$ok, stringsAsFactors = FALSE)
  saveRDS(res, outf)
  cat(sprintf("[%s u=%.2f rep %d] %.0f s | omit ok=%s | adj ok=%s\n",
      SC, U, rep, as.numeric(difftime(Sys.time(), t0, units = "secs")),
      om$ok, ad$ok))
  flush(stdout())
}
cat(sprintf("WORKER DONE %s u%d %d-%d\n", SC, UI, RF, RT))
