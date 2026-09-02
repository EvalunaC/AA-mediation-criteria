args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: Rscript worker_case.R \"<cohort key>\" <t|final>")
key <- args[1]
t   <- if (identical(args[2], "final")) Inf else as.numeric(args[2])

source("config.R")
source("cohorts_case.R")
msm_load()

dir.create(CASE_RAW, showWarnings = FALSE, recursive = TRUE)
out <- file.path(CASE_RAW, sprintf("%s__%s.rds", gsub("[^A-Za-z0-9]+", "_", key), cut_lab(t)))
if (file.exists(out)) { cat("skip (exists):", out, "\n"); quit(save = "no") }

spec <- COHORTS_CASE[[key]]
full <- case_data(key)
di   <- case_cut(full, t, spec$resp_def)

lr <- function(d, tc, ec) {
  if (sum(d[[ec]] == 1) < 5 || length(unique(d$trt)) < 2) return(NA_real_)
  s <- try(survdiff(Surv(d[[tc]], d[[ec]] == 1) ~ d$trt), silent = TRUE)
  if (inherits(s, "try-error")) return(NA_real_)
  1 - pchisq(s$chisq, 1)
}
orrp <- function(d) {
  tb <- table(d$resp, d$trt)
  if (any(dim(tb) < 2)) return(NA_real_)
  suppressWarnings(chisq.test(tb)$p.value)
}

TAU_FIX <- c(6, 12, 18, 24, 30)
tau_pri <- if (is.infinite(t)) 24 else t
TAUS    <- sort(unique(c(TAU_FIX, tau_pri)))

meta <- data.frame(
  cohort = key, label = spec$label, cut = cut_lab(t), t = t,
  resp_def = spec$resp_def, tau_primary = tau_pri,
  n = nrow(di), n0 = sum(di$trt == 0), n1 = sum(di$trt == 1),
  resp0 = sum(di$resp == 1 & di$trt == 0), resp1 = sum(di$resp == 1 & di$trt == 1),
  ORR0 = mean(di$resp[di$trt == 0]), ORR1 = mean(di$resp[di$trt == 1]),
  p_ORR_int   = orrp(di),
  os_ev_int   = sum(di$OSevent == 1),
  os_ev_final = sum(full$OSevent == 1),
  maturity    = sum(di$OSevent == 1) / sum(full$OSevent == 1),
  p_PFS_int   = lr(di, "TTP", "TTPevent"),
  p_OS_int    = lr(di, "OS", "OSevent"),
  p_ORR_final = orrp(full),
  p_OS_final  = lr(full, "OS", "OSevent"),
  stringsAsFactors = FALSE)

cat(sprintf("=== %s | cut %s | n=%d | resp %d/%d | maturity %.0f%% ===\n",
            key, cut_lab(t), meta$n, meta$resp0, meta$resp1, 100 * meta$maturity))
flush(stdout())

fit <- NULL; res <- NULL; err <- NA_character_
.ident <- (meta$resp0 > 0 && meta$resp0 < meta$n0) || (meta$resp1 > 0 && meta$resp1 < meta$n1)
if ((meta$resp0 + meta$resp1) >= 10 && .ident && meta$os_ev_int >= 10) {
  set.seed(BASE_SEED + as.integer(factor(key, levels = names(COHORTS_CASE))) * 1e5 +
             ifelse(is.infinite(t), 999, t))
  for (try_i in seq_len(FIT_MAX_TRIES)) {
    f <- try(FitJagsMstate(as.formula("Surv(TTP,TTPevent)~trt+resp"),
                           as.formula("Surv(OS,OSevent)~trt+resp"),
                           as.formula("~trt+resp"), di,
                           modelFile = MODEL_FILE,
                           include = INCLUDE_TRT, mcmc.par = MCMC), silent = TRUE)
    if (!inherits(f, "try-error")) {
      r <- try(logRRCalc(f, tau = TAUS, arm.name = "trt", med.name = "resp"), silent = TRUE)
      if (!inherits(r, "try-error")) { fit <- f; res <- r; break }
      err <- as.character(r)
    } else err <- as.character(f)
    cat("  refit attempt", try_i, "failed\n"); flush(stdout())
  }
} else err <- sprintf("not fittable: resp %d/%d (identifiable=%s), interim OS events %d",
                      meta$resp0, meta$resp1, .ident, meta$os_ev_int)

eff <- NULL
if (!is.null(res)) {
  q <- function(v) { n0 <- length(v); v <- v[is.finite(v)]
    if (length(v) < 0.5 * n0) return(c(NA, NA, NA, NA, length(v) / n0))
    c(exp(median(v)), exp(quantile(v, .025)), exp(quantile(v, .975)),
      mean(v > 0), length(v) / n0) }
  eff <- do.call(rbind, lapply(seq_along(TAUS), function(i) {
    g <- function(k) as.numeric(res$Mediation[i, , k])
    tot <- q(g(1)); nde <- q(g(2)); nie <- q(g(3))
    data.frame(cohort = key, cut = cut_lab(t), tau = TAUS[i],
               is_primary = TAUS[i] == tau_pri,
               RR_TE = tot[1],  TE_lo = tot[2],  TE_hi = tot[3],  P_TE = tot[4],
               RR_NDE = nde[1], NDE_lo = nde[2], NDE_hi = nde[3], P_NDE = nde[4],
               RR_NIE = nie[1], NIE_lo = nie[2], NIE_hi = nie[3], P_NIE = nie[4],
               finite_TE = tot[5], finite_NDE = nde[5], finite_NIE = nie[5],
               stringsAsFactors = FALSE)
  }))
}
saveRDS(list(meta = meta, eff = eff, err = err, taus = TAUS), out)
cat(if (is.null(eff)) paste("FAILED:", err, "\n") else "ok\n")
