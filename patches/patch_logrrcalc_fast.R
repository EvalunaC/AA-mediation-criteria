.patch_logrrcalc_fast <- function(pkg = MSM_DIR) {
  f <- file.path(pkg, "R/Mstate_JagsFuncs.R")
  src <- readLines(f, warn = FALSE)

  i0 <- grep("^logRRCalc<-function", src)[1]
  if (is.na(i0)) stop("logRRCalc not found in ", f)
  i1 <- grep("^[A-Za-z_.]+[[:space:]]*<-[[:space:]]*function", src)
  i1 <- min(i1[i1 > i0]) - 1
  fn <- src[i0:i1]

  sub1 <- function(fn, pat, rep, what) {
    h <- grep(pat, fn, fixed = TRUE)
    if (length(h) != 1) stop(sprintf("anchor '%s' matched %d lines (expected 1)", what, length(h)))
    fn[h] <- sub(pat, rep, fn[h], fixed = TRUE); fn
  }

  fn <- sub1(fn,
    "covar.names<-unique(c(names(mf.s12[,-1]),names(mf.s13[,-1]),names(mf.s23[,-1])))",
    "covar.names<-unique(c(names(mf.s12[,-1,drop=FALSE]),names(mf.s13[,-1,drop=FALSE]),names(mf.s23[,-1,drop=FALSE])))",
    "covar.names")

  fn <- sub1(fn, "for(ncov in 1:length(covar.names)){",
                 "for(ncov in seq_along(covar.names)){", "cov.density loop")

  fn <- sub1(fn, "  SurvPred<-array(dim=c(length(tau),nsamp,3))", paste(
    "  .key  <- do.call(paste, c(as.data.frame(cbind(datls$mat12,datls$mat13,datls$mat23)), sep=\"\\r\"))",
    "  .uidx <- which(!duplicated(.key)); .umap <- match(.key, .key[.uidx])",
    "  SurvPred<-array(dim=c(length(tau),nsamp,3))", sep = "\n"), "SurvPred allocation")

  fn <- sub1(fn, "    exb23c<-exp(beta23%*%t((model.matrix(formu23,mfc.s23)[,-1])))", paste(
    "    exb23c<-exp(beta23%*%t((model.matrix(formu23,mfc.s23)[,-1])))",
    "    .ckey  <- do.call(paste, c(as.data.frame(cbind(",
    "                 model.matrix(formu12,mfc.s12)[,-1,drop=FALSE],",
    "                 model.matrix(formu13,mfc.s13)[,-1,drop=FALSE],",
    "                 model.matrix(formu23,mfc.s23)[,-1,drop=FALSE])), sep=\"\\r\"))",
    "    .cuidx <- which(!duplicated(.ckey)); .cumap <- match(.ckey, .ckey[.cuidx])",
    sep = "\n"), "exb23c")

  h <- grep("      P01.0t.j<-numeric()", fn, fixed = TRUE)
  if (length(h) != 2) stop(sprintf("expected 2 P01 loops, found %d", length(h)))
  for (n in seq_along(h)) {
    j <- h[n]
    stopifnot(
      fn[j + 1] == "      for(kk in 1:npat) P01.0t.j[kk]<-integrate(P01Fun,0,ti,k=kk)$value",
      fn[j + 2] == "      return(P01.0t.j)")
    ix <- if (n == 1) ".uidx" else ".cuidx"
    mp <- if (n == 1) ".umap" else ".cumap"
    fn[j]     <- sprintf("      P01.u<-numeric(length(%s))", ix)
    fn[j + 1] <- sprintf("      for(kk in seq_along(%s)) P01.u[kk]<-.p01int(P01Fun,ti,%s[kk])", ix, ix)
    fn[j + 2] <- sprintf("      return(P01.u[%s])", mp)
  }

  .p01int <- function(f, ti, k) {
    v <- tryCatch(integrate(f, 0, ti, k = k, subdivisions = 2000L, rel.tol = 1e-6)$value,
                  error = function(e) NA_real_)
    if (!is.na(v)) return(v)
    n <- 2000L; x <- seq(0, ti, length.out = n + 1L)
    y <- vapply(x, f, numeric(1), k = k)
    y[!is.finite(y)] <- 0
    (ti / n) / 3 * (y[1] + y[n + 1L] +
                    4 * sum(y[seq(2, n, 2)]) + 2 * sum(y[seq(3, n - 1, 2)]))
  }
  assign(".p01int", .p01int, envir = globalenv())

  eval(parse(text = paste(fn, collapse = "\n")), envir = globalenv())
  invisible(TRUE)
}
.patch_logrrcalc_fast()
cat("logRRCalc patched (empty-covariate-safe + unique-pattern collapse).\n")
