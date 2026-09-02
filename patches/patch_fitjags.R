.patch_fitjags <- function(pkg = MSM_DIR) {
  src <- readLines(file.path(pkg, "R/Mstate_JagsFuncs.R"))
  i0 <- grep("^FitJagsMstate<-function", src)[1]
  i1 <- grep("^[A-Za-z_.]+[[:space:]]*<-[[:space:]]*function", src)
  i1 <- min(i1[i1 > i0]) - 1
  fn <- src[i0:i1]

  j_start <- grep("## initial values", fn, fixed = TRUE)[1]
  j_end   <- grep("inits<-list()", fn, fixed = TRUE)[1] - 1
  stopifnot(length(j_start) == 1, j_end >= j_start)

  new_init <- c(
    "  p12<-ncol(mat12); p13<-ncol(mat13); p23<-ncol(mat23); pp<-p12+p13+p23",
    "  nchain_<-mcmc.par$nchain",
    "  beta12.init<-matrix(rnorm(nchain_*p12,0,0.1),nchain_,p12)",
    "  beta13.init<-matrix(rnorm(nchain_*p13,0,0.1),nchain_,p13)",
    "  beta23.init<-matrix(rnorm(nchain_*p23,0,0.1),nchain_,p23)",
    "  nv12.init<-runif(nchain_,0.8,1.4); nv13.init<-runif(nchain_,0.8,1.4); nv23.init<-runif(nchain_,0.8,1.4)",
    "  .rate0<-1/max(mean(c(T1[T1>0],T2[T2>0]),na.rm=TRUE),1e-3)",
    "  loglamb12.init<-rnorm(nchain_,log(.rate0),0.2)",
    "  loglamb13.init<-rnorm(nchain_,log(.rate0),0.2)",
    "  loglamb23.init<-rnorm(nchain_,log(.rate0),0.2)",
    "  freq.est<-NULL")

  fn2 <- c(fn[seq_len(j_start - 1)], new_init, fn[(j_end + 1):length(fn)])

  i_diag <- grep("^\\s*test<-list\\(gelman\\.coef=", fn2)
  if (length(i_diag) == 1) {
    j_end2 <- grep("geweke\\.bsl=geweke\\.diag", fn2)
    j_end2 <- j_end2[j_end2 >= i_diag][1]
    if (!is.na(j_end2)) {
      fn2[i_diag] <- sub("test<-list\\(", "test<-tryCatch(list(", fn2[i_diag])
      fn2[j_end2] <- paste0(fn2[j_end2],
        ", error=function(e) list(diag_failed=conditionMessage(e)))")
    }
  }
  eval(parse(text = paste(fn2, collapse = "\n")), envir = globalenv())
  invisible(TRUE)
}
.patch_fitjags()
cat("FitJagsMstate patched (robust simple initialization).\n")
