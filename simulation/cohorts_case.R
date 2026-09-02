DAYS_PER_MO <- 30.44

.e <- new.env()
utils::data(crc_adsl, package = "MultiStateModels", envir = .e)
utils::data(crc_path, package = "MultiStateModels", envir = .e)
CRC_ADSL <- .e$crc_adsl
CRC_PATH <- .e$crc_path

.prep_adsl <- function(d) data.frame(
  trt      = as.numeric(as.character(d$ARM)),
  resp     = as.numeric(as.character(d$Response)),
  TTP      = d$PFS / DAYS_PER_MO,
  TTPevent = ifelse(d$PFSevent == 1 & d$PFS == d$OS & d$OSevent == 1, 0, d$PFSevent),
  OS       = d$OS / DAYS_PER_MO,
  OSevent  = d$OSevent,
  t_resp   = NA_real_)

.prep_path <- function() {
  p  <- CRC_PATH
  id <- unique(p$SUBJID)
  b  <- p[match(id, p$SUBJID), c("SUBJID", "ARM", "PFS", "PFSevent", "OS", "OSevent")]
  isr <- as.character(p$Resp) == "1"
  tr  <- tapply(p$DAY[isr], factor(p$SUBJID[isr], levels = id), min)
  data.frame(
    trt      = as.numeric(as.character(b$ARM)),
    resp     = as.integer(!is.na(tr[b$SUBJID])),
    TTP      = b$PFS / DAYS_PER_MO,
    TTPevent = ifelse(b$PFSevent == 1 & b$PFS == b$OS & b$OSevent == 1, 0, b$PFSevent),
    OS       = b$OS / DAYS_PER_MO,
    OSevent  = b$OSevent,
    t_resp   = as.numeric(tr[b$SUBJID]) / DAYS_PER_MO)
}

COHORTS_CASE <- list(
  `309 all` = list(
    t_star = 12, label = "309 panit+FOLFOX (all)", src = "adsl", resp_def = "final",
    filt = function(a) a[a$StudyID == "309", ]),
  `309 ECOG 0` = list(
    t_star = 12, label = "309 panit+FOLFOX (ECOG 0)", src = "adsl", resp_def = "final",
    filt = function(a) a[a$StudyID == "309" & a$B_ECOG == "0", ]),
  `263 all` = list(
    t_star = 12, label = "263 panit+FOLFIRI (all)", src = "adsl", resp_def = "final",
    filt = function(a) a[a$StudyID == "263", ]),
  `263 path` = list(
    t_star = 12, label = "263 subset with dated RECIST (n=583)", src = "path",
    resp_def = "observed_by_t", filt = NULL),
  `264+309 all` = list(
    t_star = 12, label = "264+309 panit+FOLFOX pooled (all)", src = "adsl", resp_def = "final",
    filt = function(a) a[a$StudyID %in% c("264", "309"), ]),
  `262 all` = list(
    t_star = 12, label = "262 panit+bev+chemo (all)", src = "adsl", resp_def = "final",
    filt = function(a) a[a$StudyID == "262", ]),
  `309 KRAS wt` = list(
    t_star = 12, label = "309 panit+FOLFOX (KRAS wild-type)", src = "adsl", resp_def = "final",
    filt = function(a) a[a$StudyID == "309" & !is.na(a$KRAS) & a$KRAS == "wild_type", ]),
  `310 all` = list(
    t_star = 6, label = "310 panit+BSC (all)", src = "adsl", resp_def = "final",
    filt = function(a) a[a$StudyID == "310", ]),
  `310 KRAS wt` = list(
    t_star = 6, label = "310 panit+BSC (KRAS wild-type)", src = "adsl", resp_def = "final",
    filt = function(a) a[a$StudyID == "310" & !is.na(a$KRAS) & a$KRAS == "wild_type", ]),
  `263 age>65` = list(
    t_star = 12, label = "263 panit+FOLFIRI (age>65)", src = "adsl", resp_def = "final",
    filt = function(a) a[a$StudyID == "263" & a$AGEG == ">65", ])
)

case_data <- function(key) {
  s <- COHORTS_CASE[[key]]
  if (is.null(s)) stop("unknown cohort: ", key)
  d <- if (s$src == "path") .prep_path() else .prep_adsl(s$filt(CRC_ADSL))
  d[!is.na(d$trt) & !is.na(d$resp) & !is.na(d$OS) & !is.na(d$TTP), ]
}

case_cut <- function(d, t, resp_def) {
  if (is.infinite(t)) return(d)
  o <- d
  o$TTPevent <- ifelse(o$TTP <= t, o$TTPevent, 0); o$TTP <- pmin(o$TTP, t)
  o$OSevent  <- ifelse(o$OS  <= t, o$OSevent,  0); o$OS  <- pmin(o$OS,  t)
  if (resp_def == "observed_by_t") o$resp <- as.integer(!is.na(o$t_resp) & o$t_resp <= t)
  o
}

CASE_CUTS <- c(6, 9, 12, 18, 24, Inf)
cut_lab <- function(t) ifelse(is.infinite(t), "final", sprintf("%gmo", t))
