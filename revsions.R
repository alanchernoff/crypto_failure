#!/usr/bin/env Rscript
# =====================================================================
# Cryptocurrency Failure Rates: Size, Survival, and Market Shocks
# Replication script for the resubmission
#
# Inputs  (in DATA_DIR):  coin_level3.csv   coin-level census (33,966 coins)
#                         Exhibits.xlsx     original cohort-year aggregates
# Outputs (in OUT_DIR):   tables/*.csv      one csv per table
#                         figures/*.png     figures
#                         results_tables.docx   every table and figure in Word
#
# Usage:  Rscript analysis.R [DATA_DIR] [OUT_DIR]
# Packages: readxl, dplyr, survival, sandwich, lmtest, data.table,
#           officer, flextable
# =====================================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(survival)
  library(sandwich); library(lmtest); library(data.table)
  library(officer); library(flextable)
})

args     <- commandArgs(trailingOnly = TRUE)
DATA_DIR <- "C:/Users/acher/research/JAI_JPM_PMR"
OUT_DIR  <- "C:/Users/acher/research/JAI_JPM_PMR/output"
dir.create(file.path(OUT_DIR, "tables"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)

TIER_LEVELS <- c("Below $100K", "$100K\u2013$1M", "$1M\u2013$10M", "Above $10M")
SHOCK_YEARS <- 2022:2023
save_tab <- function(df, name) {
  write.csv(df, file.path(OUT_DIR, "tables", paste0(name, ".csv")), row.names = FALSE)
  df
}
fmt_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))

# ---------------------------------------------------------------------
# PART 1. Reconcile the original cohort-year aggregates (editor / reviewers)
# ---------------------------------------------------------------------
xl <- read_excel(file.path(DATA_DIR, "Exhibits.xlsx"), sheet = "Sheet 1",
                 col_names = FALSE, .name_repair = "minimal")
agg <- data.frame(
  year        = as.integer(gsub("[^0-9]", "", as.character(unlist(xl[3:16, 1])))),
  nc_total    = as.numeric(unlist(xl[3:16, 2])),
  nc_delisted = as.numeric(unlist(xl[3:16, 3])),
  nc_failed   = as.numeric(unlist(xl[3:16, 5])),
  nc_unres    = as.numeric(unlist(xl[3:16, 7])),
  mc_total    = as.numeric(unlist(xl[3:16, 12])),
  mc_failed   = as.numeric(unlist(xl[3:16, 13]))
)
unres_all   <- as.numeric(xl[[17, 7]])          # 3,813
unres_mc    <- as.numeric(xl[[17, 16]])         # 206 (in the >= $1M panel)
match_rate  <- (sum(agg$nc_total) - sum(agg$nc_unres)) / sum(agg$nc_total)

# sub-$1M tier two ways
agg <- agg %>% mutate(
  below_failed    = nc_failed - mc_failed,
  below_total_old = nc_total - mc_total,               # as originally built
  below_total_fix = (nc_total - nc_unres) - mc_total   # unmatched excluded
)

grouped_or <- function(d, below_total_col) {
  long <- rbind(
    data.frame(year = d$year, tier = 1, failed = d$mc_failed,    total = d$mc_total),
    data.frame(year = d$year, tier = 0, failed = d$below_failed, total = d[[below_total_col]]))
  m  <- glm(cbind(failed, total - failed) ~ tier + factor(year),
            family = binomial, data = long)
  se <- sqrt(diag(vcovHC(m, type = "HC1")))["tier"]
  b  <- coef(m)["tier"]
  c(OR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se),
    p = 2 * pnorm(-abs(b / se)))
}
a1323_old <- grouped_or(filter(agg, year <= 2023), "below_total_old")
a1323_fix <- grouped_or(filter(agg, year <= 2023), "below_total_fix")
a1319_fix <- grouped_or(filter(agg, year <= 2019), "below_total_fix")
agg <- agg %>% mutate(fr_below = below_failed / below_total_fix,
                      fr_above = mc_failed / mc_total)
n_less <- sum(agg$fr_above[agg$year <= 2023] < agg$fr_below[agg$year <= 2023])

# Exhibit 4 on both samples
ex4 <- function(d, label) {
  d <- d %>% mutate(dshare = nc_delisted / nc_failed,
                    shock = as.numeric(year %in% SHOCK_YEARS),
                    yc = year - mean(year))
  m  <- lm(dshare ~ shock + yc, data = d, weights = nc_failed)
  V  <- NeweyWest(m, lag = 1, prewhite = FALSE)
  ct <- coeftest(m, vcov. = V)
  data.frame(Sample = label, N = nrow(d),
             Shock = ct["shock", 1], Shock_SE = ct["shock", 2],
             Trend = ct["yc", 1],    Trend_SE = ct["yc", 2])
}
ex4_tab <- save_tab(rbind(ex4(agg, "2013\u20132026 (as published)"),
                          ex4(filter(agg, year <= 2023), "2013\u20132023 (as stated in Methods)")),
                    "t01b_exhibit4_reproduction")

recon <- save_tab(data.frame(
  Item = c("CoinGecko match rate as stated in original manuscript",
           "CoinGecko match rate implied by workbook (matched / total, unfiltered panel)",
           "Unmatched coins in unfiltered panel",
           "  of which in the >= $1M panel",
           "  of which sub-$1M (unfiltered minus >= $1M)",
           "H1 odds ratio, 2013-2023, unmatched coins left in sub-$1M denominator",
           "H1 odds ratio, 2013-2023, unmatched coins excluded (as Methods state)",
           "H1 odds ratio, 2013-2019 only, unmatched excluded",
           "Cohorts (2013-2023) in which the >= $1M tier fails less, consistent denominator"),
  Value = c("90.18%", sprintf("%.2f%%", 100 * match_rate),
            format(unres_all, big.mark = ","), format(unres_mc, big.mark = ","),
            format(unres_all - unres_mc, big.mark = ","),
            sprintf("%.3f [%.3f, %.3f]", a1323_old["OR.tier"], a1323_old["lo.tier"], a1323_old["hi.tier"]),
            sprintf("%.3f [%.3f, %.3f]", a1323_fix["OR.tier"], a1323_fix["lo.tier"], a1323_fix["hi.tier"]),
            sprintf("%.3f [%.3f, %.3f]", a1319_fix["OR.tier"], a1319_fix["lo.tier"], a1319_fix["hi.tier"]),
            sprintf("%d of 11", n_less))), "t01_reconciliation")

# ---------------------------------------------------------------------
# PART 2. Coin-level data
# ---------------------------------------------------------------------
raw <- read.csv(file.path(DATA_DIR, "coin_level3.csv"), check.names = FALSE,
                stringsAsFactors = FALSE, na.strings = c("", "NA"))
names(raw) <- c("coin_id", "name", "symbol", "listing", "peak_cap", "peak_date",
                "cohort", "status", "cg_match", "match_conf",
                "date_delisted", "date_collapsed", "mcap90")
for (v in c("listing", "peak_date", "date_delisted", "date_collapsed"))
  raw[[v]] <- as.Date(raw[[v]])

STUDY_END <- max(c(raw$listing, raw$date_delisted, raw$date_collapsed), na.rm = TRUE)
cat("Study end (administrative censoring date):", format(STUDY_END), "\n")

d <- raw %>% mutate(
  event_date = pmin(date_collapsed, date_delisted, na.rm = TRUE),
  event      = as.integer(!is.na(event_date)),
  exit_date  = dplyr::if_else(is.na(event_date), STUDY_END, event_date),
  origin     = listing + 90,                          # size is measured at day 90; risk starts here
  at_risk90  = exit_date > (listing + 90),            # survived past the size-measurement horizon
  dur_full   = as.numeric(exit_date - listing) + 1,   # time from listing (sensitivity only)
  dur        = as.numeric(exit_date - (listing + 90)),  # days at risk after day 90
  cause      = dplyr::case_when(
    !is.na(date_collapsed) & (is.na(date_delisted) | date_collapsed <= date_delisted) ~ "collapse",
    !is.na(date_delisted) ~ "delisted_first",
    TRUE ~ "censored"),
  ly         = pmax(as.integer(format(listing, "%Y")), 2013L),   # pool 2009-2012 into 2013
  ly_raw     = as.integer(format(listing, "%Y")),
  tier90     = cut(mcap90, c(-Inf, 1e5, 1e6, 1e7, Inf), right = FALSE, labels = TIER_LEVELS),
  peak_cap0  = ifelse(is.na(peak_cap), 0, peak_cap),
  match_grp  = factor(match_conf, levels = c("exact", "symbol-only", "ambiguous", "unmatched"))
)
d$ev_collapse <- as.integer(d$cause == "collapse")
d$ev_delist   <- as.integer(d$cause == "delisted_first")
early_fail <- d %>% filter(!at_risk90)                 # failed or exited within 90 days of listing
w <- d %>% filter(!is.na(mcap90), at_risk90)            # estimation sample: size observed, alive at day 90

# ---- Sample construction and descriptive statistics -------------------
flow <- save_tab(data.frame(
  Step = c("Coins in coin-level census",
           "Less: 90-day market cap missing",
           "Less: failed or exited within 90 days of listing (size not yet measurable)",
           "Estimation sample (90-day cap observed, at risk at day 90)",
           "  of which failed after day 90 (first of collapse / delisting)",
           "  of which right-censored (no failure event by study end)"),
  N = c(nrow(d), sum(is.na(d$mcap90)), sum(!is.na(d$mcap90) & !d$at_risk90), nrow(w),
        sum(w$event), sum(1 - w$event))), "t02_sample_flow")

desc_tier <- w %>% group_by(Tier = tier90) %>% summarise(
  Coins = n(), Collapse = sum(ev_collapse), Delisted_first = sum(ev_delist),
  Censored = sum(cause == "censored"), Failed_pct = 100 * mean(event),
  Median_days_at_risk = median(dur), .groups = "drop")
desc_tier <- bind_rows(desc_tier %>% mutate(Tier = as.character(Tier)),
                       data.frame(Tier = "All estimation-sample coins", Coins = nrow(w), Collapse = sum(w$ev_collapse),
                                  Delisted_first = sum(w$ev_delist), Censored = sum(w$cause == "censored"),
                                  Failed_pct = 100 * mean(w$event), Median_days_at_risk = median(w$dur)))
desc_tier <- save_tab(desc_tier, "t03_descriptives_by_tier")

desc_cohort <- save_tab(as.data.frame.matrix(table(w$ly_raw, w$tier90)) %>%
                          tibble::rownames_to_column("Listing_year") %>%
                          mutate(Total = rowSums(across(all_of(TIER_LEVELS))),
                                 Failed = as.integer(tapply(w$event, w$ly_raw, sum)[Listing_year]),
                                 Failed_pct = round(100 * Failed / Total, 1)), "t04_coins_by_cohort_and_tier")

# ---- Matched vs. unmatched --------------------------------------------
mm <- d %>% group_by(Match = match_grp) %>% summarise(
  Coins = n(), Failed_pct = 100 * mean(event),
  Collapse_pct = 100 * mean(ev_collapse), Delisted_first_pct = 100 * mean(ev_delist),
  Pct_90d_cap_ge_1M = 100 * mean(mcap90 >= 1e6, na.rm = TRUE),
  Pct_peak_ge_1M = 100 * mean(peak_cap0 >= 1e6),
  Pct_missing_90d = 100 * mean(is.na(mcap90)),
  Median_listing_year = median(ly_raw), .groups = "drop")
d$matched <- as.integer(d$match_conf != "unmatched")
chi_fail <- chisq.test(table(d$matched, d$event))$p.value
mm <- save_tab(mm %>% mutate(Match = as.character(Match)), "t05_matched_vs_unmatched")
cat("Chi-square p (failure vs matched):", chi_fail, "\n")

# ---- Helpers for Cox and discrete-time hazard --------------------------
tidy_cox <- function(m, pattern = "^tier90") {
  s <- summary(m); ci <- s$conf.int; co <- s$coefficients
  idx <- grep(pattern, rownames(co))
  data.frame(Term = sub(pattern, "", rownames(co)[idx]),
             HR = ci[idx, "exp(coef)"], Lo = ci[idx, "lower .95"], Hi = ci[idx, "upper .95"],
             p = co[idx, "Pr(>|z|)"], row.names = NULL)
}
tidy_glm <- function(m, cluster, keep) {
  V  <- vcovCL(m, cluster = cluster, type = "HC1")
  co <- coeftest(m, vcov. = V)
  idx <- grep(keep, rownames(co))
  b <- co[idx, 1]; se <- co[idx, 2]
  data.frame(Term = rownames(co)[idx], OR = exp(b), Lo = exp(b - 1.96 * se),
             Hi = exp(b + 1.96 * se), p = co[idx, 4], row.names = NULL)
}

# ---- Cox models: any failure and cause-specific ---------------------------
w$tier90 <- relevel(factor(w$tier90, levels = TIER_LEVELS), ref = "Below $100K")
cox_any <- coxph(Surv(dur, event)       ~ tier90 + factor(ly), data = w)
cox_col <- coxph(Surv(dur, ev_collapse) ~ tier90 + factor(ly), data = w)
cox_del <- coxph(Surv(dur, ev_delist)   ~ tier90 + factor(ly), data = w)
cox_tab <- save_tab(rbind(
  cbind(Model = "Any failure",                      tidy_cox(cox_any)),
  cbind(Model = "Price collapse (delisting censored)", tidy_cox(cox_col)),
  cbind(Model = "Delisted first (collapse censored)",  tidy_cox(cox_del))), "t06_cox_models")
zph <- cox.zph(cox_any)
zph_tab <- save_tab(data.frame(Term = rownames(zph$table), chisq = zph$table[, "chisq"],
                               p = zph$table[, "p"], row.names = NULL), "t06b_ph_test")

# ---- Discrete-time (coin-year) hazard panel -----------------------------
build_panel <- function(dat) {
  dt <- as.data.table(dat)
  dt[, `:=`(L = as.integer(format(origin, "%Y")),
            E = as.integer(format(exit_date, "%Y")), cid = .I)]
  dt[, n_y := E - L + 1L]
  p <- dt[rep(seq_len(.N), n_y)]
  p[, year := L + sequence(dt$n_y) - 1L]
  p[, last := year == E]
  p[, `:=`(y_any = as.integer(last & event == 1L),
           y_col = as.integer(last & cause == "collapse"),
           y_del = as.integer(last & cause == "delisted_first"),
           age   = year - L, shock = as.integer(year %in% SHOCK_YEARS))]
  p[, age_c := factor(pmin(age, 5))]
  p[, Lg := factor(pmax(as.integer(format(listing, "%Y")), 2013L))]
  as.data.frame(p)
}
pn <- build_panel(w)
cat("Coin-year panel:", nrow(pn), "rows,", length(unique(pn$cid)), "coins\n")

dth <- function(y, rhs, data = pn, keep = "tier90|shock") {
  m <- glm(as.formula(paste(y, "~", rhs)), family = binomial, data = data)
  tidy_glm(m, data$cid, keep)
}
base <- "tier90 + shock + age_c"
dth_main <- save_tab(rbind(
  cbind(Outcome = "Any failure",      dth("y_any", base)),
  cbind(Outcome = "Price collapse",   dth("y_col", base)),
  cbind(Outcome = "Delisted first",   dth("y_del", base))), "t07_dth_main")

dth_cohortfe <- save_tab(rbind(
  cbind(Outcome = "Any failure",    dth("y_any", paste(base, "+ Lg"), keep = "^shock$")),
  cbind(Outcome = "Price collapse", dth("y_col", paste(base, "+ Lg"), keep = "^shock$")),
  cbind(Outcome = "Delisted first", dth("y_del", paste(base, "+ Lg"), keep = "^shock$"))),
  "t07b_shock_with_cohort_fe")

# tier effect in subsamples (answers "size effect is confined to the shock window")
sub_fit <- function(lbl, dat) cbind(Sample = lbl, dth("y_any", "tier90 + age_c", dat, "^tier90"))
dth_sub <- save_tab(rbind(
  sub_fit("All years", pn),
  sub_fit("Shock years (2022-23) removed", pn[pn$shock == 0, ]),
  sub_fit("Shock years only", pn[pn$shock == 1, ]),
  sub_fit("Pre-2020 only", pn[pn$year <= 2019, ])), "t08_tier_effect_subsamples")

dth_int <- save_tab(cbind(Outcome = "Any failure",
                          dth("y_any", "tier90 * shock + age_c", keep = "tier90|shock")),
                    "t09_tier_by_shock_interaction")

# ---- Placebo-window (permutation) inference for the shock effect -------------
q  <- pn[pn$year >= 2015 & pn$year <= 2025, ]
X0 <- model.matrix(~ tier90 + age_c, data = q)
pairs <- combn(2015:2025, 2)
run_pl <- function(y) apply(pairs, 2, function(pr) {
  X <- cbind(X0, s = as.numeric(q$year %in% pr))
  coef(suppressWarnings(glm.fit(X, q[[y]], family = binomial())))["s"]
})
pl <- data.frame(y1 = pairs[1, ], y2 = pairs[2, ],
                 delist = run_pl("y_del"), collapse = run_pl("y_col"), any = run_pl("y_any"))
actual_idx <- which(pl$y1 == 2022 & pl$y2 == 2023)
placebo_tab <- save_tab(data.frame(
  Outcome = c("Delisted first", "Price collapse", "Any failure"),
  Actual_OR_2022_23 = exp(c(pl$delist[actual_idx], pl$collapse[actual_idx], pl$any[actual_idx])),
  Placebo_windows = nrow(pl),
  Share_placebo_ge_actual = c(mean(pl$delist >= pl$delist[actual_idx]),
                              mean(pl$collapse >= pl$collapse[actual_idx]),
                              mean(pl$any >= pl$any[actual_idx])),
  Share_placebo_le_actual = c(mean(pl$delist <= pl$delist[actual_idx]),
                              mean(pl$collapse <= pl$collapse[actual_idx]),
                              mean(pl$any <= pl$any[actual_idx]))), "t10_placebo_windows")
invisible(save_tab(pl, "t10b_placebo_all_windows"))

# ---- Threshold sensitivity (binary cutoffs on 90-day cap) -------------------
thr <- do.call(rbind, lapply(c(1e5, 1e6, 1e7), function(cut) {
  w$hi <- as.integer(w$mcap90 >= cut)
  f <- function(y) tidy_cox(coxph(as.formula(paste0("Surv(dur,", y, ") ~ hi + factor(ly)")), data = w), "^hi")
  data.frame(Cutoff = paste0("$", format(cut, big.mark = ",", scientific = FALSE)),
             Coins_above = sum(w$hi), Pct_above = 100 * mean(w$hi),
             HR_any = f("event")$HR, Lo_any = f("event")$Lo, Hi_any = f("event")$Hi,
             HR_collapse = f("ev_collapse")$HR, HR_delisted_first = f("ev_delist")$HR)
}))
thr_tab <- save_tab(thr, "t11_threshold_sensitivity")

# ---- Mechanical-link diagnostics ------------------------------------------------
d$hi_peak <- as.integer(d$peak_cap0 >= 1e6)
d$hi_90   <- as.integer(!is.na(d$mcap90) & d$mcap90 >= 1e6)
w$hi_peak <- as.integer(w$peak_cap0 >= 1e6); w$hi_90 <- as.integer(w$mcap90 >= 1e6)
cx <- function(lbl, f, dat, term) cbind(Measure = lbl, N = nrow(dat), tidy_cox(coxph(f, data = dat), term))
mech <- rbind(
  cx("A. Peak cap >= $1M, hazard from listing, full census (original definition)",
     Surv(dur_full, event) ~ hi_peak + factor(ly), d, "^hi_peak"),
  cx("B. Peak cap >= $1M, hazard from day 90, estimation sample",
     Surv(dur, event) ~ hi_peak + factor(ly), w, "^hi_peak"),
  cx("C. 90-day cap >= $1M, hazard from day 90, estimation sample (preferred)",
     Surv(dur, event) ~ hi_90 + factor(ly), w, "^hi_90"),
  cx("D. 90-day cap >= $1M, hazard from listing (failures in first 90 days retained)",
     Surv(dur_full, event) ~ hi_90 + factor(ly), filter(d, !is.na(mcap90)), "^hi_90"))
mech_tab <- save_tab(mech, "t12_mechanical_link_diagnostic")

# ---- Tier effect by age band (non-proportional hazards check) ------------------
age_fit <- function(lbl, keepage) cbind(Age_since_day90 = lbl,
                                        dth("y_any", "tier90 + factor(Lg)", pn[pn$age %in% keepage, ], "^tier90"))
age_tab <- save_tab(rbind(age_fit("Year 0", 0), age_fit("Year 1", 1), age_fit("Year 2", 2),
                          age_fit("Years 3+", 3:20)), "t08b_tier_effect_by_age")

# ---- Robustness: match-confidence subsamples and missing 90-day cap ---------------
sub_cox <- function(lbl, dat) cbind(Sample = lbl, N = nrow(dat), tidy_cox(
  coxph(Surv(dur, event) ~ tier90 + factor(ly), data = dat)))
rb <- rbind(
  sub_cox("All coins (main)", w),
  sub_cox("Exact CoinGecko match only", filter(w, match_conf == "exact")),
  sub_cox("Exact + symbol-only + ambiguous (any match)", filter(w, match_conf != "unmatched")),
  sub_cox("Unmatched only", filter(w, match_conf == "unmatched")),
  sub_cox("Positive 90-day cap only (zero-cap coins excluded)", filter(w, mcap90 > 0)))
d2 <- d %>% filter(at_risk90) %>% mutate(tier90 = as.character(tier90)); d2$tier90[is.na(d2$mcap90)] <- "Below $100K"
d2$tier90 <- factor(d2$tier90, levels = TIER_LEVELS)
rb <- rbind(rb, sub_cox("Missing 90-day cap treated as below $100K", d2))
rb_ctrl <- coxph(Surv(dur, event) ~ tier90 + match_grp + factor(ly), data = w)
rb <- rbind(rb, cbind(Sample = "Main + match-confidence controls", N = nrow(w), tidy_cox(rb_ctrl)))
rb_tab <- save_tab(rb, "t13_robustness_samples")

# ---------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------
fig_km <- file.path(OUT_DIR, "figures", "fig1_km_by_tier.png")
png(fig_km, width = 1800, height = 1100, res = 220)
par(mar = c(4.6, 4.6, 0.8, 0.8))
km <- survfit(Surv(dur / 365.25, event) ~ tier90, data = w)
cols <- c("#B23A48", "#E08E45", "#3D7EA6", "#1B3A57")
plot(km, col = cols, lwd = 2.2, conf.int = FALSE, xlim = c(0, 10), ylim = c(0, 1),
     xlab = "Years since day 90 after listing (start of risk period)", ylab = "Share of coins with no failure event",
     main = "", bty = "l", las = 1)
legend("topright", TIER_LEVELS, col = cols, lwd = 2.2, bty = "n", title = "90-day market cap")
dev.off()

fig_pl <- file.path(OUT_DIR, "figures", "fig2_placebo_windows.png")
png(fig_pl, width = 1800, height = 1100, res = 220)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.2, 1))
hist(exp(pl$delist), breaks = 14, col = "grey80", border = "white", main = "Delisted-first hazard",
     xlab = "Odds ratio for placebo two-year window", las = 1)
abline(v = exp(pl$delist[actual_idx]), col = "#B23A48", lwd = 2.5)
hist(exp(pl$any), breaks = 14, col = "grey80", border = "white", main = "Any failure",
     xlab = "Odds ratio for placebo two-year window", las = 1)
abline(v = exp(pl$any[actual_idx]), col = "#B23A48", lwd = 2.5)
dev.off()

# ---------------------------------------------------------------------
# Word document with every table and figure
# ---------------------------------------------------------------------
mk_ft <- function(df, digits = 3) {
  df <- as.data.frame(df)
  for (j in seq_along(df)) if (is.numeric(df[[j]])) {
    x <- df[[j]]; nm <- names(df)[j]
    df[[j]] <- if (nm %in% c("p", "Share_placebo_ge_actual", "Share_placebo_le_actual")) {
      ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))
    } else if (all(is.na(x) | abs(x - round(x)) < 1e-9)) {
      format(round(x), big.mark = ",", trim = TRUE)
    } else sprintf(paste0("%.", digits, "f"), x)
  }
  ft <- flextable(df); ft <- fontsize(ft, size = 8.5, part = "all")
  ft <- font(ft, fontname = "Times New Roman", part = "all"); ft <- autofit(ft)
  ft <- theme_booktabs(ft); bold(ft, part = "header")
}
doc <- read_docx()
bold_par <- function(doc, txt, size = 11) body_add_fpar(doc, fpar(ftext(txt, fp_text(bold = TRUE, font.size = size, font.family = "Times New Roman"))))
add_h <- function(doc, txt) bold_par(doc, txt, 13)
add_t <- function(doc, df, title, note = NULL, digits = 3) {
  doc <- bold_par(doc, title, 10.5)
  doc <- body_add_flextable(doc, mk_ft(df, digits))
  if (!is.null(note)) doc <- body_add_par(doc, note, style = "Normal")
  body_add_par(doc, "", style = "Normal")
}
doc <- bold_par(doc, "Replication output: Cryptocurrency Failure Rates (resubmission)", 15)
doc <- body_add_par(doc, paste0("Generated by analysis.R. Administrative censoring date: ", format(STUDY_END),
                                ". Coin-level census: ", format(nrow(d), big.mark = ","), " coins."), style = "Normal")
doc <- add_h(doc, "A. Reconciliation of the original aggregate exhibits")
doc <- add_t(doc, recon, "Table A1. Checks on the original submission")
doc <- add_t(doc, ex4_tab, "Table A2. Original Exhibit 4 on the published and the stated sample",
             "Weighted least squares, weights = annual failure count, Newey-West (lag 1) standard errors.", 4)
doc <- add_h(doc, "B. Sample and descriptive statistics")
doc <- add_t(doc, flow, "Table B1. Sample construction")
doc <- add_t(doc, desc_tier, "Table B2. Coins and failures by 90-day market-cap tier", NULL, 1)
doc <- add_t(doc, desc_cohort, "Table B3. Coins by listing year and 90-day market-cap tier", NULL, 1)
doc <- add_t(doc, mm, "Table B4. Matched versus unmatched coins (full census)",
             paste0("Chi-square test of failure share, matched vs unmatched: p = ", fmt_p(chi_fail), "."), 1)
doc <- add_h(doc, "C. Survival models")
doc <- add_t(doc, cox_tab, "Table C1. Cox proportional hazards models (hazard ratios vs. below $100K)",
             "Clock starts 90 days after listing; coins that failed within 90 days are excluded. Listing-year fixed effects.")
doc <- add_t(doc, zph_tab, "Table C2. Proportional-hazards test for the any-failure Cox model")
doc <- bold_par(doc, "Figure C1. Kaplan-Meier curves by 90-day market-cap tier", 10.5)
doc <- body_add_img(doc, fig_km, width = 6, height = 3.7)
doc <- add_t(doc, dth_main, "Table C3. Discrete-time hazard: tier and shock effects (odds ratios)",
             "Coin-year panel; age fixed effects (0-5+ years); standard errors clustered by coin.")
doc <- add_t(doc, dth_sub, "Table C4. Tier effect in subsamples (any failure)")
doc <- add_t(doc, age_tab, "Table C4b. Tier effect by years since day 90 (any failure; addresses non-proportional hazards)")
doc <- add_t(doc, dth_int, "Table C5. Tier-by-shock interaction (any failure)")
doc <- add_t(doc, dth_cohortfe, "Table C6. Shock effect with listing-cohort fixed effects")
doc <- add_h(doc, "D. Placebo-window inference for the 2022-2023 shock")
doc <- add_t(doc, placebo_tab, "Table D1. Actual 2022-23 effect against all 55 placebo two-year windows (2015-2025)")
doc <- bold_par(doc, "Figure D1. Distribution of placebo-window odds ratios", 10.5)
doc <- body_add_img(doc, fig_pl, width = 6, height = 3.7)
doc <- add_h(doc, "E. Robustness")
doc <- add_t(doc, thr_tab, "Table E1. Alternative size cutoffs (binary indicator, Cox hazard ratios)")
doc <- add_t(doc, mech_tab, "Table E2. Peak-based versus fixed-horizon size, and the role of the day-90 risk origin (hazard ratios for >= $1M)")
doc <- add_t(doc, rb_tab, "Table E3. Tier effect by match confidence and treatment of missing 90-day cap")
print(doc, target = file.path(OUT_DIR, "results_tables.docx"))
cat("Done. Outputs in", OUT_DIR, "\n")
