## =============================================================================
## Crypto Failure Rate Study: Hypothesis Testing Script
##
## H1: Failure risk declines with peak market-cap size ("size/credibility effect")
##     - Tests whether coins that ever crossed the $1M peak-cap threshold have
##       systematically lower failure rates than coins that never did,
##       controlling for cohort year.
##
## H2: The COMPOSITION of failure (delisting vs. price-collapse) shifts in
##     shock years, independent of cohort age.
##     - Tests whether the share of failures attributable to delisting spikes
##       in known market-shock years (2022 Terra/FTX collapse and 2023
##       regulatory fallout), controlling for a linear year trend.
##
## Data inputs (place in your working directory / research folder):
##   milcap.csv  -- coins that ever peaked >= $1M market cap
##                  columns: Year, Total Coins, Failed Total, Active,
##                           Failure Rate (%), Unresolved Matches
##   nocap.csv   -- unfiltered universe of all tracked coins
##                  columns: Year, Total Coins, Delisted/Untracked,
##                           Price-Collapsed, Failed Total, Active,
##                           Unresolved Matches, Failure Rate (%)
##
## Output: crypto_failure_hypothesis_tests.docx (tables + charts + captions)
## =============================================================================

## ---- 0. Setup ---------------------------------------------------------------

# setwd("~/path/to/your/research_folder")   # <-- set this to your data folder

packages <- c("readr", "dplyr", "tidyr", "ggplot2", "sandwich", "lmtest",
              "officer", "flextable", "scales", "broom", "tibble")
to_install <- packages[!packages %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(packages, library, character.only = TRUE))

SOURCE_NOTE <- "Source: CoinMarketCap historical snapshots (2013\u20132026) and CoinGecko API cross-check; author's own calculations."

## ---- 1. Load and clean data --------------------------------------------------

clean_year <- function(y) {
  # "2026 (ytd)" -> 2026, flagged separately as an incomplete/censored cohort
  y_chr <- trimws(as.character(y))
  is_ytd <- grepl("ytd", y_chr, ignore.case = TRUE)
  y_num <- as.numeric(gsub("[^0-9]", "", y_chr))
  list(year = y_num, ytd = is_ytd)
}

milcap_raw <- read_csv("milcap.csv", show_col_types = FALSE)
nocap_raw  <- read_csv("nocap.csv",  show_col_types = FALSE)

names(milcap_raw) <- c("Year", "TotalCoins", "FailedTotal", "Active",
                        "FailureRatePct", "UnresolvedMatches")
names(nocap_raw)  <- c("Year", "TotalCoins", "Delisted", "PriceCollapsed",
                        "FailedTotal", "Active", "UnresolvedMatches",
                        "FailureRatePct")

# Drop the "All years" summary row before parsing years; keep it aside if needed
milcap <- milcap_raw %>% filter(!grepl("all years", Year, ignore.case = TRUE))
nocap  <- nocap_raw  %>% filter(!grepl("all years", Year, ignore.case = TRUE))

parse_block <- function(df) {
  yr_info <- clean_year(df$Year)
  df$YearNum <- yr_info$year
  df$IsYTD   <- yr_info$ytd
  df
}

milcap <- parse_block(milcap)
nocap  <- parse_block(nocap)

## ---- 2. Derive the sub-$1M tier by differencing the two datasets ------------
## milcap = coins that ever peaked >= $1M. nocap = all tracked coins.
## nocap - milcap, matched by year, approximates coins that never crossed $1M peak.

tier_data <- nocap %>%
  select(YearNum, IsYTD, TotalCoins_all = TotalCoins,
         FailedTotal_all = FailedTotal, Active_all = Active) %>%
  inner_join(
    milcap %>% select(YearNum, TotalCoins_mil = TotalCoins,
                       FailedTotal_mil = FailedTotal, Active_mil = Active),
    by = "YearNum"
  ) %>%
  mutate(
    TotalCoins_sub1m = TotalCoins_all - TotalCoins_mil,
    FailedTotal_sub1m = FailedTotal_all - FailedTotal_mil,
    Active_sub1m = Active_all - Active_mil
  )

tier_long <- bind_rows(
  tier_data %>% transmute(YearNum, IsYTD, Tier = "At/above $1M peak",
                           Total = TotalCoins_mil, Failed = FailedTotal_mil,
                           Active = Active_mil),
  tier_data %>% transmute(YearNum, IsYTD, Tier = "Below $1M peak",
                           Total = TotalCoins_sub1m, Failed = FailedTotal_sub1m,
                           Active = Active_sub1m)
) %>%
  mutate(FailureRate = Failed / Total,
         Tier = factor(Tier, levels = c("Below $1M peak", "At/above $1M peak")))

## Main analysis sample excludes immature/right-censored cohorts (2024-2026 ytd),
## since young coins have not had time to fail yet and will mechanically show
## lower failure rates. These are retained separately for a robustness check.
tier_mature   <- tier_long %>% filter(YearNum < 2024)
tier_immature <- tier_long %>% filter(YearNum >= 2024)

## ---- 3. H1: logistic regression of failure on market-cap tier ---------------
## Weighted (grouped) logistic regression: cbind(Failed, Active) ~ Tier + year FE
## Robust (HC) standard errors used given known heteroskedasticity across
## cohorts of very different sizes.

h1_model <- glm(cbind(Failed, Active) ~ Tier + factor(YearNum),
                 family = binomial(link = "logit"), data = tier_mature)

h1_robust <- coeftest(h1_model, vcov. = vcovHC(h1_model, type = "HC1"))

h1_tier_row <- h1_robust["TierAt/above $1M peak", ]
h1_or <- exp(h1_tier_row["Estimate"])
h1_or_ci <- exp(h1_tier_row["Estimate"] + c(-1, 1) * 1.96 * h1_tier_row["Std. Error"])

h1_table_df <- data.frame(
  Term = c("At/above $1M peak (vs. below $1M)"),
  `Odds Ratio` = round(as.numeric(h1_or), 3),
  `95% CI Lower` = round(as.numeric(h1_or_ci[1]), 3),
  `95% CI Upper` = round(as.numeric(h1_or_ci[2]), 3),
  `p-value` = signif(as.numeric(h1_tier_row["Pr(>|z|)"]), 3),
  check.names = FALSE
)

## ---- 4. H1 chart: failure rate by year, split by tier ------------------------

h1_chart <- ggplot(tier_long, aes(x = YearNum, y = FailureRate, color = Tier,
                                   linetype = IsYTD)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "dashed"),
                         guide = "none") +
  labs(title = "Annual Failure Rate by Peak Market-Cap Tier",
       x = "Year", y = "Failure Rate",
       color = "Market-Cap Tier") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave("h1_chart.png", h1_chart, width = 7, height = 4.2, dpi = 300)

## ---- 5. H2: composition of failure (delisting share) vs. shock years --------

nocap_h2 <- nocap %>%
  mutate(DelistShare = Delisted / FailedTotal,
         ShockYear = as.integer(YearNum %in% c(2022, 2023)),
         YearTrend = YearNum - min(YearNum))

## Weighted OLS (weights = FailedTotal, i.e. more failures = more precisely
## estimated share that year). Newey-West HAC SEs used given annual time-series
## structure and potential serial correlation.
h2_model <- lm(DelistShare ~ ShockYear + YearTrend, data = nocap_h2,
                weights = FailedTotal)

h2_robust <- coeftest(h2_model, vcov. = NeweyWest(h2_model, lag = 1, prewhite = FALSE))

h2_mat <- unclass(h2_robust)
h2_table_df <- data.frame(
  Term = c("Intercept", "Shock year (2022\u201323)", "Year trend"),
  Estimate = signif(h2_mat[, "Estimate"], 4),
  `Std. Error` = signif(h2_mat[, "Std. Error"], 4),
  `t-value` = signif(h2_mat[, "t value"], 4),
  `p-value` = signif(h2_mat[, "Pr(>|t|)"], 4),
  check.names = FALSE,
  row.names = NULL
)

## ---- 6. H2 chart: delisting share by year, shock years highlighted -----------

shock_years <- c(2022, 2023)

h2_chart <- ggplot(nocap_h2, aes(x = YearNum, y = DelistShare)) +
  annotate("rect", xmin = min(shock_years) - 0.5, xmax = max(shock_years) + 0.5,
           ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "red") +
  geom_line(linewidth = 1, color = "steelblue") +
  geom_point(size = 2, color = "steelblue") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "Delisting Share of Total Failures by Year",
       subtitle = "Shaded band = shock years (2022 Terra/FTX collapse and 2023 regulatory fallout)",
       x = "Year", y = "Delisted / Total Failures") +
  theme_minimal(base_size = 11)

ggsave("h2_chart.png", h2_chart, width = 7, height = 4.2, dpi = 300)

## ---- 7. Assemble Word document ------------------------------------------------

doc <- read_docx()

doc <- doc %>%
  body_add_par("Cryptocurrency Failure Rate: Hypothesis Test Results", style = "heading 1") %>%
  body_add_par(" ")

## --- H1 section ---
doc <- doc %>%
  body_add_par("Hypothesis 1: Peak market-cap size and survival", style = "heading 2") %>%
  body_add_par(paste("H1 posits that coins which ever crossed a $1M peak market",
                      "cap fail at a lower rate than coins that never did, controlling",
                      "for cohort-year fixed effects. Immature cohorts (2024\u20132026 ytd)",
                      "are excluded from this test due to right-censoring."))

doc <- doc %>%
  body_add_par(" ") %>%
  body_add_img("h1_chart.png", width = 6.2, height = 3.7) %>%
  body_add_par(paste("Figure 1. Annual failure rate by peak market-cap tier, 2013\u20132026.",
                      "Dashed segments denote incomplete (right-censored) cohorts.",
                      SOURCE_NOTE),
               style = "Caption")

h1_flex <- flextable(h1_table_df) %>% autofit()
doc <- doc %>%
  body_add_par(" ") %>%
  body_add_flextable(h1_flex) %>%
  body_add_par(paste("Table 1. Weighted logistic regression of coin failure on peak",
                      "market-cap tier, controlling for cohort-year fixed effects.",
                      "Robust (HC1) standard errors.", SOURCE_NOTE),
               style = "Caption")

## --- H2 section ---
doc <- doc %>%
  body_add_par(" ") %>%
  body_add_par("Hypothesis 2: Shock years and the composition of failure", style = "heading 2") %>%
  body_add_par(paste("H2 posits that the share of failures attributable to delisting",
                      "(as opposed to gradual price collapse) rises disproportionately",
                      "in market-shock years, independent of cohort age."))

doc <- doc %>%
  body_add_par(" ") %>%
  body_add_img("h2_chart.png", width = 6.2, height = 3.7) %>%
  body_add_par(paste("Figure 2. Delisting share of total failures by year, with 2022\u201323",
                      "shock years highlighted.", SOURCE_NOTE),
               style = "Caption")

h2_flex <- flextable(h2_table_df) %>% autofit()
doc <- doc %>%
  body_add_par(" ") %>%
  body_add_flextable(h2_flex) %>%
  body_add_par(paste("Table 2. Weighted OLS regression of annual delisting share on a",
                      "shock-year dummy (2022\u201323) and a linear year trend, weighted by",
                      "total failures per year. Newey-West HAC standard errors.",
                      SOURCE_NOTE),
               style = "Caption")

print(doc, target = "crypto_failure_hypothesis_tests.docx")

cat("Done. Output written to crypto_failure_hypothesis_tests.docx\n")
cat("(along with h1_chart.png and h2_chart.png in the working directory)\n")
