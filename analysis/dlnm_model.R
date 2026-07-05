# dlnm_model.R
# -----------------------------------------------------------------------------
# Publication-grade distributed-lag non-linear model (DLNM) for the
# climate-dengue association in East Java, Indonesia, using the canonical
# `dlnm` R package (Gasparrini 2011). This mirrors the modelling approach in
# Lowe et al. 2021 (Lancet Planetary Health).
#
# Install once:
#   install.packages(c("dlnm","splines","MASS","ggplot2","readr","dplyr"))
#
# Run:
#   Rscript analysis/dlnm_model.R
# -----------------------------------------------------------------------------

library(dlnm)
library(splines)
library(MASS)
library(readr)
library(dplyr)

df <- read_csv("data/processed/analysis_eastjava.csv") |>
  arrange(date) |>
  filter(!is.na(cases), !is.na(temp_mean), !is.na(precip_total), !is.na(population),
         !is.na(oni), !is.na(dmi))

maxlag <- 6

# ---- Cross-basis for temperature -------------------------------------------
cb.temp <- crossbasis(
  df$temp_mean, lag = maxlag,
  argvar = list(fun = "ns", df = 4),          # exposure-response: natural spline
  arglag = list(fun = "ns", df = 4)           # lag-response:      natural spline
)

# ---- Cross-basis for rainfall ----------------------------------------------
cb.rain <- crossbasis(
  df$precip_total, lag = maxlag,
  argvar = list(fun = "ns", df = 4),
  arglag = list(fun = "ns", df = 4)
)

# ---- Cross-bases for large-scale drivers: ENSO (ONI) and IOD (DMI) ----------
cb.oni <- crossbasis(
  df$oni, lag = maxlag,
  argvar = list(fun = "ns", df = 4),
  arglag = list(fun = "ns", df = 4)
)
cb.dmi <- crossbasis(
  df$dmi, lag = maxlag,
  argvar = list(fun = "ns", df = 4),
  arglag = list(fun = "ns", df = 4)
)

# ---- Confounders: seasonality (cyclic) + long-term trend -------------------
df$time <- seq_len(nrow(df))
seas  <- harmonic(df$month, nfreq = 2, period = 12)   # annual + semi-annual
trend <- ns(df$time, df = 6)

# ---- Negative-binomial GLM (overdispersed counts) --------------------------
model <- glm.nb(
  cases ~ cb.temp + cb.rain + cb.oni + cb.dmi + seas + trend + offset(log(population)),
  data = df
)
summary(model)

# ---- Predictions -----------------------------------------------------------
# Cumulative exposure-response, centred at the median exposure
pred.temp <- crosspred(cb.temp, model, cen = median(df$temp_mean), cumul = TRUE)
pred.rain <- crosspred(cb.rain, model, cen = median(df$precip_total), cumul = TRUE)

dir.create("figures", showWarnings = FALSE)

# Overall cumulative exposure-response curves
png("figures/R_er_temperature.png", width = 1400, height = 900, res = 180)
plot(pred.temp, "overall", xlab = "Mean temperature (C)",
     ylab = "Cumulative RR (0-6 mo)", main = "Temperature-dengue (East Java)")
dev.off()

png("figures/R_er_rainfall.png", width = 1400, height = 900, res = 180)
plot(pred.rain, "overall", xlab = "Monthly rainfall (mm)",
     ylab = "Cumulative RR (0-6 mo)", main = "Rainfall-dengue (East Java)")
dev.off()

# 3-D exposure-lag-response surface for rainfall
png("figures/R_surface_rainfall.png", width = 1400, height = 1100, res = 180)
plot(pred.rain, xlab = "Rainfall (mm)", zlab = "RR", ylab = "Lag (months)",
     theta = 210, phi = 30, main = "Rainfall exposure-lag-response surface")
dev.off()

cat("\nDone. Figures written to figures/ (R_*.png)\n")
cat("Model AIC:", AIC(model), "\n")
