"""
03_dlnm_model.py
----------------
Distributed-lag non-linear model (DLNM) of the association between monthly
climate (temperature, rainfall) and dengue incidence in East Java, Indonesia.

Framework (Gasparrini 2011; Lowe et al. 2021, Lancet Planet Health):
  - A "cross-basis" jointly models the non-linear exposure-response AND the
    lagged (delayed) effect of each climate variable over 0-L months.
  - Counts are modelled with a negative-binomial GLM to allow overdispersion,
    adjusting for seasonality (cyclic harmonics) and long-term trend.

This is a self-contained Python implementation of the crossbasis so it runs
without R. For the publication-grade version see analysis/dlnm_model.R, which
uses the canonical `dlnm` R package.

Outputs (in figures/):
  fig1_timeseries.png          observed dengue + climate
  fig2_crosscorr.png           cross-correlation (dengue vs lagged climate)
  fig3_er_temperature.png      cumulative exposure-response, temperature
  fig4_er_rainfall.png         cumulative exposure-response, rainfall
  fig5_fit.png                 observed vs fitted
Also writes writeup/model_summary.txt
"""

import pathlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from patsy import dmatrix, build_design_matrices
import statsmodels.api as sm

ROOT = pathlib.Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)

MAXLAG = 6        # months of lag to consider
VAR_DF = 4        # spline df for the exposure-response dimension
LAG_DF = 4        # spline df for the lag-response dimension

# ----------------------------------------------------------------------------
# Crossbasis machinery (faithful to dlnm's construction)
# ----------------------------------------------------------------------------

def lag_matrix(x, maxlag):
    """Return N x (maxlag+1) matrix; column l holds x shifted down by l months."""
    n = len(x)
    Q = np.full((n, maxlag + 1), np.nan)
    for l in range(maxlag + 1):
        Q[l:, l] = x[: n - l] if l > 0 else x
    return Q


def make_value_basis(x, df):
    """Natural cubic spline basis with knots locked to x's distribution."""
    design = dmatrix(f"cr(v, df={df}) - 1", {"v": x}, return_type="dataframe")
    info = design.design_info

    def transform(col):
        return np.asarray(build_design_matrices([info], {"v": col})[0])

    return transform, design.shape[1]


def crossbasis(x, maxlag, var_df, lag_df):
    """Build the DLNM cross-basis for exposure series x."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    Q = lag_matrix(x, maxlag)

    # value (exposure-response) basis, knots from the observed distribution
    vtrans, nv = make_value_basis(x, var_df)

    # lag-response basis over lags 0..L
    lags = np.arange(maxlag + 1)
    lag_design = dmatrix(f"cr(l, df={lag_df}) - 1", {"l": lags}, return_type="dataframe")
    C = np.asarray(lag_design)            # (L+1) x nl
    nl = C.shape[1]

    # value basis applied to each lagged column (fill NaN with a neutral value
    # first, then re-mask so the early rows are dropped consistently)
    filled = np.where(np.isnan(Q), np.nanmean(x), Q)
    vb = [vtrans(filled[:, l]) for l in range(maxlag + 1)]   # each n x nv

    cb = np.zeros((n, nv * nl))
    for j in range(nv):
        for k in range(nl):
            col = np.zeros(n)
            for l in range(maxlag + 1):
                col += vb[l][:, j] * C[l, k]
            cb[:, j * nl + k] = col

    cb[:maxlag, :] = np.nan   # first L rows have incomplete lag history
    meta = dict(nv=nv, nl=nl, C=C, vtrans=vtrans, maxlag=maxlag,
                xref=float(np.median(x)), xrange=(float(np.nanmin(x)), float(np.nanmax(x))))
    return cb, meta


def crosspred_cumulative(meta, beta, vcov, values):
    """Cumulative (summed over lags) exposure-response relative to xref.

    Returns RR and 95% CI at each value in `values`.
    """
    nv, nl, C = meta["nv"], meta["nl"], meta["C"]
    vtrans = meta["vtrans"]
    lag_sum = C.sum(axis=0)                       # length nl (sum of lag basis)

    # basis of the reference value
    bref = vtrans(np.array([meta["xref"]]))[0]    # length nv

    logrr, se = [], []
    for v in values:
        bv = vtrans(np.array([float(v)]))[0]      # length nv
        # cumulative cross-basis row = (bv - bref) outer lag_sum, flattened
        row = np.outer(bv - bref, lag_sum).reshape(-1)   # length nv*nl
        logrr.append(row @ beta)
        se.append(np.sqrt(row @ vcov @ row))
    logrr = np.array(logrr); se = np.array(se)
    return np.exp(logrr), np.exp(logrr - 1.96 * se), np.exp(logrr + 1.96 * se)


# ----------------------------------------------------------------------------
# Load data
# ----------------------------------------------------------------------------

df = pd.read_csv(ROOT / "data" / "processed" / "analysis_eastjava.csv", parse_dates=["date"])
df = df.dropna(subset=["cases", "temp_mean", "precip_total", "population", "oni", "dmi"]).reset_index(drop=True)
synthetic = "SYNTHETIC" in df.columns and bool(df["SYNTHETIC"].iloc[0]) if "SYNTHETIC" in df.columns else False

print(f"Analysis window: {df.date.min().date()} -> {df.date.max().date()}  ({len(df)} months)")
print(f"Climate source : {'SYNTHETIC DEMO' if synthetic else 'REAL (Open-Meteo)'}")

# ----------------------------------------------------------------------------
# Build cross-bases and design matrix
# ----------------------------------------------------------------------------

cb_temp, meta_temp = crossbasis(df["temp_mean"], MAXLAG, VAR_DF, LAG_DF)
cb_rain, meta_rain = crossbasis(df["precip_total"], MAXLAG, VAR_DF, LAG_DF)
cb_oni, meta_oni = crossbasis(df["oni"], MAXLAG, VAR_DF, LAG_DF)
cb_dmi, meta_dmi = crossbasis(df["dmi"], MAXLAG, VAR_DF, LAG_DF)

# seasonality (cyclic) + long-term trend
seas = dmatrix("cc(month, df=4) - 1", {"month": df["month"].values}, return_type="dataframe")
trend = dmatrix("cr(t, df=6) - 1", {"t": df["time_index"].values}, return_type="dataframe")

Xparts = [np.ones((len(df), 1)), cb_temp, cb_rain, cb_oni, cb_dmi,
          np.asarray(seas), np.asarray(trend)]
X = np.hstack(Xparts)
y = df["cases"].values.astype(float)
offset = np.log(df["population"].values.astype(float))   # incidence model offset

# drop rows with NaN (first MAXLAG months)
mask = ~np.isnan(X).any(axis=1)
Xf, yf, off_f = X[mask], y[mask], offset[mask]

# column index bookkeeping
i0 = 1
i_temp = slice(i0, i0 + cb_temp.shape[1]); i0 += cb_temp.shape[1]
i_rain = slice(i0, i0 + cb_rain.shape[1]); i0 += cb_rain.shape[1]
i_oni = slice(i0, i0 + cb_oni.shape[1]); i0 += cb_oni.shape[1]
i_dmi = slice(i0, i0 + cb_dmi.shape[1]); i0 += cb_dmi.shape[1]

# ----------------------------------------------------------------------------
# Fit negative-binomial GLM (estimate alpha via a quick Poisson->NB step)
# ----------------------------------------------------------------------------

poi = sm.GLM(yf, Xf, family=sm.families.Poisson(), offset=off_f).fit()
mu = poi.fittedvalues
# method-of-moments dispersion for NB alpha
alpha = max(((yf - mu) ** 2 - mu).sum() / (mu ** 2).sum(), 1e-3)
nb = sm.GLM(yf, Xf, family=sm.families.NegativeBinomial(alpha=alpha), offset=off_f).fit()

print(f"\nNegative-binomial GLM fitted. alpha={alpha:.3f}, "
      f"pseudo-R2(dev)={1 - nb.deviance / nb.null_deviance:.3f}")

beta = nb.params
vcov = nb.cov_params()

# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------
plt.rcParams.update({"figure.dpi": 130, "font.size": 10, "axes.grid": True,
                     "grid.alpha": 0.3, "axes.spines.top": False, "axes.spines.right": False})
TAG = "  [SYNTHETIC DEMO CLIMATE]" if synthetic else ""

# Fig 1: time series
fig, ax = plt.subplots(3, 1, figsize=(9, 7), sharex=True)
ax[0].plot(df.date, df.cases, color="#b3123f"); ax[0].set_ylabel("Dengue cases")
ax[0].set_title("East Java: monthly dengue (OpenDengue v1.3) and climate" + TAG)
ax[1].plot(df.date, df.temp_mean, color="#d1651a"); ax[1].set_ylabel("Mean temp (C)")
ax[2].plot(df.date, df.precip_total, color="#1f6fb3"); ax[2].set_ylabel("Rainfall (mm)")
ax[2].set_xlabel("Date")
fig.tight_layout(); fig.savefig(FIG / "fig1_timeseries.png"); plt.close(fig)

# Fig 2: cross-correlation of dengue with lagged climate
def ccf_lagged(a, b, maxlag):
    a = (a - np.mean(a)) / np.std(a); b = (b - np.mean(b)) / np.std(b)
    return [np.corrcoef(a[l:], b[:len(b) - l])[0, 1] if l > 0 else np.corrcoef(a, b)[0, 1]
            for l in range(maxlag + 1)]

lags = np.arange(MAXLAG + 1)
cc_t = ccf_lagged(df.cases.values, df.temp_mean.values, MAXLAG)
cc_r = ccf_lagged(df.cases.values, df.precip_total.values, MAXLAG)
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(lags, cc_t, "o-", label="Temperature", color="#d1651a")
ax.plot(lags, cc_r, "s-", label="Rainfall", color="#1f6fb3")
ax.axhline(0, color="k", lw=0.8)
ax.set_xlabel("Lag (months): dengue vs climate that many months earlier")
ax.set_ylabel("Correlation"); ax.set_title("Cross-correlation" + TAG); ax.legend()
fig.tight_layout(); fig.savefig(FIG / "fig2_crosscorr.png"); plt.close(fig)

# Fig 3 & 4: cumulative exposure-response
def er_plot(meta, colslice, values, label, unit, color, fname):
    b = beta[colslice]; V = vcov[colslice, colslice]
    rr, lo, hi = crosspred_cumulative(meta, b, V, values)
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.fill_between(values, lo, hi, color=color, alpha=0.2)
    ax.plot(values, rr, color=color, lw=2)
    ax.axhline(1, color="k", lw=0.8, ls="--")
    ax.axvline(meta["xref"], color="grey", lw=0.8, ls=":")
    ax.set_xlabel(f"{label} ({unit})")
    ax.set_ylabel("Cumulative RR (0-%d mo)" % MAXLAG)
    ax.set_title(f"Exposure-response: {label} vs dengue" + TAG)
    fig.tight_layout(); fig.savefig(FIG / fname); plt.close(fig)

tvals = np.linspace(*meta_temp["xrange"], 60)
rvals = np.linspace(*meta_rain["xrange"], 60)
er_plot(meta_temp, i_temp, tvals, "Mean temperature", "C", "#d1651a", "fig3_er_temperature.png")
er_plot(meta_rain, i_rain, rvals, "Monthly rainfall", "mm", "#1f6fb3", "fig4_er_rainfall.png")

# ENSO (ONI) and IOD (DMI) exposure-response: large-scale drivers
def er_index_plot(meta, colslice, label, color, fname, xlab):
    b = beta[colslice]; V = vcov[colslice, colslice]
    values = np.linspace(*meta["xrange"], 60)
    rr, lo, hi = crosspred_cumulative(meta, b, V, values)
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.fill_between(values, lo, hi, color=color, alpha=0.2)
    ax.plot(values, rr, color=color, lw=2)
    ax.axhline(1, color="k", lw=0.8, ls="--")
    ax.axvline(0, color="grey", lw=0.8, ls=":")
    ax.set_xlabel(xlab); ax.set_ylabel("Cumulative RR (0-%d mo)" % MAXLAG)
    ax.set_title(f"Exposure-response: {label} vs dengue" + TAG)
    fig.tight_layout(); fig.savefig(FIG / fname); plt.close(fig)

er_index_plot(meta_oni, i_oni, "ENSO (ONI)", "#7a2f8f", "fig6_er_oni.png",
              "Oceanic Nino Index (< -0.5 La Nina, > 0.5 El Nino)")
er_index_plot(meta_dmi, i_dmi, "IOD (DMI)", "#2f8f6b", "fig7_er_dmi.png",
              "Dipole Mode Index (+ positive IOD)")

# Fig 5: observed vs fitted
fitted = np.full(len(df), np.nan); fitted[mask] = nb.fittedvalues
fig, ax = plt.subplots(figsize=(9, 4))
ax.plot(df.date, df.cases, color="#333", lw=1.2, label="Observed")
ax.plot(df.date, fitted, color="#b3123f", lw=1.6, label="Fitted (NB-DLNM)")
ax.set_ylabel("Dengue cases"); ax.set_xlabel("Date")
ax.set_title("Model fit" + TAG); ax.legend()
fig.tight_layout(); fig.savefig(FIG / "fig5_fit.png"); plt.close(fig)

# ----------------------------------------------------------------------------
# Text summary
# ----------------------------------------------------------------------------
peak_t_lag = int(np.argmax(cc_t)); peak_r_lag = int(np.argmax(cc_r))
with open(ROOT / "writeup" / "model_summary.txt", "w") as f:
    f.write("DLNM negative-binomial model - East Java dengue vs climate\n")
    f.write("=" * 60 + "\n")
    f.write(f"Climate source: {'SYNTHETIC DEMO' if synthetic else 'REAL (Open-Meteo)'}\n")
    f.write(f"Months analysed: {int(mask.sum())} "
            f"({df.date.min().date()} to {df.date.max().date()})\n")
    f.write(f"Max lag: {MAXLAG} months | var df: {VAR_DF} | lag df: {LAG_DF}\n")
    f.write(f"NB dispersion alpha: {alpha:.3f}\n")
    f.write(f"Deviance pseudo-R2: {1 - nb.deviance / nb.null_deviance:.3f}\n")
    f.write(f"Strongest temp cross-correlation at lag {peak_t_lag} mo "
            f"(r={cc_t[peak_t_lag]:.2f})\n")
    f.write(f"Strongest rain cross-correlation at lag {peak_r_lag} mo "
            f"(r={cc_r[peak_r_lag]:.2f})\n")

print("\nFigures written to figures/. Summary -> writeup/model_summary.txt")
print(f"Strongest temp cc at lag {peak_t_lag} mo (r={cc_t[peak_t_lag]:.2f}); "
      f"rain at lag {peak_r_lag} mo (r={cc_r[peak_r_lag]:.2f})")
