"""
04_multiprovince.py
-------------------
Pooled multi-province distributed-lag non-linear model (DLNM) across seven
Indonesian provinces (Java + Bali), modelling dengue INCIDENCE with a
population offset. A shared climate cross-basis (temperature, rainfall) is
estimated while allowing province-specific baseline risk, seasonality and
long-term trend.

This is the natural extension of the single-province analysis toward a
hierarchical early-warning framework. (A fully hierarchical Bayesian version
would use INLA / mvmeta; here we use a pooled NB-GLM with province fixed
effects, which is transparent and fast.)

Outputs (figures/):
  mp_fig1_incidence_panel.png   monthly incidence by province
  mp_fig2_er_rainfall.png       pooled cumulative exposure-response, rainfall
  mp_fig3_er_temperature.png    pooled cumulative exposure-response, temperature
  mp_fig4_heterogeneity.png     province-specific peak rainfall-lag correlation
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
FIG = ROOT / "figures"; FIG.mkdir(exist_ok=True)

MAXLAG, VAR_DF, LAG_DF = 6, 4, 4

# ---- crossbasis helpers (shared value-basis knots; per-province lag blocks) --

def make_value_basis(x, df):
    d = dmatrix(f"cr(v, df={df}) - 1", {"v": x}, return_type="dataframe")
    info = d.design_info
    return (lambda col: np.asarray(build_design_matrices([info], {"v": col})[0])), d.shape[1]


def lag_block(x, vtrans, C, maxlag):
    """Cross-basis block for a single province's series x."""
    x = np.asarray(x, float); n = len(x); nv = vtrans(x[:1]).shape[1]; nl = C.shape[1]
    Q = np.full((n, maxlag + 1), np.nan)
    for l in range(maxlag + 1):
        Q[l:, l] = x[: n - l] if l > 0 else x
    filled = np.where(np.isnan(Q), np.nanmean(x), Q)
    vb = [vtrans(filled[:, l]) for l in range(maxlag + 1)]
    cb = np.zeros((n, nv * nl))
    for j in range(nv):
        for k in range(nl):
            col = np.zeros(n)
            for l in range(maxlag + 1):
                col += vb[l][:, j] * C[l, k]
            cb[:, j * nl + k] = col
    cb[:maxlag, :] = np.nan
    return cb


def build_cb(df, var, maxlag, var_df, lag_df):
    vtrans, nv = make_value_basis(df[var].values, var_df)
    lags = np.arange(maxlag + 1)
    C = np.asarray(dmatrix(f"cr(l, df={lag_df}) - 1", {"l": lags}, return_type="dataframe"))
    blocks = [lag_block(g[var].values, vtrans, C, maxlag) for _, g in df.groupby("province", sort=False)]
    cb = np.vstack(blocks)
    meta = dict(nv=nv, nl=C.shape[1], C=C, vtrans=vtrans,
                xref=float(np.median(df[var])), xrange=(float(df[var].min()), float(df[var].max())))
    return cb, meta


def crosspred_cumulative(meta, beta, vcov, values):
    C, vtrans = meta["C"], meta["vtrans"]
    lag_sum = C.sum(axis=0)
    bref = vtrans(np.array([meta["xref"]]))[0]
    rr, lo, hi = [], [], []
    for v in values:
        bv = vtrans(np.array([float(v)]))[0]
        row = np.outer(bv - bref, lag_sum).reshape(-1)
        m = row @ beta; s = np.sqrt(row @ vcov @ row)
        rr.append(np.exp(m)); lo.append(np.exp(m - 1.96 * s)); hi.append(np.exp(m + 1.96 * s))
    return np.array(rr), np.array(lo), np.array(hi)


# ---- load -------------------------------------------------------------------
df = pd.read_csv(ROOT / "data/processed/analysis_provinces.csv", parse_dates=["date"])
df = df.dropna(subset=["cases", "temp_mean", "precip_total", "population"]).copy()
df = df.sort_values(["province", "date"]).reset_index(drop=True)
synthetic = "SYNTHETIC" in df.columns and bool(df["SYNTHETIC"].iloc[0])
provs = list(df.province.unique())
print(f"Provinces: {provs}")
print(f"Climate source: {'SYNTHETIC DEMO' if synthetic else 'REAL (Open-Meteo)'}")

# ---- design -----------------------------------------------------------------
cb_temp, meta_temp = build_cb(df, "temp_mean", MAXLAG, VAR_DF, LAG_DF)
cb_rain, meta_rain = build_cb(df, "precip_total", MAXLAG, VAR_DF, LAG_DF)
seas = np.asarray(dmatrix("cc(month, df=4) - 1", {"month": df.month.values}, return_type="dataframe"))
prov_fe = np.asarray(dmatrix("C(province) - 1", {"province": df.province.values}, return_type="dataframe"))
trend = np.asarray(dmatrix("cr(t, df=6) - 1", {"t": df.time_index.values}, return_type="dataframe"))

X = np.hstack([cb_temp, cb_rain, seas, prov_fe, trend])
y = df.cases.values.astype(float)
offset = np.log(df.population.values.astype(float))
mask = ~np.isnan(X).any(axis=1)
Xf, yf, off_f = X[mask], y[mask], offset[mask]

i0 = 0
i_temp = slice(i0, i0 + cb_temp.shape[1]); i0 += cb_temp.shape[1]
i_rain = slice(i0, i0 + cb_rain.shape[1]); i0 += cb_rain.shape[1]

poi = sm.GLM(yf, Xf, family=sm.families.Poisson(), offset=off_f).fit()
mu = poi.fittedvalues
alpha = max(((yf - mu) ** 2 - mu).sum() / (mu ** 2).sum(), 1e-3)
nb = sm.GLM(yf, Xf, family=sm.families.NegativeBinomial(alpha=alpha), offset=off_f).fit()
beta, vcov = nb.params, nb.cov_params()
print(f"Pooled NB-DLNM fitted on {int(mask.sum())} province-months. alpha={alpha:.3f}")

# ---- figures ----------------------------------------------------------------
plt.rcParams.update({"figure.dpi": 130, "font.size": 9, "axes.grid": True,
                     "grid.alpha": 0.3, "axes.spines.top": False, "axes.spines.right": False})
TAG = "  [SYNTHETIC DEMO CLIMATE]" if synthetic else ""

# Fig 1: incidence small multiples
fig, axes = plt.subplots(4, 2, figsize=(11, 10), sharex=True)
for ax, prov in zip(axes.ravel(), provs):
    g = df[df.province == prov]
    ax.plot(g.date, g.incidence_per_100k, color="#b3123f", lw=1)
    ax.set_title(prov, fontsize=9); ax.set_ylabel("per 100k")
for ax in axes.ravel()[len(provs):]:
    ax.axis("off")
fig.suptitle("Monthly dengue incidence by province (OpenDengue v1.3)" + TAG, y=1.0)
fig.tight_layout(); fig.savefig(FIG / "mp_fig1_incidence_panel.png"); plt.close(fig)

# Fig 2 & 3: pooled exposure-response
def er_plot(meta, sl, values, label, unit, color, fname):
    rr, lo, hi = crosspred_cumulative(meta, beta[sl], vcov[sl, sl], values)
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.fill_between(values, lo, hi, color=color, alpha=0.2)
    ax.plot(values, rr, color=color, lw=2)
    ax.axhline(1, color="k", lw=0.8, ls="--"); ax.axvline(meta["xref"], color="grey", lw=0.8, ls=":")
    ax.set_xlabel(f"{label} ({unit})"); ax.set_ylabel(f"Pooled cumulative RR (0-{MAXLAG} mo)")
    ax.set_title(f"Pooled exposure-response: {label}" + TAG)
    fig.tight_layout(); fig.savefig(FIG / fname); plt.close(fig)

er_plot(meta_rain, i_rain, np.linspace(*meta_rain["xrange"], 60),
        "Monthly rainfall", "mm", "#1f6fb3", "mp_fig2_er_rainfall.png")
er_plot(meta_temp, i_temp, np.linspace(*meta_temp["xrange"], 60),
        "Mean temperature", "C", "#d1651a", "mp_fig3_er_temperature.png")

# Fig 4: heterogeneity in peak rainfall-lag correlation
def peak_lag(a, b, maxlag):
    a = (a - a.mean()) / a.std(); b = (b - b.mean()) / b.std()
    cc = [np.corrcoef(a[l:], b[:len(b) - l])[0, 1] if l else np.corrcoef(a, b)[0, 1]
          for l in range(maxlag + 1)]
    return int(np.argmax(cc)), max(cc)

rows = []
for prov in provs:
    g = df[df.province == prov]
    lag, r = peak_lag(g.cases.values.astype(float), g.precip_total.values.astype(float), MAXLAG)
    rows.append((prov, lag, r))
het = pd.DataFrame(rows, columns=["province", "peak_lag", "r"]).sort_values("r")
fig, ax = plt.subplots(figsize=(7, 4.2))
bars = ax.barh(het.province, het.r, color="#1f6fb3")
for b_, lag in zip(bars, het.peak_lag):
    ax.text(b_.get_width() + 0.01, b_.get_y() + b_.get_height() / 2,
            f"lag {lag} mo", va="center", fontsize=8)
ax.set_xlabel("Peak rainfall-dengue correlation")
ax.set_title("Between-province heterogeneity in rainfall lag" + TAG)
fig.tight_layout(); fig.savefig(FIG / "mp_fig4_heterogeneity.png"); plt.close(fig)

het.to_csv(ROOT / "writeup" / "province_heterogeneity.csv", index=False)
print("Multi-province figures written to figures/ (mp_*.png)")
print(het.to_string(index=False))
