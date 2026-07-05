"""
05_forecast_eval.py
-------------------
Out-of-sample forecast-skill evaluation, turning the association model into an
early-warning prototype.

Design (temporal hold-out):
  - Train on 2004-2019, forecast 2020-2024 (never seen in fitting).
  - Compare TWO models on the held-out period:
      * CLIMATE : NB-DLNM with lagged temperature + rainfall + seasonality + trend
      * BASELINE: seasonal "climatology" (same, but WITHOUT climate)
    The question an early-warning system must answer is exactly this: does
    knowing the climate improve the forecast over simply knowing the season?
  - Skill metrics: RMSE, MAE, correlation, and RMSE skill score vs baseline.
  - Outbreak alerts: flag months whose incidence exceeds the training 75th
    percentile; evaluate hit rate / false-alarm rate for both models.

Note on lead time: climate enters at 0-6 month lags, so the climate driving a
given month is already observed. In operational use this yields genuine lead
time; here we use reanalysis climate to isolate the value of that lagged signal.

Caveat: 2020-2021 dengue surveillance was disrupted by COVID-19; the test period
is therefore a hard, realistic stress test.

Outputs:
  figures/fc_fig1_eastjava_forecast.png   observed vs forecast (train/test)
  figures/fc_fig2_alerts.png              outbreak-alert performance
  writeup/forecast_metrics.txt            skill table (East Java + all provinces)
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

MAXLAG, VAR_DF, LAG_DF = 6, 4, 3
SPLIT = pd.Timestamp("2020-01-01")


# ---- crossbasis (knots locked on TRAINING exposure to avoid leakage) --------

def value_basis_fitter(x_train, df):
    d = dmatrix(f"cr(v, df={df}) - 1", {"v": x_train}, return_type="dataframe")
    info = d.design_info
    return (lambda col: np.asarray(build_design_matrices([info], {"v": col})[0])), d.shape[1]


def crossbasis(x_full, train_mask, maxlag, var_df, lag_df):
    x_full = np.asarray(x_full, float)
    vtrans, nv = value_basis_fitter(x_full[train_mask], var_df)
    lags = np.arange(maxlag + 1)
    C = np.asarray(dmatrix(f"cr(l, df={lag_df}) - 1", {"l": lags}, return_type="dataframe"))
    nl = C.shape[1]
    n = len(x_full)
    Q = np.full((n, maxlag + 1), np.nan)
    for l in range(maxlag + 1):
        Q[l:, l] = x_full[: n - l] if l > 0 else x_full
    filled = np.where(np.isnan(Q), np.nanmean(x_full[train_mask]), Q)
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


def seasonal_trend(df):
    seas = np.asarray(dmatrix("cc(month, df=4) - 1", {"month": df.month.values}, return_type="dataframe"))
    trend = df.time_index.values.reshape(-1, 1).astype(float)   # LINEAR trend (extrapolatable)
    trend = (trend - trend.mean()) / trend.std()
    return seas, trend


def fit_predict(df, use_climate, use_indices=False):
    """Fit on train, predict full series. Returns predicted counts (full)."""
    train = df.date < SPLIT
    seas, trend = seasonal_trend(df)
    parts = [np.ones((len(df), 1)), seas, trend]
    if use_climate:
        cbt = crossbasis(df.temp_mean.values, train.values, MAXLAG, VAR_DF, LAG_DF)
        cbr = crossbasis(df.precip_total.values, train.values, MAXLAG, VAR_DF, LAG_DF)
        parts += [cbt, cbr]
    if use_indices:
        # Parsimonious: large-scale indices enter as linear terms at a few lags
        # (0, 3, 6 months) rather than full cross-bases, to add lead-time signal
        # without the overfitting that many spline parameters would cause.
        def lagged(v, L):
            out = np.full(len(v), np.nan); out[L:] = v[:len(v) - L] if L else v
            return out.reshape(-1, 1)
        oni = df.oni.values.astype(float); dmi = df.dmi.values.astype(float)
        idx_cols = [lagged(oni, L) for L in (0, 3, 6)] + [lagged(dmi, L) for L in (0, 3, 6)]
        parts += idx_cols
    X = np.hstack(parts)
    off = np.log(df.population.values.astype(float))
    y = df.cases.values.astype(float)

    valid = ~np.isnan(X).any(axis=1)
    tr = valid & train.values
    Xtr, ytr, otr = X[tr], y[tr], off[tr]

    poi = sm.GLM(ytr, Xtr, family=sm.families.Poisson(), offset=otr).fit()
    mu = poi.fittedvalues
    alpha = max(((ytr - mu) ** 2 - mu).sum() / (mu ** 2).sum(), 1e-3)
    nb = sm.GLM(ytr, Xtr, family=sm.families.NegativeBinomial(alpha=alpha), offset=otr).fit()

    pred = np.full(len(df), np.nan)
    pred[valid] = nb.predict(X[valid], offset=off[valid])
    return pred


def skill(df, pred_clim, pred_base):
    test = (df.date >= SPLIT).values & ~np.isnan(pred_clim) & ~np.isnan(pred_base)
    obs = df.cases.values[test].astype(float)
    pc, pb = pred_clim[test], pred_base[test]

    def rmse(a, b): return float(np.sqrt(np.mean((a - b) ** 2)))
    def mae(a, b): return float(np.mean(np.abs(a - b)))
    r_c = float(np.corrcoef(obs, pc)[0, 1]); r_b = float(np.corrcoef(obs, pb)[0, 1])
    rmse_c, rmse_b = rmse(obs, pc), rmse(obs, pb)
    return dict(n=int(test.sum()), rmse_clim=rmse_c, rmse_base=rmse_b,
                mae_clim=mae(obs, pc), mae_base=mae(obs, pb),
                r_clim=r_c, r_base=r_b,
                rmse_skill=1 - rmse_c / rmse_b)


def corr_rmse(df, pred):
    test = (df.date >= SPLIT).values & ~np.isnan(pred)
    obs = df.cases.values[test].astype(float); p = pred[test]
    return float(np.corrcoef(obs, p)[0, 1]), float(np.sqrt(np.mean((obs - p) ** 2)))


# ---- East Java deep-dive ----------------------------------------------------
ej = pd.read_csv(ROOT / "data/processed/analysis_eastjava.csv", parse_dates=["date"])
ej = ej.dropna(subset=["cases", "temp_mean", "precip_total", "population", "oni", "dmi"]).reset_index(drop=True)
synthetic = "SYNTHETIC" in ej.columns and bool(ej["SYNTHETIC"].iloc[0])
TAG = "  [SYNTHETIC DEMO CLIMATE]" if synthetic else ""

pred_base = fit_predict(ej, use_climate=False)
pred_clim = fit_predict(ej, use_climate=True)
pred_full = fit_predict(ej, use_climate=True, use_indices=True)   # climate + ENSO/IOD
m = skill(ej, pred_clim, pred_base)
r_full, rmse_full = corr_rmse(ej, pred_full)

print(f"East Java forecast skill (test {SPLIT.year}-2024, n={m['n']} months)")
print(f"  RMSE  baseline={m['rmse_base']:.0f}  climate={m['rmse_clim']:.0f}  "
      f"climate+indices={rmse_full:.0f}")
print(f"  corr  baseline={m['r_base']:.2f}  climate={m['r_clim']:.2f}  "
      f"climate+indices={r_full:.2f}")

# ---- Figure 1: observed vs forecast ----------------------------------------
plt.rcParams.update({"figure.dpi": 130, "font.size": 10, "axes.grid": True,
                     "grid.alpha": 0.3, "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(figsize=(10, 4.3))
ax.axvspan(SPLIT, ej.date.max(), color="#f0f0f0", label="Test (held-out)")
ax.plot(ej.date, ej.cases, color="#333", lw=1.3, label="Observed")
ax.plot(ej.date, pred_base, color="#8aa0b3", lw=1.4, ls="--", label="Baseline (seasonal)")
ax.plot(ej.date, pred_clim, color="#1f6fb3", lw=1.5, label="Climate DLNM")
ax.plot(ej.date, pred_full, color="#b3123f", lw=1.8, label="Climate + ENSO/IOD")
ax.axvline(SPLIT, color="k", lw=0.8)
ax.set_ylabel("Dengue cases"); ax.set_xlabel("Date")
ax.set_title(f"East Java: out-of-sample forecast, trained on 2004-2019{TAG}")
ax.legend(fontsize=8, ncol=2)
fig.tight_layout(); fig.savefig(FIG / "fc_fig1_eastjava_forecast.png"); plt.close(fig)

# ---- Figure 2: outbreak-alert performance ----------------------------------
thr = np.percentile(ej.loc[ej.date < SPLIT, "incidence_per_100k"], 75)
inc = ej.incidence_per_100k.values
pop = ej.population.values
test = (ej.date >= SPLIT).values & ~np.isnan(pred_clim)
obs_alert = inc[test] > thr
def alert_from_pred(pred):
    pinc = 100_000 * pred[test] / pop[test]
    return pinc > thr

def rates(pred_alert):
    tp = np.sum(pred_alert & obs_alert); fn = np.sum(~pred_alert & obs_alert)
    fp = np.sum(pred_alert & ~obs_alert); tn = np.sum(~pred_alert & ~obs_alert)
    sens = tp / (tp + fn) if (tp + fn) else np.nan
    spec = tn / (tn + fp) if (tn + fp) else np.nan
    return sens, spec

sens_c, spec_c = rates(alert_from_pred(pred_clim))
sens_b, spec_b = rates(alert_from_pred(pred_base))

fig, ax = plt.subplots(figsize=(7, 4.2))
labels = ["Sensitivity\n(hit rate)", "Specificity"]
x = np.arange(2); w = 0.35
ax.bar(x - w/2, [sens_c, spec_c], w, label="Climate DLNM", color="#b3123f")
ax.bar(x + w/2, [sens_b, spec_b], w, label="Baseline", color="#8aa0b3")
ax.set_xticks(x); ax.set_xticklabels(labels); ax.set_ylim(0, 1.05)
ax.set_ylabel("Rate"); ax.set_title(f"Outbreak-alert performance, test period{TAG}")
ax.legend(fontsize=9)
for i, (a, b) in enumerate(zip([sens_c, spec_c], [sens_b, spec_b])):
    ax.text(i - w/2, a + 0.02, f"{a:.2f}", ha="center", fontsize=8)
    ax.text(i + w/2, b + 0.02, f"{b:.2f}", ha="center", fontsize=8)
fig.tight_layout(); fig.savefig(FIG / "fc_fig2_alerts.png"); plt.close(fig)

# ---- Multi-province skill table --------------------------------------------
allp = pd.read_csv(ROOT / "data/processed/analysis_provinces.csv", parse_dates=["date"])
rows = []
for prov, g in allp.groupby("province"):
    g = g.dropna(subset=["cases", "temp_mean", "precip_total", "population", "oni", "dmi"]).reset_index(drop=True)
    try:
        pb = fit_predict(g, False)
        pc = fit_predict(g, True)
        pf = fit_predict(g, True, True)
        s = skill(g, pc, pb)
        r_f, _ = corr_rmse(g, pf)
        rows.append((prov, s["r_base"], s["r_clim"], r_f, s["rmse_skill"] * 100))
    except Exception:
        rows.append((prov, np.nan, np.nan, np.nan, np.nan))
tab = pd.DataFrame(rows, columns=["province", "r_baseline", "r_climate", "r_clim+idx", "rmse_skill_%"])

with open(ROOT / "writeup" / "forecast_metrics.txt", "w") as f:
    f.write("Out-of-sample forecast skill (train 2004-2019, test 2020-2024)\n")
    f.write("=" * 62 + "\n")
    f.write(f"Climate source: {'SYNTHETIC DEMO' if synthetic else 'REAL (Open-Meteo)'}\n\n")
    f.write("East Java:\n")
    f.write(f"  RMSE baseline={m['rmse_base']:.0f}  climate={m['rmse_clim']:.0f}  "
            f"climate+indices={rmse_full:.0f}\n")
    f.write(f"  corr baseline={m['r_base']:.2f}  climate={m['r_clim']:.2f}  "
            f"climate+indices={r_full:.2f}\n")
    f.write(f"  outbreak sensitivity climate={sens_c:.2f} baseline={sens_b:.2f}; "
            f"specificity climate={spec_c:.2f} baseline={spec_b:.2f}\n\n")
    f.write("Per-province correlation with held-out cases (baseline / climate / climate+ENSO-IOD)\n")
    f.write("and RMSE skill of climate vs baseline (+ve = climate helps):\n")
    f.write(tab.to_string(index=False, float_format=lambda v: f"{v:.2f}"))
    f.write("\n")

print("\nPer-province skill (RMSE skill %, +ve = climate improves forecast):")
print(tab.to_string(index=False, float_format=lambda v: f"{v:.2f}"))
print("\nFigures -> figures/fc_*.png ; metrics -> writeup/forecast_metrics.txt")
