"""
00_make_synthetic_climate.py
----------------------------
Generate a CLEARLY-LABELLED SYNTHETIC per-province monthly climate series so the
full pipeline (single- and multi-province) runs without any network access.

!!! THIS IS NOT REAL CLIMATE DATA. !!!
It reproduces plausible seasonality for each province (wet season ~Nov-Apr, mean
temperature modulated by latitude/elevation) plus an ENSO-like slow oscillation.

For a REAL analysis run scripts/01_fetch_climate.py (Open-Meteo, no API key),
which overwrites data/processed/climate_provinces_monthly.csv with real data.
"""
import pathlib
import numpy as np
import pandas as pd
from provinces import PROVINCES

rng = np.random.default_rng(42)
dates = pd.date_range("2004-01-01", "2024-07-01", freq="MS")
n = len(dates); month = dates.month.values; t = np.arange(n)
enso = np.sin(2 * np.pi * t / 42 + 0.7)

frames = []
for i, (prov, (lat, lon)) in enumerate(PROVINCES.items()):
    wet = np.cos(2 * np.pi * (month - 1) / 12)     # peaks ~January
    season = np.sin(2 * np.pi * (month - 1) / 12)
    base_temp = 27.6 + 0.15 * (lat + 7.0)          # mild latitude effect
    amp = 0.9 + 0.05 * i
    temp = base_temp - amp * wet + 0.4 * season + 0.5 * enso + rng.normal(0, 0.3, n)
    precip = np.clip(210 + (160 + 8 * i) * wet - 40 * season - 60 * enso
                     + rng.normal(0, 35, n), 5, None)
    rh = np.clip(80 + 6 * wet - 2 * season - 2 * enso + rng.normal(0, 2, n), 55, 98)
    frames.append(pd.DataFrame({
        "province": prov, "date": dates,
        "temp_mean": np.round(temp, 2),
        "temp_max": np.round(temp + 4.5, 2),
        "temp_min": np.round(temp - 3.5, 2),
        "precip_total": np.round(precip, 1),
        "rh_mean": np.round(rh, 1),
        "SYNTHETIC": True,
    }))

df = pd.concat(frames, ignore_index=True)
out = pathlib.Path(__file__).resolve().parents[1] / "data" / "processed" / "climate_provinces_monthly.csv"
out.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(out, index=False)
print(f"Wrote SYNTHETIC demo climate for {df.province.nunique()} provinces "
      f"({len(df)} rows) -> {out}")

# ---- Large-scale climate indices (national, one series): ENSO + IOD ----------
# Synthetic ONI reuses the ENSO oscillation (typical range ~ -2.5..+2.5);
# synthetic DMI (IOD) is a distinct oscillation (~ -1..+1). CLEARLY LABELLED.
oni = np.round(2.1 * enso + rng.normal(0, 0.15, n), 2)
dmi = np.round(0.6 * np.sin(2 * np.pi * t / 30 + 2.1) + rng.normal(0, 0.08, n), 2)
idx = pd.DataFrame({"date": dates, "oni": oni, "dmi": dmi, "SYNTHETIC": True})
idx_out = out.parent / "climate_indices_monthly.csv"
idx.to_csv(idx_out, index=False)
print(f"Wrote SYNTHETIC ENSO (ONI) + IOD (DMI) indices ({len(idx)} rows) -> {idx_out}")
print("Replace all with real data via scripts/01_fetch_climate.py and 01b_fetch_climate_indices.py")
