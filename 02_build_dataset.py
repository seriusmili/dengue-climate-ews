"""
02_build_dataset.py
-------------------
Merge OpenDengue cases + population + climate into analysis-ready tables:
  data/processed/analysis_eastjava.csv     (single-province, with incidence)
  data/processed/analysis_provinces.csv    (all provinces, with incidence)
Population (BPS) is interpolated to monthly and used as a model offset.
"""
import pathlib
import pandas as pd
from provinces import PRIMARY

ROOT = pathlib.Path(__file__).resolve().parents[1]

dengue = pd.read_csv(ROOT / "data/raw/dengue_provinces_opendengue_v1.3.csv", parse_dates=["date"])
pop = pd.read_csv(ROOT / "data/raw/population_provinces.csv")
clim = pd.read_csv(ROOT / "data/processed/climate_provinces_monthly.csv", parse_dates=["date"])
idx = pd.read_csv(ROOT / "data/processed/climate_indices_monthly.csv", parse_dates=["date"])

dengue["year"] = dengue.date.dt.year
df = dengue.merge(pop, on=["province", "year"], how="left")
df = df.merge(clim, on=["province", "date"], how="inner")
df = df.merge(idx[["date", "oni", "dmi"]], on="date", how="left")   # national indices
df["month"] = df.date.dt.month
df["incidence_per_100k"] = 100_000 * df.cases / df.population
df = df.sort_values(["province", "date"]).reset_index(drop=True)
df["time_index"] = df.groupby("province").cumcount()

full = ROOT / "data/processed/analysis_provinces.csv"
df.to_csv(full, index=False)
print(f"Multi-province dataset: {df.province.nunique()} provinces, {len(df)} rows -> {full}")

ej = df[df.province == PRIMARY].reset_index(drop=True)
ejp = ROOT / "data/processed/analysis_eastjava.csv"
ej.to_csv(ejp, index=False)
print(f"{PRIMARY}: {len(ej)} months -> {ejp}")
print(ej[["date", "cases", "population", "incidence_per_100k", "temp_mean", "precip_total"]].head())
