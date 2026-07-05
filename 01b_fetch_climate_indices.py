"""
01b_fetch_climate_indices.py
----------------------------
Download REAL large-scale climate indices used in dengue early-warning:
  - ONI  (Oceanic Nino Index, ENSO)  -- NOAA CPC
  - DMI  (Dipole Mode Index, Indian Ocean Dipole / IOD) -- NOAA PSL

Both are free and require no API key. Overwrites
data/processed/climate_indices_monthly.csv with real monthly values.

Sources:
  ONI: https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt
       (columns: SEAS YR TOTAL ANOM; ANOM is the ONI. The 3-month season code
        is mapped to its centre month, e.g. DJF -> January.)
  DMI: https://psl.noaa.gov/gcos_wgsp/Timeseries/Data/dmi.had.long.data
       (year followed by 12 monthly values; missing flagged, e.g. -9999.)

Note: as of Feb 2026 NOAA uses RONI for official ENSO monitoring; ONI remains
published and is standard for historical analysis. Swap the URL if you prefer
RONI. IOD (DMI) matters as much as ENSO for Indonesian rainfall.

Usage:  python scripts/01b_fetch_climate_indices.py
"""
import io
import pathlib
import requests
import numpy as np
import pandas as pd

OUT = pathlib.Path(__file__).resolve().parents[1] / "data" / "processed" / "climate_indices_monthly.csv"

SEASON_CENTRE = {  # 3-month season code -> centre month number
    "DJF": 1, "JFM": 2, "FMA": 3, "MAM": 4, "AMJ": 5, "MJJ": 6,
    "JJA": 7, "JAS": 8, "ASO": 9, "SON": 10, "OND": 11, "NDJ": 12,
}


def fetch_oni():
    url = "https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt"
    txt = requests.get(url, timeout=60).text
    df = pd.read_csv(io.StringIO(txt), sep=r"\s+")
    df.columns = [c.upper() for c in df.columns]  # SEAS YR TOTAL ANOM
    df["month"] = df["SEAS"].map(SEASON_CENTRE)
    # NDJ centre month is Dec of the SEAS year's span; keep YR as the year the
    # centre month falls in (DJF centre = Jan of YR for CPC's labelling).
    df["date"] = pd.to_datetime(dict(year=df["YR"], month=df["month"], day=1))
    return df[["date", "ANOM"]].rename(columns={"ANOM": "oni"})


def fetch_dmi():
    url = "https://psl.noaa.gov/gcos_wgsp/Timeseries/Data/dmi.had.long.data"
    txt = requests.get(url, timeout=60).text
    rows = []
    for line in txt.splitlines():
        parts = line.split()
        if len(parts) == 13 and parts[0].isdigit() and len(parts[0]) == 4:
            year = int(parts[0])
            for m, v in enumerate(parts[1:], start=1):
                val = float(v)
                if val > -900:  # missing-value flag guard
                    rows.append((pd.Timestamp(year=year, month=m, day=1), val))
    return pd.DataFrame(rows, columns=["date", "dmi"])


def main():
    oni = fetch_oni()
    dmi = fetch_dmi()
    df = oni.merge(dmi, on="date", how="outer").sort_values("date")
    df = df[(df.date >= "2003-01-01") & (df.date <= "2024-12-01")].reset_index(drop=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False)
    print(f"Saved {len(df)} monthly index records -> {OUT}")
    print(df.dropna().tail())


if __name__ == "__main__":
    main()
