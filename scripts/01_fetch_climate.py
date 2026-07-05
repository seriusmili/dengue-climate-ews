"""
01_fetch_climate.py
-------------------
Download REAL monthly climate for each analysed province from the Open-Meteo
Historical Weather API (ERA5-based). Free, NO API key required.
Docs: https://open-meteo.com/en/docs/historical-weather-api

Overwrites data/processed/climate_provinces_monthly.csv with real data.
Usage:  python scripts/01_fetch_climate.py

Note: the free API rate-limits bursts of requests (HTTP 429). This script
therefore waits between provinces and retries automatically with backoff.
"""
import pathlib, time
import requests
import pandas as pd
from provinces import PROVINCES, START, END

DAILY = ["temperature_2m_mean", "temperature_2m_max", "temperature_2m_min",
         "precipitation_sum", "relative_humidity_2m_mean"]
OUT = pathlib.Path(__file__).resolve().parents[1] / "data" / "processed" / "climate_provinces_monthly.csv"

PAUSE = 12          # seconds to wait between provinces (be gentle to the free API)
MAX_RETRIES = 6     # how many times to retry one province if rate-limited


def fetch(lat, lon):
    """Fetch one point, retrying with growing waits if the API says 429."""
    url = "https://archive-api.open-meteo.com/v1/archive"
    params = {"latitude": lat, "longitude": lon, "start_date": START, "end_date": END,
              "daily": ",".join(DAILY), "timezone": "Asia/Jakarta"}
    for attempt in range(1, MAX_RETRIES + 1):
        r = requests.get(url, timeout=60, params=params)
        if r.status_code == 429:
            wait = 20 * attempt   # 20s, 40s, 60s, ... progressively longer
            print(f"   ...rate-limited (429). Waiting {wait}s then retrying "
                  f"(attempt {attempt}/{MAX_RETRIES})")
            time.sleep(wait)
            continue
        r.raise_for_status()
        d = pd.DataFrame(r.json()["daily"]); d["time"] = pd.to_datetime(d["time"])
        return d
    raise RuntimeError("Still rate-limited after several retries. "
                       "Please wait a few minutes and run this script again.")


def main():
    frames = []
    provs = list(PROVINCES.items())
    for i, (prov, (lat, lon)) in enumerate(provs, start=1):
        print(f"[{i}/{len(provs)}] Fetching {prov} ({lat},{lon}) ...")
        d = fetch(lat, lon)
        d["date"] = d["time"].dt.to_period("M").dt.to_timestamp()
        m = d.groupby("date").agg(
            temp_mean=("temperature_2m_mean", "mean"),
            temp_max=("temperature_2m_max", "mean"),
            temp_min=("temperature_2m_min", "mean"),
            precip_total=("precipitation_sum", "sum"),
            rh_mean=("relative_humidity_2m_mean", "mean")).reset_index()
        m.insert(0, "province", prov)
        frames.append(m)
        if i < len(provs):
            print(f"   done. Pausing {PAUSE}s before the next province...")
            time.sleep(PAUSE)
    out = pd.concat(frames, ignore_index=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)
    print(f"\nSaved {len(out)} rows for {out.province.nunique()} provinces -> {OUT}")
    print("All 7 provinces fetched successfully. You can now run 01b, then 02-05.")


if __name__ == "__main__":
    main()
