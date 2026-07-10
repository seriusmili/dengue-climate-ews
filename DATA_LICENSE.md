# Data provenance, attribution, and licensing

This repository combines material with **different origins and different terms**.
Please read this before reusing anything in `data/`.

| Material | Terms |
|----------|-------|
| **Code** — `scripts/`, `analysis/`, `run_all.sh` | MIT (see `LICENSE`) |
| **Dengue data** — `data/raw/dengue_*.csv`, and the OpenDengue-derived columns in `data/processed/analysis_*.csv` | Redistributed under the terms set by **OpenDengue**. See below. |
| **Population table** — `data/raw/population_provinces.csv` | Compiled by the author from BPS Statistics Indonesia figures |
| **Climate data** — `data/processed/climate_*.csv` | Open-Meteo **CC BY 4.0** (over Copernicus/ERA5 **CC-BY**); NOAA indices are US Government works, not copyrighted. See below. |
| **Figures, preprint, and written analysis** | © the author; reuse with attribution |

---

## Dengue data: OpenDengue

Source: **OpenDengue database, version 1.3**, LSHTM Dengue Mapping and Modelling Group.

### A note on the licence

OpenDengue states its licence in two places, and **the statements differ**:

| Where | Stated licence |
|-------|----------------|
| Peer-reviewed data paper (Sci Data, 2024) | Creative Commons **CC BY-SA** |
| Project website footer, including the data page | **CC BY 4.0** |

This repository **does not attempt to resolve that ambiguity, and does not
relicense OpenDengue's data.** The extracts here are redistributed under whatever
terms OpenDengue in fact applies. Both candidate licences expressly permit
redistribution and adaptation with attribution, so republishing these extracts is
permitted under either reading; they differ only in whether adaptations must
carry the same licence (ShareAlike).

**If you intend to redistribute or relicense these files yourself, consult the
OpenDengue project directly** rather than relying on this repository.

### Attribution

> Clarke J, Lim A, Gupte P, Pigott DM, van Panhuis WG, Brady OJ.
> A global dataset of publicly available dengue case count data.
> *Scientific Data*. 2024;11:296. https://doi.org/10.1038/s41597-024-03120-7

> Clarke J, Lim A, Gupte P, Pigott DM, van Panhuis WG, Brady OJ.
> OpenDengue: data from the OpenDengue database. Version 1.3. figshare; 2025.
> https://doi.org/10.6084/m9.figshare.24259573

OpenDengue collates data from publicly available sources, including official
ministry of health reporting. Original sources for individual records can be
traced via the record UUID in OpenDengue's source data file.

### Changes made to the original data

Both CC BY and CC BY-SA require that modifications be indicated. The changes here:

- Filtered to Indonesia, then to seven provinces (East Java, West Java, Central
  Java, DKI Jakarta, Banten, Bali, DI Yogyakarta).
- Filtered to monthly temporal resolution (`T_res == "Month"`).
- Yogyakarta name variants merged into a single label.
- Records aggregated by province and month; `dengue_total` summed.
- In `data/processed/`, joined with population and climate data and augmented with
  derived columns (incidence per 100,000, month, time index).

No case counts were altered. The transformation is fully reproducible from
`scripts/02_build_dataset.py`.

---

## Climate data

### Temperature and rainfall — Open-Meteo (ERA5 reanalysis)

`data/processed/climate_provinces_monthly.csv`, retrieved by
`scripts/01_fetch_climate.py` from the Open-Meteo Historical Weather API.

**Two licences apply in layers:**

1. **Open-Meteo API data** is offered under **CC BY 4.0**. Attribution is required,
   including a link to Open-Meteo wherever its data is displayed. Open-Meteo's free
   API tier is for **non-commercial use**; public research at public institutions is
   listed by Open-Meteo as qualifying non-commercial use.

2. **The underlying reanalysis is ERA5**, produced by the **Copernicus Climate Change
   Service (C3S)** at ECMWF. Since 2 July 2025 the Climate Data Store distributes
   ERA5 under a **CC-BY licence**. Copernicus requires a specific attribution notice,
   and because this repository aggregates the data (daily → monthly, averaged across
   five points per province), the *modified* form of the notice applies.

**Required notices** — reproduced here to satisfy both licences:

> Weather data by [Open-Meteo.com](https://open-meteo.com/), licensed under
> [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

> Contains modified Copernicus Climate Change Service information 2026.
> Neither the European Commission nor ECMWF is responsible for any use that may be
> made of the Copernicus information or data it contains.

Suggested citations:

- Zippenfenig, P. (2023). *Open-Meteo.com Weather API* [Computer software]. Zenodo.
- Hersbach, H. et al. (2023). *ERA5 hourly data on single levels from 1940 to present.*
  Copernicus Climate Change Service (C3S) Climate Data Store (CDS).
  DOI: 10.24381/cds.adbb2d47

**Changes made:** daily values retrieved for five representative coordinates per
province; averaged across points; aggregated to monthly (means for temperature and
humidity, sums for precipitation). No values were otherwise altered. Fully
reproducible from `scripts/01_fetch_climate.py`.

### ENSO and IOD indices — NOAA

`data/processed/climate_indices_monthly.csv`, retrieved by
`scripts/01b_fetch_climate_indices.py`:

- **ONI** (Oceanic Niño Index) — NOAA Climate Prediction Center
- **DMI** (Dipole Mode Index) — NOAA Physical Sciences Laboratory

Information created by the U.S. Government and presented on U.S. Government websites
is **not subject to copyright in the United States** and may be used freely for
lawful purposes. Attribution is nonetheless expected as a scholarly norm, and NOAA
requests acknowledgement as the source. Note also that NOAA's name and emblem are
protected by trademark and must not be used in a way implying endorsement.

**Changes made:** ONI season codes mapped to their centre month; DMI parsed from
fixed-width annual rows and missing-value flags dropped; both restricted to
2003–2024 and merged into one monthly table. No values were altered.

---

## Disclaimer

This document records provenance and attribution in good faith. It is not legal
advice, and it is not a licence grant over material the author does not own. For
authoritative terms, consult each data provider.
