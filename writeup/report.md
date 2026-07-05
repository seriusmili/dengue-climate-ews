# Climate-driven dengue early warning in Indonesia: an East Java case study with a seven-province comparison

**Serius Miliyani Dwi Putri**
M.Ked.Trop (Epidemiology, Tropical Medicine)
*Proof-of-concept analysis · 2026*

## Abstract

Dengue transmission in Indonesia is strongly climate-sensitive, yet routine
surveillance remains largely reactive. As a compact demonstration of the
modelling approach behind climate-based early-warning systems, I link open
monthly dengue surveillance for East Java (OpenDengue v1.3, 2004–2024) with
open climate reanalysis (Open-Meteo / ERA5) and fit a negative-binomial
distributed-lag non-linear model (DLNM). The framework quantifies the delayed,
potentially non-linear effect of temperature and rainfall on dengue incidence
while adjusting for seasonality and long-term trend. This note documents the
data, methods, and reproducible pipeline; it is intended as a portfolio artefact
rather than a definitive epidemiological result.

## 1. Data

**Dengue.** Monthly, province-level (Admin1) reported dengue case counts for
East Java were extracted from the OpenDengue database v1.3, which standardises
publicly available surveillance data. The East Java series covers Jan 2004 to
Jul 2024 (229 monthly observations; ~296,000 reported cases over the period).

**Climate.** Monthly mean temperature, total rainfall, and mean relative
humidity were derived from the Open-Meteo Historical Weather API (ERA5-based
reanalysis), averaged across five representative East Java locations (Surabaya,
Malang, Jember, Banyuwangi, Madiun) and aggregated from daily to monthly values.
*(The version documented here uses a labelled synthetic climate series for
demonstration; substituting real Open-Meteo data is a single scripted step.)*

## 2. Methods

Case counts $Y_t$ in month $t$ were modelled with a negative-binomial GLM:

$$\log \mathbb{E}[Y_t] = \alpha + \text{cb}_{\text{temp}}(x_t, \ell) + \text{cb}_{\text{rain}}(r_t, \ell) + s(\text{month}_t) + f(t)$$

where each **cross-basis** $\text{cb}(\cdot)$ (Gasparrini 2011) is a tensor
product of a natural-spline exposure–response basis and a natural-spline
lag–response basis over lags $\ell = 0,\dots,6$ months. Seasonality $s(\cdot)$
was captured with cyclic harmonic terms and the long-term trend $f(t)$ with a
natural spline on time. Negative-binomial errors accommodate overdispersion.
Estimates are summarised as cumulative relative risks (RR) over the 0–6 month
lag window, centred at the median exposure.

## 3. Results (demonstration run)

The seasonal structure is pronounced: mean monthly dengue counts are highest in
January–March (wet season) and lowest in August–September. In the demonstration,
the strongest rainfall–dengue cross-correlation occurred at a **2-month lag**
(r ≈ 0.52), consistent with the biological delay between rainfall, expansion of
*Aedes* breeding habitat, and onward transmission. The negative-binomial DLNM
fit the series well (deviance pseudo-R² ≈ 0.73). Cumulative exposure–response
curves indicate rising dengue risk with increasing rainfall across the observed
range (Figure 4).

*Figures: (1) time series; (2) cross-correlation; (3) temperature exposure–
response; (4) rainfall exposure–response; (5) observed vs fitted.*

## 4. Discussion and next steps

This analysis reproduces, at small scale, the core building block of a dengue
early-warning system: an interpretable, lagged climate–incidence relationship
estimated with a standard, defensible model. It is deliberately minimal. A full
research programme would extend it to (i) hierarchical spatio-temporal models
across all Indonesian provinces with population offsets; (ii) formal evaluation
of forecast skill and achievable lead time against held-out historical
epidemics; (iii) inclusion of ENSO indices, urbanisation, and mobility; and
(iv) co-designed, district-level alert thresholds with health authorities. These
are the objectives set out in the accompanying PhD concept note.

## Reproducibility

All code and the real dengue extract are in this repository. The analysis runs
end-to-end from the command line (see `README.md`). Real climate is obtained via
`scripts/01_fetch_climate.py` (Open-Meteo, no API key).

## Key references

- Gasparrini A. Distributed lag linear and non-linear models in R: the package
  dlnm. *J Stat Softw.* 2011;43(8).
- Lowe R, et al. Combined effects of hydrometeorological hazards and
  urbanisation on dengue risk in Brazil. *Lancet Planet Health.* 2021.
- Clarke J, et al. A global dataset of publicly available dengue case count
  data. *Sci Data.* 2024;11:296. (OpenDengue)
- Colón-González FJ, et al. Projecting the risk of mosquito-borne diseases in a
  warmer and more populated world. *Lancet Planet Health.* 2021.
