# Two analysis tracks

This pipeline supports two distinct questions that should **not** share a model.
They are kept separate on purpose: the spatial smooth that makes Track A work is
the same term that would bias Track B, and the longitudinal structure Track B
needs is irrelevant to Track A. They share only the predictor-extraction front end.

| | **Track A — interpolation** | **Track B — degradation / re-occupation** |
|---|---|---|
| Question | What is orangutan abundance across the whole landscape in a given year? | Does logging/degradation drive a temporary change in abundance, and do orangutans return? |
| Unit | grid cell | transect **segment** × year |
| Model | NB GAM with `s(x_km, y_km)` spatial smooth | NB GAM, recovery smooth `s(years_since_last_logging, by = ever_logged)`, transect + year random effects |
| Role of `s(x,y)` | **central** — interpolates across unsurveyed space | **absent** (or at most a nuisance mop-up if residuals are spatially clustered) |
| Output | wall-to-wall density surface, summed to a landscape total | degradation contrast + recovery trajectory in occupied habitat |
| Scripts | `prepare_predictors_*`, `best_model/`, `*_predict_*` | `prepare_segments_track_b.R`, `fit_track_b.R` |

## Why they must stay separate

Track A's `s(x_km, y_km)` smooth exists to borrow strength across space and fill
gaps. That same smooth, in Track B, would absorb the spatially structured
environmental signal (logged vs unlogged areas) that *is* the thing we want to
measure. Conversely, Track B's transect random effect and recovery smooth answer
a change-over-time question that Track A's per-year prediction surface cannot.

---

# Track B: design, model, and what it can honestly claim

## Survey design (the constraints that shape everything)

- Baseline survey on a systematic grid of 1 km transects in **2017–2018**.
- Resurveys in **2022, 2023, 2024, 2025**, but **only transects that had
  orangutans at baseline were revisited.** Revisit selection was
  **independent of disturbance.**
- Nest GPS points carry a transect id; transects are split into fixed-length
  **segments** (≈250 m) so within-transect degradation gradients are preserved.
- Per-segment covariates are extracted from a **tight buffer** (buffer ≤ segment
  length) so adjacent segments do not share most of their covariate values.

## Model

```
nr_nests ~ ever_logged                                  # logged vs never-logged baseline
         + s(years_since_last_logging, by = ever_logged) # recovery curve, LOGGED segments only
         + prop_degraded_1yr_prior                       # recent/standing degradation intensity
         + baseline_nests_c                              # partial regression-to-the-mean control
         + dem + slope + distance_road                   # environmental covariates
         + s(transect_id, bs = "re")                     # repeated visits to same transect
         + s(year_f,      bs = "re")                     # survey-wide year variation
         + offset(offset_term)                           # nest -> density conversion (log scale)
```

Key design choices:

- **Negative binomial**, counts per segment, zeros retained.
- **`ever_logged` is an ordered factor** so the `by =` recovery smooth is zero on
  never-logged segments; those act as the flat reference baseline and the curve is
  estimated only from logged segments.
- **`years_since_last_logging`** is the re-occupation axis. Because the *forest*
  is not expected to recover structurally within the study window, this term
  measures behavioural re-occupation of degraded-but-habitable forest, not
  habitat regrowth. The hypothesis: abundance dips just after logging, then climbs
  back toward the never-logged baseline.
- **`prop_degraded_1yr_prior`** measures standing/acute pressure. Because
  orangutans tolerate degraded forest, this coefficient may be weak; the signal is
  expected in the *timing* term, not the *amount* term.
- **`offset_term`** carries the nest→density constants (effective strip width,
  proportion of builders, nest production rate, nest decay time) on the log scale.

## What Track B can and cannot claim

**Can claim (defensible):**
- The **degradation/recovery contrast**. Because revisit selection was
  independent of disturbance, the logged-vs-intact comparison *within* the
  revisited set is not confounded by the selection mechanism. This is the
  headline result.

**Cannot claim without heavy caveats (selection effects):**
- **Absolute trajectory level.** Transects were selected on high baseline
  occupancy, so the set as a whole drifts downward on remeasurement by regression
  to the mean, independent of any real change. `baseline_nests_c` partially
  corrects this; the `fit_track_b.R` sensitivity check (fit with vs without it)
  quantifies how much. Report the with-baseline version and state the limitation.
- **Recolonisation / recovery beyond baseline.** Transects empty at baseline were
  dropped, so colonisation of previously unoccupied habitat is **unobservable**.
  Recovery is therefore censored at the baseline-occupied level: we can see a
  segment dip and return toward its own baseline, but not orangutans appearing in
  forest that was empty before.

## Recommended framing line

> We characterise the trajectory of abundance in **known-occupied** habitat, and
> test whether degradation timing is associated with steeper declines and, where
> conditions allow, partial re-occupation — while acknowledging that selection on
> baseline occupancy censors recolonisation and biases absolute trends toward
> apparent decline.

## Diagnostics built into `fit_track_b.R`

- `gam.check` (is the recovery smooth's `k` high enough?).
- **Recovery curve plot** — the headline visual; dip-then-return = re-occupation.
- **Between-transect Moran's I** on residuals — if strongly positive, add
  `s(x_km, y_km)` as a *nuisance* term and refit (this is the only legitimate use
  of a spatial smooth in Track B, and it is not the Track A interpolation use).
- **Regression-to-the-mean sensitivity** — refit without `baseline_nests_c`.
- **distance_road sensitivity** — refit without it; if degradation terms shrink
  when road is included, roads (e.g. logging roads) may be a *mediator* of the
  degradation effect rather than a confounder, so report both specifications.

## Caveats to watch in the data

- **Nest decay time.** If decay differs between degraded and intact forest, or
  between wetter/drier years, a constant decay term in the offset can manufacture
  the very trend Track B tests for. Scrutinise this parameter first.
- **Buffer size.** Keep the buffer ≤ segment length. Wide buffers (e.g. 3 km)
  around segments make neighbours share covariate values, manufacturing spatial
  autocorrelation and eroding the within-transect contrast.
