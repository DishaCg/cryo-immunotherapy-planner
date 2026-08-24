# Module A — Calibration Notes (v1)

This documents a calibration pass of the geometry model against
**published, citable benchmark data**, done before any external
outreach — per the reasoning that a model with no connection to real
data isn't ready to show a potential collaborator.

## Benchmark data used (all public, no patient data required)

- Galil IceRod cryoprobe (17G, single probe, standard clinical
  protocol): published ice-ball short diameter = 40mm, corresponding
  lethal-zone short diameter = 27mm.
  Source: Oncohema Key, *Cryoablation: Mechanism of Action and Devices*.
- General principle: cytotoxic temperatures (−20 to −40°C) occur
  3–5mm behind the visualized 0°C ice margin.
  Source: ClinicalTrials.gov protocol documents (NCT01853618, NCT02821754).
- Manufacturer bench guideline: single-probe ice-ball maximum
  dimension range 3.7–5.3cm.
  Source: PMC10289866 (Visual Ice / Boston Scientific system study).

## Calibration steps taken, in order, with reasoning

| # | Change | Why | Effect |
|---|--------|-----|--------|
| 1 | Elongated probe active-freeze zone (17mm) instead of point sink | Commercial cryoprobes have an active freezing section, not a point tip; point-sink drastically under-predicts total heat extraction | Iceball radius roughly doubled |
| 2 | Switched explicit FDM → implicit finite-volume solver | Explicit scheme was numerically stiff at the latent-heat front and stalled growth (~3x under-prediction) | Growth curve became physically sensible (no stalling) |
| 3 | Temperature-dependent frozen thermal conductivity (ice conductivity genuinely increases at colder sub-zero temps, not constant) | Physical property, not a tuning knob — supported by ice thermal-conductivity literature (e.g. Slack 1980) | Closed roughly another 20% of the gap |
| 4 | Tested double freeze-thaw-freeze protocol (matches standard clinical practice) vs. single freeze | Manufacturer specs are usually quoted for the standard two-cycle protocol | Negligible effect in this model — useful negative finding, not the source of the gap |
| 5 | Tested no-perfusion ("gel phantom") vs. in-vivo perfused tissue | Manufacturer bench specs are typically measured in gel/agar phantoms with no blood flow, not perfused tissue | Gel-phantom prediction (30.6mm) tracks much closer to the 40mm bench spec than the perfused in-vivo prediction (25.6mm) — consistent with the literature distinction between bench and in-vivo ice-ball size |

## Current state (v1)

| Quantity | Model (gel-phantom mode) | Published benchmark | Gap |
|---|---|---|---|
| Ice-ball diameter, single probe, 10 min | 30.6mm | 40mm (IceRod) | ~24% low |
| Lethal-zone diameter | 13.2mm | 27mm (IceRod) | ~51% low |

**Order of magnitude: correct. Absolute calibration: not yet exact.**
The remaining gap is most likely in: (a) the −40°C lethal-isotherm
threshold possibly being stricter than the criterion implicit in the
manufacturer's "lethal zone" figure, and (b) generic (non-probe-
specific) tissue property assumptions. This is exactly the kind of
gap that real cryoablation modeling papers close using their own
bench/gel experiments — which is a legitimate, honest thing to name
explicitly as a next step rather than something to hide.

## Why this matters for outreach

This is no longer "an idea with no connection to reality." It's a
working model that has been checked against real published numbers,
found to be directionally and qualitatively correct, and has an
identified, specific, well-reasoned remaining calibration gap. That
last part is actually useful in an email to a potential advisor — it
gives them a concrete, bounded, collaborative next step (e.g., "would
your lab have access to gel-phantom or clinical ice-ball measurements
I could calibrate against?") rather than an open-ended ask.
