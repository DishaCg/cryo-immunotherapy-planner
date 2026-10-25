# Module A calibration - notes

Before trusting any numbers out of the geometry model, I checked it
against published data. Here's what I used and what I found.

## What I compared against

- Galil IceRod probe (17G, single probe, standard protocol): published
  ice-ball short diameter of 40mm, lethal-zone short diameter of 27mm
  (from Oncohema Key's cryoablation device overview).
- The general rule that cytotoxic temperature (-20 to -40C) shows up
  roughly 3-5mm behind the visible 0C margin (from ClinicalTrials.gov
  protocol docs, NCT01853618 and NCT02821754).
- Boston Scientific's Visual Ice bench guideline: single-probe ice
  balls between 3.7 and 5.3cm max diameter (PMC10289866).

## What I tried, in order

My first version used a point-source probe and an explicit
finite-difference solver, and it was badly wrong - the ice ball
basically stalled after a couple of minutes instead of growing. Two
things were going on: real cryoprobes have an active freezing zone
about 15-20mm long, not a single point, so a point source massively
underestimates how much heat they can pull out. And the explicit
solver couldn't handle the stiff latent-heat term near the freezing
front without a wildly small time step, so on a normal grid it just
got stuck.

Fixing the probe geometry (giving it a proper 17mm active length) and
switching to an implicit finite-volume solver got growth behaving
correctly - no more stalling, and the size roughly doubled. Still
short of the benchmark though.

Two more things helped close the gap:

- Ice conductivity actually goes up as it gets colder (it's not
  constant), so I made frozen-tissue conductivity temperature
  dependent instead of fixed. Picked up another ~20% or so.
- I tested a double freeze-thaw-freeze cycle since that's the real
  clinical protocol, expecting it to close a lot of the remaining gap.
  It barely mattered - final radius came out almost the same as a
  single freeze. Worth recording as a negative result even though it
  didn't help.
- I also tried turning perfusion off entirely, on the theory that
  manufacturer specs are usually measured on gel phantoms with no
  blood flow, not real perfused tissue. That made a bigger difference
  than the double-freeze test - gel-mode prediction landed at 30.6mm
  vs. the 40mm benchmark, noticeably closer than the perfused in-vivo
  number (25.6mm).

## Where it stands right now (SUPERSEDED - see the Oct 2026 addendum below)

> The 40mm / 27mm benchmark and the "~24% / ~51% low" figures in this section came from a secondary source and are not in Boston Scientific's isotherm guide. Use the addendum.

| | model (gel-phantom mode) | published | gap |
|---|---|---|---|
| ice-ball diameter, 10 min, single probe | 30.6mm | 40mm | ~24% low |
| lethal-zone diameter | 13.2mm | 27mm | ~51% low |

So: right order of magnitude, right shape of growth curve, lethal zone
correctly nested inside the visible ice ball - but not dialed in
exactly, especially the lethal zone. My best guess is that either my
-40C cutoff for "lethal" is stricter than whatever criterion the
manufacturer numbers are actually built on, or the generic tissue
properties I'm using are off for whatever tissue their bench test
used. Actually pinning that down would need real gel-phantom or
clinical ice-ball measurements to fit against, which I don't have
access to right now.

## Addendum (Oct 2026): benchmark protocol, grid convergence, and what the data say

**Benchmark protocol.** Boston Scientific's isotherm guide (PI-814202-AB) states its
isotherms were measured in 37C gel after two 10-minute freezes separated by a 5-minute
passive thaw. My earlier comparison used a single 10-minute freeze.

**Published values.** The guide's isotherm table is jumbled in the PDF text layer but its
structure is consistent (each family lists 4, 3, 2, 1-needle panels in decreasing size).
Read that way, the IceRod 1.5 single-needle side-view numbers are about
width 34-36 / 25-26 / 16 mm and length 53-55 / 45 / 39-40 mm at 0 / -20 / -40C
(tolerance +-3 mm width, +-4 mm length). **Confirmed** against Boston Scientific's CX Needle Isotherm Guide (PI-754604) and the Cryoablation Product Catalogue (PI-929802): IceRod 1.5 CX single needle 36x53 / 26x45 / 16x39 mm; IceRod 1.5 PLUS 25x45 / 16x40 mm at -20 / -40C. Independent gel data (Shah et al., Urology 91:234-240, 2016; manufacturer-funded): 17G single-needle ice ball 3.4-3.6 cm, -40C core ~1.8 cm. The 40/27 mm figures
used earlier (Oncohema Key) are not in this guide; do not rely on the old "24% / 51% low".

**Grid convergence (single 600 s freeze, gel, widths 0/-20/-40C; constant 0.8 mm probe radius):**
dr=0.8: 29.5/19.3/12.3; dr=0.4: 31.2/21.1/13.8; dr=0.2: 32.0/21.8/14.5 mm.
Roughly first-order; extrapolated limit ~32.5/22.4/15.1 mm. So dr=1.0 (the old default)
under-predicts widths by ~6-13%; dr=0.4 by ~4-9%. Time step is converged (<1% from dt=4 to 0.5 s).
The optimizer default is now dr=0.4 mm (pass dr_mm=1.0 for quick runs).

**Benchmark protocol at dr=0.4 mm (gel, 17 mm active length, -140C):**
widths 37.0 / 23.7 / 15.0 mm; lengths 42.4 / 30.6 / 24.2 mm.
dr=0.2 mm: widths 37.9 / 24.6 / 15.8; lengths 43.3 / 31.4 / 24.8.
=> widths agree with published within ~3-9% (0C slightly high, -20/-40C slightly low),
but axial lengths are 21-40% too short: the model's ice ball is too round
(length/width 1.15 vs ~1.5 published).

**Active-length sensitivity (dr=1.0, benchmark protocol):** L=17: widths 36.9/23.4/14.7,
lengths 42.9/30.8/24.5. L=30: widths 41.7/27.8/17.7, lengths 55.2/43.3/36.8.
Longer active zones fix the lengths but overshoot the widths (0C by ~16-23%), so no single
active length matches both. A uniform -140C line source is probably the wrong boundary
condition for the aspect ratio; a probe model with finite cooling power / non-uniform tip
temperature is the obvious next step.

Solver assembly is vectorised; it matches the earlier loop version to ~1e-13 and is ~10-30x faster.
