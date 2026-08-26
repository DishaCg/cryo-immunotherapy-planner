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

## Where it stands right now

| | model (gel-phantom mode) | published | gap |
|---|---|---|---|
| ice-ball diameter, 10 min, single probe | 30.6mm | 40mm | ~24% low |
| lethal-zone diameter | 13.2mm | 27mm | ~51% low |

So: right order of magnitude, right shape of growth curve, lethal zone
correctly nested inside the visible ice ball -- but not dialed in
exactly, especially the lethal zone. My best guess is that either my
-40C cutoff for "lethal" is stricter than whatever criterion the
manufacturer numbers are actually built on, or the generic tissue
properties I'm using are off for whatever tissue their bench test
used. Actually pinning that down would need real gel-phantom or
clinical ice-ball measurements to fit against, which I don't have
access to right now.
