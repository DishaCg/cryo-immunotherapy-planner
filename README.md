# Cryo-Immunotherapy Joint Planning Prototype

This is the working prototype behind my proposal *Joint Prediction of
Ablation Geometry and Immunogenic Response for Personalized
Cryo-Immunotherapy Planning*.

## The idea

Cryoablation planning software today only optimizes for one thing:
does the ice ball cover the tumor. That's been true for a while and
it's a solved-enough problem. But cryoablation is increasingly used
alongside checkpoint-inhibitor immunotherapy, and the freeze protocol
you choose (how fast you cool, how cold you get, how many freeze-thaw
cycles) also changes how strong the resulting immune response is. As
far as I could find, nothing currently plans for that second part.
This project does both at once.

## Layout

```
core/
  geometry_model.py   Bioheat solver. Predicts ice-ball (0C) and
                       lethal-zone (-40C) radius given probe placement
                       and freeze protocol.
  immune_model.py      Immune-activation scoring function, built from
                       literature rather than trained on data (more on
                       why below).
  optimizer.py          Ties the two together and searches candidate
                       protocols for the Pareto-efficient set.
tests/
  test_geometry_model.py
CALIBRATION_NOTES.md  How I checked Module A against published numbers.
```

## Where each piece actually stands

**Geometry model.** This solves the Pennes bioheat equation on an
axisymmetric grid, implicit finite-volume, with latent heat handled
through an enthalpy method. It passes a small regression suite (the
ice ball has to grow monotonically, the lethal isotherm has to stay
inside the visible one, two probes have to do better than one) and
it's in the right ballpark against published Galil IceRod numbers --
not exact, and I've written up exactly where the gap is and my best
guess at why in CALIBRATION_NOTES.md.

**Immune model.** I initially wanted to train something on real
outcome data, but that data doesn't really exist publicly at the
granularity you'd need -- studies report trial-level response rates,
not per-case freeze-parameter-to-immune-marker pairs. So instead this
is a scoring function where every factor's weight comes from a
specific paper or mechanism, and the output always comes with an
evidence-quality tag rather than pretending to be more certain than
it is.

**Optimizer.** Runs both modules across a set of candidate protocols
and returns the Pareto front instead of mashing the two scores
together into one fake "best" number -- a physics prediction and a
literature guess aren't the same kind of evidence and I didn't want to
hide that.

**Not built yet:** any kind of UI, a wider search over multi-probe
configurations, and real calibration data for the geometry model.

## Running it

```bash
pip install numpy scipy
python3 -c "
from core.optimizer import search_protocols
results = search_protocols(tumor_radius_mm=10.0)
for e in results:
    if e.is_pareto_optimal:
        print(e.protocol, '-> lethal radius:', round(e.lethal_radius_mm, 1),
              'mm, immune score:', e.immune_score)
"
```

## One thing to know if you're poking at this

Each protocol evaluation takes 10-15 seconds because of the solver, so
a full grid search isn't fast. Fine for testing a handful of protocols
offline, not fine for an interactive tool yet -- that'd need a faster
or cached solver mode first.
