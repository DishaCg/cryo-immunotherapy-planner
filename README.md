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
analysis/
  sensitivity_analysis.py  Joint-optimisation example and weight / ablation-
                       fraction sensitivity analysis used in the paper.
  convergence_study.py Grid / time-step convergence and benchmark comparison.
CALIBRATION_NOTES.md  How I checked Module A against published numbers.
```

## Where each piece actually stands

**Geometry model.** This solves the Pennes bioheat equation on an
axisymmetric grid, implicit finite-volume, with latent heat handled
through an enthalpy method. It passes a small regression suite (the
ice ball has to grow monotonically, the lethal isotherm has to stay
inside the visible one, the -20C isotherm sits between them, and the
published gel-phantom numbers stay pinned) and it's in the right ballpark
against Boston Scientific's published IceRod isotherm data: widths agree
within ~3-9% at the benchmark protocol, but axial lengths come out 21-40%
too short (ice ball too round). The full comparison and my best guess at
the cause are in CALIBRATION_NOTES.md.

**Immune model.** I initially wanted to train something on real
outcome data, but that data doesn't really exist publicly at the
granularity you'd need - studies report trial-level response rates,
not per-case freeze-parameter-to-immune-marker pairs. So instead this
is a scoring function where every factor's weight comes from a
specific paper or mechanism, and the output always comes with an
evidence-quality tag rather than pretending to be more certain than
it is.

**Optimizer.** Runs both modules across a set of candidate protocols
and returns the Pareto front instead of mashing the two scores
together into one fake "best" number - a physics prediction and a
literature guess aren't the same kind of evidence and I didn't want to
hide that.

**Not built yet:** any kind of UI, multi-probe configurations (the solver is
axisymmetric, so an off-axis probe is a ring source, not a second probe - this
needs a 3-D discretisation), and real calibration data for the geometry model.

## Running it

```bash
pip install numpy scipy
python3 -c "
from core.optimizer import search_protocols
results = search_protocols(tumor_radius_mm=6.8,
                           freeze_time_options_s=(300.0, 450.0, 600.0, 900.0),
                           checkpoint_options=(True,))
for e in results:
    if e.is_pareto_optimal:
        print(e.protocol, '-> lethal radius:', round(e.lethal_radius_mm, 1),
              'mm, immune score:', e.immune_score)
"
```

To reproduce the paper's joint-optimisation example and sensitivity analysis:

```bash
python -m analysis.sensitivity_analysis
```

## One thing to know if you're poking at this

The solver's matrix assembly is vectorised. At the default grid spacing
(dr_mm=0.4, which is what the paper uses) a protocol evaluation takes roughly 5-40 seconds
(two-cycle protocols are the slowest); pass dr_mm=1.0 for quick, coarser runs (about
1-3 seconds, but widths are ~5-10% smaller - see CALIBRATION_NOTES.md). A handful of
protocols is fine; a large grid is still not interactive-speed.

The immune score's only interaction with the geometry is the ablation-fraction
factor, and the Pareto trade-off in the paper's example comes entirely from its
contested "antigen depot" penalty above 95% coverage (see the paper, Section V-D).
