# Cryo-Immunotherapy Joint Planning Prototype

Research prototype accompanying: *Joint Prediction of Ablation Geometry
and Immunogenic Response for Personalized Cryo-Immunotherapy Planning*
(see `Cryo_Immunotherapy_Planning_Proposal.docx`).

## What this is

Existing cryoablation planning tools optimize freeze protocols for one
objective: geometric tumor coverage. This prototype optimizes for two,
jointly: geometric coverage AND predicted immune-activation potential
-- because a growing share of clinical cryoablation is now delivered
in combination with checkpoint-inhibitor immunotherapy, and no
existing tool accounts for that second objective when recommending a
freeze protocol.

## Structure

```
core/
  geometry_model.py   Module A -- implicit finite-volume bioheat solver.
                       Predicts ice-ball (0C) and lethal-zone (-40C)
                       radius from probe placement and freeze protocol.
  immune_model.py      Module B -- literature-parameterized immune-
                       activation scoring function. Explicitly NOT a
                       trained ML classifier (see file docstring for
                       why -- no adequate public training dataset
                       currently exists).
  optimizer.py         Integration layer -- searches candidate freeze
                       protocols, evaluates both modules, returns the
                       Pareto-efficient set (protocols where you can't
                       improve one objective without giving up the
                       other).
tests/
  test_geometry_model.py   Sanity + regression tests for Module A.
CALIBRATION_NOTES.md  Honest record of Module A's validation against
                       published benchmark data (Galil IceRod specs),
                       what was tried, what worked, and the documented
                       remaining gap.
```

## Status (be upfront about this with anyone you show it to)

- **Module A**: numerically validated (passes regression tests),
  order-of-magnitude checked against published cryoprobe benchmark
  data, with a documented ~20-50% remaining calibration gap. Not yet
  fitted to gel-phantom or clinical ground truth -- that's the
  natural next step, and a good, concrete thing to ask a potential
  advisor for help accessing.
- **Module B**: an interpretable, literature-cited scoring heuristic,
  not a trained model, because the outcome-level dataset needed to
  train one doesn't exist in usable public form. This is stated
  explicitly in the code and in its output, not hidden.
- **Integration layer**: functionally complete and tested end-to-end;
  correctly identifies Pareto-optimal protocols rather than collapsing
  two different kinds of evidence (physics-grounded vs. literature
  heuristic) into one falsely-precise number.
- **Not done yet**: an actual demo interface (planned next), broader
  multi-probe configurations, and real calibration data.

## Running it

```bash
pip install numpy scipy
python3 -c "
from core.optimizer import search_protocols
results = search_protocols(tumor_radius_mm=10.0)
for e in results:
    if e.is_pareto_optimal:
        print(e.protocol, '-> lethal radius:', round(e.lethal_radius_mm,1),
              'mm, immune score:', e.immune_score)
"
```

## Known limitation to flag proactively

The current grid search in `optimizer.py` takes roughly 10-15s per
protocol evaluation (finite-volume solve is the bottleneck). Fine for
a research prototype and a handful of candidate protocols; would need
a coarser/faster solver mode or caching before it could back an
interactive demo UI. Worth mentioning as a known next step, not a
hidden flaw.
