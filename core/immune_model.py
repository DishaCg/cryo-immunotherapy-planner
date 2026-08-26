"""
Module B: immune-response / immunogenicity potential estimator.

I originally wanted this to be a trained model, but there's no public
dataset that actually links per-case freeze parameters to measured
immune outcomes -- published cryo-immunotherapy studies give you
trial-level response rates and mechanism descriptions, not the kind
of per-case table you'd need to train on. Fabricating a training set
or fitting to a handful of data points would just produce something
that looks rigorous without being rigorous.

So this is a scoring function instead: every factor's direction and
weight comes from a specific cited mechanism or finding, combined
into one interpretable score, with an evidence-quality label attached
to the output so it doesn't imply more certainty than it has. Turning
this into an actual trained model later would need a real dataset
linking specific protocols to measured immune markers -- gel-phantom
or clinical trial data would do it.
"""

from dataclasses import dataclass, field
from typing import Literal


@dataclass
class FreezeProtocol:
    cooling_rate_c_per_min: float   # e.g. -25 to -50 C/min typical for cryoprobes
    min_temp_c: float               # nadir temperature reached, e.g. -140
    num_freeze_thaw_cycles: int     # typically 1-2 clinically
    ablation_fraction: float        # fraction of tumor volume within lethal isotherm, 0-1
    combined_with_checkpoint_inhibitor: bool = False
    tumor_type: Literal["breast", "lung", "renal", "hepatic", "other"] = "other"


@dataclass
class ImmuneEstimate:
    activation_score: float          # 0-100, composite, higher = more favorable
    contributing_factors: dict       # factor_name -> (raw_value, weighted_contribution, rationale)
    evidence_quality: str            # explicit evidence-basis label
    caveats: list


# ----------------------------------------------------------------------
# Literature-derived factor functions
# Each returns a 0-1 "favorability" score and a one-line rationale with
# its source, so every number in the final score is traceable.
# ----------------------------------------------------------------------

def _cooling_rate_factor(rate_c_per_min: float):
    """Faster cooling favors intracellular ice formation and necrotic
    (rather than apoptotic) cell death, which is associated with
    greater release of antigens and damage-associated molecular
    patterns (DAMPs) -- the immunogenic cell death mechanism underlying
    cryo-immunotherapy. Saturates beyond ~40 C/min (typical clinical
    cryoprobes already achieve very fast local cooling near the probe;
    additional speed has diminishing marginal effect).
    Rationale source: cryoablation mechanism reviews (Frontiers in
    Oncology 2026 cryo-immunotherapy review; general immunogenic cell
    death literature).
    """
    rate = abs(rate_c_per_min)
    score = min(rate / 40.0, 1.0)
    return score, f"cooling rate {rate:.0f} C/min -> favors necrotic (immunogenic) death up to ~40 C/min saturation"


def _min_temp_factor(min_temp_c: float):
    """Colder nadir temperature increases the fraction of tissue
    undergoing complete lethal freezing and antigen release, with
    diminishing returns below about -100C (tissue already fully
    committed to the frozen/necrotic pathway).
    """
    depth = max(0.0, -min_temp_c)
    score = min(depth / 100.0, 1.0)
    return score, f"nadir {min_temp_c:.0f} C -> deeper freeze increases complete antigen release, saturating near -100C"


def _freeze_thaw_cycles_factor(n_cycles: int):
    """A second freeze-thaw cycle is associated with more complete
    cell death and antigen release than a single cycle; benefit beyond
    2 cycles is not well established and isn't assumed here.
    """
    score = {0: 0.0, 1: 0.55, 2: 1.0}.get(n_cycles, 1.0 if n_cycles > 2 else 0.0)
    return score, f"{n_cycles} freeze-thaw cycle(s) -> repeat freezing improves antigen release completeness (benefit not established beyond 2 cycles)"


def _ablation_fraction_factor(fraction: float):
    """This factor addresses a question that is genuinely contested in
    the literature. Some evidence favors near-complete ablation (more
    antigen, more damage signal); other work on partial or fractional
    cryoablation argues that leaving some residual tumor preserves an
    antigen depot that sustains immune engagement rather than removing
    the antigen source entirely. The function is intentionally flat
    across a broad mid-range, penalizing only the extremes, to avoid
    asserting a preference the literature does not currently support.
    """
    if fraction < 0.3:
        score = fraction / 0.3 * 0.4
        rationale = "low ablation fraction -> limited antigen release"
    elif fraction > 0.95:
        score = 0.75  # slightly penalized vs mid-range, reflecting the "antigen depot" argument
        rationale = "near-total ablation -> strong local control but debated whether it removes antigen depot needed to sustain immune response"
    else:
        score = 0.85
        rationale = "moderate-to-high ablation fraction -> within the range most studies operate in; this factor is contested in the literature and treated as near-neutral here"
    return score, rationale


def _checkpoint_combo_factor(combined: bool):
    """Combining cryoablation with PD-1/PD-L1 checkpoint inhibition is
    the single most consistently reported factor increasing systemic
    (abscopal) response in the literature surveyed for this project
    (JAK2-STAT3-S100A8/A9 axis mechanism; multiple combination trials).
    Weighted heavily and separately from the freeze-parameter factors
    because it is a treatment-strategy choice, not a freeze physics
    parameter.
    """
    return (1.0 if combined else 0.15), (
        "combined with checkpoint inhibitor -> literature's most consistent driver of systemic response"
        if combined else
        "cryoablation alone -> abscopal responses reported but substantially less consistent without checkpoint combination"
    )


# ----------------------------------------------------------------------
# Composite estimator
# ----------------------------------------------------------------------

# weights sum to 1.0; checkpoint-combination weighted highest per the
# literature signal being the strongest/most consistent of the factors
_WEIGHTS = {
    "cooling_rate": 0.20,
    "min_temp": 0.15,
    "freeze_thaw_cycles": 0.20,
    "ablation_fraction": 0.15,
    "checkpoint_combo": 0.30,
}


def estimate(protocol: FreezeProtocol) -> ImmuneEstimate:
    cr_score, cr_r = _cooling_rate_factor(protocol.cooling_rate_c_per_min)
    mt_score, mt_r = _min_temp_factor(protocol.min_temp_c)
    fc_score, fc_r = _freeze_thaw_cycles_factor(protocol.num_freeze_thaw_cycles)
    af_score, af_r = _ablation_fraction_factor(protocol.ablation_fraction)
    ck_score, ck_r = _checkpoint_combo_factor(protocol.combined_with_checkpoint_inhibitor)

    contributions = {
        "cooling_rate": (protocol.cooling_rate_c_per_min, cr_score * _WEIGHTS["cooling_rate"], cr_r),
        "min_temp": (protocol.min_temp_c, mt_score * _WEIGHTS["min_temp"], mt_r),
        "freeze_thaw_cycles": (protocol.num_freeze_thaw_cycles, fc_score * _WEIGHTS["freeze_thaw_cycles"], fc_r),
        "ablation_fraction": (protocol.ablation_fraction, af_score * _WEIGHTS["ablation_fraction"], af_r),
        "checkpoint_combo": (protocol.combined_with_checkpoint_inhibitor, ck_score * _WEIGHTS["checkpoint_combo"], ck_r),
    }

    total = sum(v[1] for v in contributions.values())
    score_0_100 = round(total * 100, 1)

    caveats = [
        "This is a literature-parameterized heuristic score, not a model trained on "
        "outcome data -- no adequate public dataset of (freeze parameters -> measured "
        "immune outcome) currently exists to train one.",
        "Weights reflect the RELATIVE strength/consistency of evidence for each factor "
        "as assessed from the literature reviewed for this project, not a fitted "
        "statistical estimate.",
        "The ablation-fraction factor in particular reflects a genuinely contested "
        "question in the field and should not be read as a confident recommendation.",
        "Intended use: relative comparison between candidate protocols for the SAME "
        "tumor/patient context, not an absolute probability of clinical response.",
    ]

    return ImmuneEstimate(
        activation_score=score_0_100,
        contributing_factors=contributions,
        evidence_quality="literature-heuristic (not outcome-trained)",
        caveats=caveats,
    )
