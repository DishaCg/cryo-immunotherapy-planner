"""
Integration layer: joint protocol optimizer.

This module implements the core novel contribution of the associated
proposal. Existing planning tools optimize freeze protocols for a
single objective -- geometric tumor coverage (Module A's scope alone).
This layer searches candidate protocols and evaluates the trade-off
between:
  (a) ablation completeness/margin (from Module A, physics-grounded)
  (b) predicted immune activation potential (from Module B, a
      literature-derived heuristic; see its evidence-quality label)

Rather than combining these into a single composite score -- which
would obscure the difference in evidence basis between a
physics-grounded prediction and a literature heuristic behind one
number -- this module returns the Pareto-efficient set: protocols for
which no improvement in one objective is possible without sacrificing
the other. The trade-off is surfaced explicitly; the selection between
Pareto-optimal candidates is left to clinical or research judgment.
"""

from dataclasses import dataclass
from typing import List
import numpy as np

from .geometry_model import simulate, ProbeConfig, TissueProperties
from .immune_model import estimate as immune_estimate, FreezeProtocol


@dataclass
class CandidateProtocol:
    freeze_time_s: float
    num_freeze_thaw_cycles: int
    probe_temp_c: float
    combined_with_checkpoint_inhibitor: bool


@dataclass
class ProtocolEvaluation:
    protocol: CandidateProtocol
    lethal_radius_mm: float
    immune_score: float
    is_pareto_optimal: bool = False


def evaluate_protocol(
    protocol: CandidateProtocol,
    tumor_radius_mm: float,
    probe_r_mm: float = 0.0,
    probe_z_mm: float = 0.0,
    domain_r_mm: float = 35.0,
    domain_z_mm: float = 60.0,
) -> ProtocolEvaluation:
    """Run Module A for the geometric outcome and Module B for the
    immune-activation estimate on a single candidate protocol."""
    probe = ProbeConfig(r_mm=probe_r_mm, z_mm=probe_z_mm, temp_c=protocol.probe_temp_c)

    result = simulate(
        [probe],
        freeze_time_s=protocol.freeze_time_s,
        dt_s=3.0,
        record_every_s=protocol.freeze_time_s,
        domain_r_mm=domain_r_mm,
        domain_z_mm=domain_z_mm,
    )

    if protocol.num_freeze_thaw_cycles >= 2:
        thawed = simulate([probe], freeze_time_s=300.0, dt_s=3.0, record_every_s=300.0,
                           domain_r_mm=domain_r_mm, domain_z_mm=domain_z_mm,
                           initial_T=result.final_T, probe_active=False)
        result = simulate([probe], freeze_time_s=protocol.freeze_time_s, dt_s=3.0,
                           record_every_s=protocol.freeze_time_s,
                           domain_r_mm=domain_r_mm, domain_z_mm=domain_z_mm,
                           initial_T=thawed.final_T)

    lethal_radius = result.lethal_radius_mm[-1]
    ablation_fraction = min(lethal_radius / max(tumor_radius_mm, 1e-6), 1.5)  # allow >1 = full coverage + margin
    ablation_fraction_clamped = min(ablation_fraction, 1.0)

    # approximate cooling rate from probe temp / a nominal ramp time (~2 min to nadir, typical for argon systems)
    approx_cooling_rate = protocol.probe_temp_c / 2.0  # C per min, crude but documented approximation

    fp = FreezeProtocol(
        cooling_rate_c_per_min=approx_cooling_rate,
        min_temp_c=protocol.probe_temp_c,
        num_freeze_thaw_cycles=protocol.num_freeze_thaw_cycles,
        ablation_fraction=ablation_fraction_clamped,
        combined_with_checkpoint_inhibitor=protocol.combined_with_checkpoint_inhibitor,
    )
    immune = immune_estimate(fp)

    return ProtocolEvaluation(
        protocol=protocol,
        lethal_radius_mm=lethal_radius,
        immune_score=immune.activation_score,
    )


def _pareto_front(evals: List[ProtocolEvaluation]) -> List[ProtocolEvaluation]:
    """Mark protocols as Pareto-optimal: no other protocol beats them
    on BOTH objectives simultaneously (higher lethal_radius_mm AND
    higher immune_score)."""
    for e in evals:
        dominated = any(
            (o.lethal_radius_mm >= e.lethal_radius_mm and o.immune_score >= e.immune_score
             and (o.lethal_radius_mm > e.lethal_radius_mm or o.immune_score > e.immune_score))
            for o in evals if o is not e
        )
        e.is_pareto_optimal = not dominated
    return evals


def search_protocols(
    tumor_radius_mm: float,
    freeze_time_options_s: List[float] = (300.0, 600.0, 900.0),
    cycle_options: List[int] = (1, 2),
    checkpoint_options: List[bool] = (False, True),
    probe_temp_c: float = -140.0,
) -> List[ProtocolEvaluation]:
    """Grid search over a small, clinically meaningful protocol space
    and return every candidate annotated with whether it's on the
    Pareto-efficient frontier."""
    candidates = []
    for ft in freeze_time_options_s:
        for cycles in cycle_options:
            for combo in checkpoint_options:
                candidates.append(CandidateProtocol(
                    freeze_time_s=ft,
                    num_freeze_thaw_cycles=cycles,
                    probe_temp_c=probe_temp_c,
                    combined_with_checkpoint_inhibitor=combo,
                ))

    evals = [evaluate_protocol(c, tumor_radius_mm) for c in candidates]
    return _pareto_front(evals)
