"""
Sanity + regression tests for Module A (geometry_model).

These are NOT clinical validation (that requires gel-phantom / clinical
ground truth per the project proposal, Section 8). They check that the
model behaves physically sensibly and catch the kind of silent
numerical breakage (e.g. the earlier explicit-scheme stalling bug)
that isn't obvious just from the code compiling and running.
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from core.geometry_model import simulate, ProbeConfig, TissueProperties


def test_single_probe_monotonic_growth():
    """Ice-ball and lethal radius must grow monotonically during an
    active freeze -- if this fails, something is numerically wrong
    (e.g. the front stalling, as happened with the first explicit
    implementation)."""
    probes = [ProbeConfig(r_mm=0.0, z_mm=0.0)]
    result = simulate(probes, freeze_time_s=300.0, dt_s=3.0,
                       record_every_s=30.0, domain_r_mm=30, domain_z_mm=50)
    assert np.all(np.diff(result.iceball_radius_mm) >= -1e-6), \
        "ice-ball radius should never shrink during active freeze"
    assert np.all(np.diff(result.lethal_radius_mm) >= -1e-6), \
        "lethal radius should never shrink during active freeze"
    assert result.iceball_radius_mm[-1] > result.iceball_radius_mm[0] + 3, \
        "ice-ball should grow meaningfully over 5 minutes, not stall"


def test_lethal_radius_smaller_than_iceball():
    """The -40C isotherm must always be strictly inside the 0C isotherm
    -- if not, the isotherm-extraction logic is broken."""
    probes = [ProbeConfig(r_mm=0.0, z_mm=0.0)]
    result = simulate(probes, freeze_time_s=180.0, dt_s=3.0,
                       record_every_s=60.0, domain_r_mm=30, domain_z_mm=50)
    assert np.all(result.lethal_radius_mm <= result.iceball_radius_mm + 1e-9)


def test_order_of_magnitude_vs_clinical_benchmark():
    """Single-probe ice-ball diameter after ~8-10 min is clinically
    cited around 2.5-3.5 cm (radius ~12.5-17.5mm) for commercial
    cryoprobes. This model is not yet calibrated to match that exactly
    (see proposal Section 8), but a 10-minute single-probe freeze
    should land within the same order of magnitude, not off by 3x+
    as the earlier buggy explicit version was."""
    probes = [ProbeConfig(r_mm=0.0, z_mm=0.0)]
    result = simulate(probes, freeze_time_s=600.0, dt_s=3.0,
                       record_every_s=600.0, domain_r_mm=35, domain_z_mm=60)
    final_radius = result.iceball_radius_mm[-1]
    assert 6.0 <= final_radius <= 25.0, (
        f"single-probe 10-min ice-ball radius ({final_radius:.1f}mm) is "
        "outside a plausible clinical order-of-magnitude range (6-25mm); "
        "check for a regression before trusting downstream results"
    )


def test_two_probes_produce_larger_combined_zone():
    """Two probes placed close together should produce a larger fused
    ablation zone than either probe alone -- a basic multi-probe
    synergy sanity check."""
    single = simulate([ProbeConfig(r_mm=0.0, z_mm=0.0)],
                       freeze_time_s=300.0, dt_s=3.0, record_every_s=300.0,
                       domain_r_mm=30, domain_z_mm=50)

    two_probes = [
        ProbeConfig(r_mm=0.0, z_mm=0.0),
        ProbeConfig(r_mm=12.0, z_mm=0.0),
    ]
    result = simulate(two_probes, freeze_time_s=300.0, dt_s=3.0,
                       record_every_s=300.0, domain_r_mm=30, domain_z_mm=50)
    # radius measured from the FIRST probe's axis should be at least as
    # large with a second probe assisting nearby
    assert result.iceball_radius_mm[-1] >= single.iceball_radius_mm[-1] - 1e-6


if __name__ == "__main__":
    test_single_probe_monotonic_growth()
    print("PASS: test_single_probe_monotonic_growth")
    test_lethal_radius_smaller_than_iceball()
    print("PASS: test_lethal_radius_smaller_than_iceball")
    test_order_of_magnitude_vs_clinical_benchmark()
    print("PASS: test_order_of_magnitude_vs_clinical_benchmark")
    test_two_probes_produce_larger_combined_zone()
    print("PASS: test_two_probes_produce_larger_combined_zone")
    print("\nAll geometry model sanity tests passed.")
