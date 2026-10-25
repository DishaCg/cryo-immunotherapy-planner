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


def test_gel_phantom_regression_pin():
    """Regression pin (NOT a paper number): single 17G probe, -140C, perfusion off,
    single 10 min freeze at the default 1.0 mm grid / 2 s step gives an ice-ball
    diameter ~30.6 mm and -40C diameter ~13.2 mm. Guards against silent solver
    changes. The paper's Sec. V-A uses the two-freeze benchmark at dr = 0.4 mm."""
    gel = TissueProperties(perfusion_rate=0.0)
    r = simulate([ProbeConfig(r_mm=0.0, z_mm=0.0)], tissue=gel, freeze_time_s=600.0,
                 dt_s=2.0, record_every_s=600.0)
    assert abs(2 * r.iceball_radius_mm[-1] - 30.6) < 0.3
    assert abs(2 * r.lethal_radius_mm[-1] - 13.2) < 0.3


def test_extra_isotherm_nesting():
    """-20C isotherm must lie between the 0C and -40C isotherms."""
    r = simulate([ProbeConfig(r_mm=0.0, z_mm=0.0)], freeze_time_s=300.0, dt_s=3.0,
                 record_every_s=300.0, domain_r_mm=30, domain_z_mm=50,
                 extra_isotherms=(-20.0,))
    m20 = r.extra_radii[-20.0][-1]
    assert r.lethal_radius_mm[-1] <= m20 <= r.iceball_radius_mm[-1]


# NOTE: the earlier "two probes produce a larger zone" test was removed. The solver
# is axisymmetric, so a probe at r > 0 is a RING source, not a second probe, and the
# old test only asserted ">=" against the single-probe case. Multi-probe behaviour
# needs a 3-D discretisation before it can be tested meaningfully.


if __name__ == "__main__":
    test_single_probe_monotonic_growth()
    print("PASS: test_single_probe_monotonic_growth")
    test_lethal_radius_smaller_than_iceball()
    print("PASS: test_lethal_radius_smaller_than_iceball")
    test_order_of_magnitude_vs_clinical_benchmark()
    print("PASS: test_order_of_magnitude_vs_clinical_benchmark")
    test_gel_phantom_regression_pin()
    print("PASS: test_gel_phantom_regression_pin")
    test_extra_isotherm_nesting()
    print("PASS: test_extra_isotherm_nesting")
    print("\nAll geometry model sanity tests passed.")
