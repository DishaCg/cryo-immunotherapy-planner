"""
Grid / time-step convergence and benchmark comparison for Module A (paper Sec. V-A).

    python -m analysis.convergence_study --quick     # ~1 min: dr = 0.8, 0.4 mm and the dt study
    python -m analysis.convergence_study             # adds dr = 0.2 mm and the benchmark protocol (~10+ min)

Notes
* dr = 0.8, 0.4, 0.2 mm all give the same effective probe radius (0.8 mm) so the
  sequence isolates discretisation error. (dr = 1.0 mm gives 1.0 mm; dr = 0.25 gives 0.75 mm.)
* Published values are the IceRod 1.5 single-needle side-view widths/lengths read from the
  isotherm table of Boston Scientific PI-814202-AB (jumbled text layer; verify against the figure).
"""
import sys
import numpy as np
from core.geometry_model import simulate, ProbeConfig, TissueProperties

GEL = TissueProperties(perfusion_rate=0.0)
PUBLISHED = {"width": {0: (34, 36), -20: (25, 26), -40: (16, 16)},
             "length": {0: (53, 55), -20: (45, 45), -40: (39, 40)}}

def axial_len(res, target):
    z = res.z_centers_mm; T = res.final_T[0, :]
    idx = np.where(T <= target)[0]
    if len(idx) == 0: return 0.0
    lo, hi = idx[0], idx[-1]
    def edge(a, b):
        return z[a] if T[a] == T[b] else z[a] + (target - T[a]) / (T[b] - T[a]) * (z[b] - z[a])
    return (edge(hi, hi + 1) if hi < len(z) - 1 else z[hi]) - (edge(lo, lo - 1) if lo > 0 else z[lo])

def single(dr, dt, L=17.0, dom=(40, 70)):
    r = simulate([ProbeConfig(0, 0, active_length_mm=L)], tissue=GEL, domain_r_mm=dom[0], domain_z_mm=dom[1],
                 dr_mm=dr, dt_s=dt, freeze_time_s=600, record_every_s=600, extra_isotherms=(-20.0,))
    return 2 * r.iceball_radius_mm[-1], 2 * r.extra_radii[-20.0][-1], 2 * r.lethal_radius_mm[-1]

def benchmark(dr, L=17.0, dom=(40, 70), dt=2.0):
    kw = dict(tissue=GEL, domain_r_mm=dom[0], domain_z_mm=dom[1], dr_mm=dr, dt_s=dt, extra_isotherms=(-20.0,))
    P = [ProbeConfig(0, 0, active_length_mm=L)]
    a = simulate(P, freeze_time_s=600, record_every_s=600, **kw)
    b = simulate(P, freeze_time_s=300, record_every_s=300, initial_T=a.final_T, probe_active=False, **kw)
    c = simulate(P, freeze_time_s=600, record_every_s=600, initial_T=b.final_T, **kw)
    return dict(w0=2 * c.iceball_radius_mm[-1], w20=2 * c.extra_radii[-20.0][-1], w40=2 * c.lethal_radius_mm[-1],
                l0=axial_len(c, 0.0), l20=axial_len(c, -20.0), l40=axial_len(c, -40.0))

if __name__ == "__main__":
    quick = "--quick" in sys.argv
    print("Single 600 s freeze, gel phantom, 17G probe: widths (0C, -20C, -40C) in mm")
    for dr in ((0.8, 0.4) if quick else (0.8, 0.4, 0.2)):
        print("  dr=%.1f mm  dt=2 s :" % dr, tuple(round(float(x), 2) for x in single(dr, 2.0)))
    for dt in (4.0, 2.0, 1.0, 0.5):
        print("  dr=1.0 mm  dt=%.1f s:" % dt, tuple(round(float(x), 2) for x in single(1.0, dt)))
    if not quick:
        print("\\nBenchmark protocol (2 x 600 s freeze, 300 s passive thaw), dr = 0.4 mm:")
        r = benchmark(0.4); print({k: round(float(v), 1) for k, v in r.items()})
        print("published widths", PUBLISHED["width"], "lengths", PUBLISHED["length"])
