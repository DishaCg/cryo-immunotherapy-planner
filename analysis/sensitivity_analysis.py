"""
Joint-optimisation example (Table 4) and sensitivity analysis (Section V.D).
Run from the repo root:
    python -m analysis.sensitivity_analysis
"""
import functools, json
import numpy as np
from core.geometry_model import simulate, ProbeConfig
from core.immune_model import estimate, FreezeProtocol

TUMOR_R_MM = 6.8   # chosen so candidate coverages straddle the 95% threshold; see paper Sec. V-C
DR_MM = 0.4        # grid spacing; see convergence study (paper Sec. V-A)
FREEZE_S = (300, 450, 600, 900); CYCLES = (1, 2); PROBE_T = -140.0
W0 = dict(cooling=0.20, nadir=0.15, cycles=0.20, ablation=0.15, checkpoint=0.30)

@functools.lru_cache(maxsize=None)
def lethal_radius(freeze_s, cycles, temp_c, domain=(35.0, 60.0), dr=DR_MM):
    p = ProbeConfig(0, 0, temp_c=temp_c); kw = dict(domain_r_mm=domain[0], domain_z_mm=domain[1], dt_s=3.0, dr_mm=dr)
    res = simulate([p], freeze_time_s=freeze_s, record_every_s=freeze_s, **kw)
    if cycles >= 2:
        th = simulate([p], freeze_time_s=300.0, record_every_s=300.0, initial_T=res.final_T, probe_active=False, **kw)
        res = simulate([p], freeze_time_s=freeze_s, record_every_s=freeze_s, initial_T=th.final_T, **kw)
    return float(res.lethal_radius_mm[-1])

def score(rate, nadir, cycles, frac, ck, w=W0, depot=0.75, mid=0.85, hi=0.95):
    """Same factor functions as core/immune_model.py, with the contested
    ablation-fraction parameters (depot, mid, hi) exposed for sensitivity runs."""
    f = dict(cooling=min(abs(rate) / 40.0, 1.0), nadir=min(max(0.0, -nadir) / 100.0, 1.0),
             cycles={0: 0.0, 1: 0.55, 2: 1.0}.get(cycles, 1.0 if cycles > 2 else 0.0),
             ablation=(frac / 0.3 * 0.4 if frac < 0.3 else depot if frac > hi else mid),
             checkpoint=1.0 if ck else 0.15)
    return round(sum(f[k] * w[k] for k in w) * 100, 1)

def _selfcheck():
    """score() above must agree with core.immune_model.estimate() at default parameters."""
    for cyc in (1, 2):
        for frac in (0.2, 0.6, 0.97):
            for ck in (False, True):
                ref = estimate(FreezeProtocol(-70.0, -140.0, cyc, frac, ck)).activation_score
                assert abs(ref - score(-70.0, -140.0, cyc, frac, ck)) < 1e-9, (cyc, frac, ck)

def pareto(points):
    return [not any(l2 >= l and s2 >= s and (l2 > l or s2 > s) for j, (l2, s2) in enumerate(points) if j != i)
            for i, (l, s) in enumerate(points)]

if __name__ == "__main__":
    _selfcheck()
    rng = np.random.default_rng(42)
    cands = [(f, c) for f in FREEZE_S for c in CYCLES]
    L = {k: lethal_radius(k[0], k[1], PROBE_T) for k in cands}   # checkpoint inhibitor held fixed (= True)
    def front(w=W0, **kw):
        pts = [(L[k], score(PROBE_T / 2, PROBE_T, k[1], min(L[k] / TUMOR_R_MM, 1.0), True, w=w, **kw)) for k in cands]
        return pts, pareto(pts)
    pts, pf = front(); base = frozenset(k for k, p in zip(cands, pf) if p)
    print("Table 4 (tumor R = %.1f mm, checkpoint inhibitor fixed):" % TUMOR_R_MM)
    for k, p, o in zip(cands, pts, pf): print(k, round(p[0], 2), p[1], "PARETO" if o else "")
    keys = list(W0); w0 = np.array([W0[k] for k in keys])
    def run(n, conc, depot=None, hi=None):
        same = single = 0
        for _ in range(n):
            w = dict(zip(keys, rng.dirichlet(conc * w0))); kw = {}
            if depot: kw["depot"] = rng.uniform(*depot)
            if hi: kw["hi"] = rng.uniform(*hi)
            _, pf_ = front(w, **kw); s = frozenset(k for k, p in zip(cands, pf_) if p)
            same += s == base; single += len(s) == 1
        return dict(P_front_unchanged=same / n, P_single_point=single / n)
    print("weights only, conc=50 :", run(5000, 50))
    print("weights only, conc=10 :", run(5000, 10))
    print("+ depot in U(.55,.84) :", run(5000, 10, depot=(0.55, 0.84)))
    print("+ depot U(.55,1), hi U(.90,.99):", run(5000, 10, depot=(0.55, 1.0), hi=(0.90, 0.99)))
    print("depot sweep (value, #Pareto):", [(round(float(d), 2), sum(front(depot=float(d))[1])) for d in np.arange(0.5, 1.001, 0.05)])
