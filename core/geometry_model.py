"""
Module A: ice-ball / lethal-isotherm geometry prediction.

Solves the axisymmetric (r, z) Pennes bioheat equation, using an
enthalpy (apparent heat capacity) treatment of latent heat, with an
implicit finite-volume scheme:

    rho * c_eff(T) * dT/dt = div(k(T) grad T)
                              - w_b(T) * rho_b * c_b * (T - T_arterial)
                              + Q_metabolic

Why implicit finite-volume: latent heat is concentrated in a narrow
(~7C) mushy zone, which makes the apparent heat capacity spike sharply
right at the freezing front. That makes the PDE stiff there, and an
explicit scheme needs an impractically small time step to stay stable
-- on a normal grid this basically shows up as the ice front stalling
instead of growing. Implicit time-stepping gets rid of that
restriction. Cell-centered finite volumes were the natural choice on
top of that since they handle r=0 without any special-casing and
conserve energy exactly.

Two isotherms get tracked at the probe's mid-plane:
  - 0C   : the ice-ball boundary you'd actually see on imaging
  - -40C : the "lethal" isotherm, i.e. the temperature associated with
           complete tumor cell death, which is what actually matters
           for margin planning

Passes the regression tests in tests/test_geometry_model.py and is in
the right order of magnitude against published single-probe benchmark
data -- see CALIBRATION_NOTES.md for the full comparison. Not yet
fitted to gel-phantom or clinical measurements.
"""

from dataclasses import dataclass
import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import spsolve


# ----------------------------------------------------------------------
# Tissue and probe parameter models
# ----------------------------------------------------------------------

@dataclass
class TissueProperties:
    """Thermal properties. Defaults are representative soft-tissue /
    liver-like values assembled from standard cryoablation bioheat
    literature (Rossi & Rabin 2007; Baust et al. 2014 mechanism
    reviews) -- NOT patient-specific and should be swapped for
    tissue-specific values (renal, breast, etc.) before real use.
    """
    rho: float = 1060.0            # tissue density, kg/m^3
    c_unfrozen: float = 3600.0     # specific heat, unfrozen, J/kg/K
    c_frozen: float = 1800.0       # specific heat, frozen, J/kg/K
    k_unfrozen: float = 0.512      # thermal conductivity, unfrozen, W/m/K
    k_frozen: float = 2.0          # thermal conductivity, frozen, W/m/K
    latent_heat: float = 250000.0  # latent heat of fusion, J/kg
    T_liquidus: float = -1.0       # freezing onset, deg C
    T_solidus: float = -8.0        # fully frozen, deg C (mushy zone width)
    T_body: float = 37.0           # baseline tissue / arterial temp, deg C
    perfusion_rate: float = 0.008  # blood perfusion, 1/s (unfrozen only)
    rho_blood: float = 1060.0
    c_blood: float = 3600.0
    Q_metabolic: float = 400.0     # W/m^3, small vs. freezing terms


@dataclass
class ProbeConfig:
    r_mm: float                     # radial position of probe shaft, mm
    z_mm: float                     # axial center of the active freeze zone, mm
    diameter_mm: float = 1.7        # cryoprobe shaft diameter
    active_length_mm: float = 17.0  # active freeze-zone length along shaft
    # Commercial cryoprobes have an active freezing section on the order
    # of 1.5-3 cm, not a point tip -- a point-sink approximation
    # drastically under-predicts total heat extraction capacity.
    temp_c: float = -140.0          # tip temperature during active freeze


@dataclass
class SimulationResult:
    times_s: np.ndarray
    iceball_radius_mm: np.ndarray      # 0C isotherm at active-zone mid-plane
    lethal_radius_mm: np.ndarray       # -40C isotherm, same cross-section
    final_T: np.ndarray                # (nr, nz) grid, deg C
    r_centers_mm: np.ndarray
    z_centers_mm: np.ndarray
    extra_radii: dict = None           # optional {temp_C: radius_mm array} for extra_isotherms


# ----------------------------------------------------------------------
# Enthalpy-method effective properties (evaluated at previous time step
# -- this is what makes the scheme "semi-implicit": nonlinear material
# properties are linearized/lagged, but diffusion itself is implicit)
# ----------------------------------------------------------------------

def _effective_c(T, tp: TissueProperties):
    c = np.where(T >= tp.T_liquidus, tp.c_unfrozen, tp.c_frozen)
    mushy = (T < tp.T_liquidus) & (T > tp.T_solidus)
    width = tp.T_liquidus - tp.T_solidus
    latent_contribution = tp.latent_heat / width
    c = np.where(mushy, c + latent_contribution, c)
    return c


def _effective_k(T, tp: TissueProperties):
    frozen_frac = np.clip((tp.T_liquidus - T) / (tp.T_liquidus - tp.T_solidus), 0, 1)
    # Ice conductivity is not constant with temperature -- it increases
    # substantially as temperature drops further below 0C (roughly
    # k_ice(T) ~ 2.2 W/m/K at 0C rising toward ~3.5-4.0 W/m/K near -100C;
    # e.g. Slack 1980 ice thermal conductivity data). tp.k_frozen is
    # treated as the value at T_solidus; conductivity increases further
    # for fully-frozen tissue colder than that, up to a capped ceiling.
    k_frozen_local = tp.k_frozen + np.clip((tp.T_solidus - T) / 100.0, 0, 1) * (3.8 - tp.k_frozen)
    return tp.k_unfrozen + frozen_frac * (k_frozen_local - tp.k_unfrozen)


def _perfusion_active(T, tp: TissueProperties):
    return (T >= tp.T_liquidus).astype(float) * tp.perfusion_rate


# ----------------------------------------------------------------------
# Solver (cell-centered axisymmetric finite volume, implicit in time)
# ----------------------------------------------------------------------

def simulate(
    probes: list[ProbeConfig],
    tissue: TissueProperties = TissueProperties(),
    domain_r_mm: float = 40.0,
    domain_z_mm: float = 70.0,
    dr_mm: float = 1.0,
    freeze_time_s: float = 600.0,
    dt_s: float = 2.0,
    record_every_s: float = 30.0,
    initial_T: np.ndarray = None,
    probe_active: bool = True,
    extra_isotherms: tuple = (),
) -> SimulationResult:
    """
    initial_T: optional (nr, nz) starting temperature field, for chaining
        multi-phase protocols (e.g. freeze -> thaw -> freeze) by passing
        the previous phase's `result.final_T` in as the next phase's start.
    probe_active: if False, probes are NOT clamped to temp_c and instead
        left to evolve passively (used to simulate a thaw phase without
        active warming -- tissue re-warms via perfusion/conduction from
        surrounding body-temperature tissue, which is how passive thaw
        actually works clinically).
    extra_isotherms: optional extra isotherm temperatures (deg C), e.g.
        (-20.0,); their mid-plane radii come back in result.extra_radii.

    Matrix assembly is vectorised (no Python double loop). It is numerically
    identical to the earlier loop version (max |dT| ~ 1e-13 on every path
    tested) but runs in seconds instead of 10-15 s per protocol.
    """
    dr = dr_mm / 1000.0
    nr = int(round(domain_r_mm / dr_mm))
    nz = int(round(domain_z_mm / dr_mm))

    r = (np.arange(nr) + 0.5) * dr
    z = (np.arange(nz) + 0.5) * dr - domain_z_mm / 2000.0
    R, Z = np.meshgrid(r, z, indexing="ij")

    r_face = np.arange(nr + 1) * dr
    cell_vol = np.pi * (r_face[1:] ** 2 - r_face[:-1] ** 2) * dr
    cell_vol_grid = np.repeat(cell_vol[:, None], nz, axis=1)

    face_area_r = 2 * np.pi * r_face * dr
    face_area_z = 2 * np.pi * r * dr

    T = np.full((nr, nz), tissue.T_body, dtype=float) if initial_T is None else initial_T.copy()

    probe_mask = np.zeros((nr, nz), dtype=bool)
    probe_temp_grid = np.zeros((nr, nz), dtype=float)
    for p in probes:
        pr = p.r_mm / 1000.0
        pz = p.z_mm / 1000.0
        rad = max(p.diameter_mm / 2000.0, dr)
        half_len = max(p.active_length_mm / 2000.0, dr / 2)
        m = (np.abs(R - pr) <= rad) & (np.abs(Z - pz) <= half_len)
        if probe_active:
            probe_mask |= m
            probe_temp_grid[m] = p.temp_c
    if probe_active:
        T[probe_mask] = probe_temp_grid[probe_mask]

    n_steps = int(round(freeze_time_s / dt_s))
    record_every = max(1, int(round(record_every_s / dt_s)))

    times_record, iceball_record, lethal_record = [], [], []
    extra_record = {x: [] for x in extra_isotherms}
    z0_idx = np.argmin(np.abs(z))

    N = nr * nz

    def idx(i, j):
        return i * nz + j

    for step in range(n_steps + 1):
        t = step * dt_s

        if step % record_every == 0 or step == n_steps:
            profile = T[:, z0_idx]
            ice_r = _isotherm_radius(r, profile, 0.0)
            lethal_r = _isotherm_radius(r, profile, -40.0)
            times_record.append(t)
            iceball_record.append(ice_r * 1000.0)
            lethal_record.append(lethal_r * 1000.0)
            for x in extra_isotherms:
                extra_record[x].append(_isotherm_radius(r, profile, x) * 1000.0)

        if step == n_steps:
            break

        c_eff = _effective_c(T, tissue)
        k_eff = _effective_k(T, tissue)
        w_b = _perfusion_active(T, tissue)

        # --- vectorised assembly of the implicit system (same stencil as before) ---
        vol = cell_vol_grid
        accum = tissue.rho * c_eff * vol / dt_s
        diag = accum.copy()
        rhs = accum * T
        perf_coeff = w_b * tissue.rho_blood * tissue.c_blood * vol
        diag += perf_coeff
        rhs += perf_coeff * tissue.T_body + tissue.Q_metabolic * vol

        ID = np.arange(N).reshape(nr, nz)
        rows, cols, vals = [], [], []

        # radial face between i-1 and i (i > 0)
        trans = 0.5 * (k_eff[1:, :] + k_eff[:-1, :]) * face_area_r[1:nr, None] / dr
        diag[1:, :] += trans
        rows.append(ID[1:, :].ravel()); cols.append(ID[:-1, :].ravel()); vals.append(-trans.ravel())

        # radial face between i and i+1 (outer boundary held at body temperature)
        k_next = np.empty_like(k_eff)
        k_next[:-1] = k_eff[1:]
        k_next[-1] = tissue.k_unfrozen
        trans = 0.5 * (k_eff + k_next) * face_area_r[1:, None] / dr
        diag += trans
        rows.append(ID[:-1, :].ravel()); cols.append(ID[1:, :].ravel()); vals.append(-trans[:-1].ravel())
        rhs[-1, :] += trans[-1] * tissue.T_body

        # axial faces (bottom/top boundaries held at body temperature, half-cell distance)
        trans = 0.5 * (k_eff[:, 1:] + k_eff[:, :-1]) * face_area_z[:, None] / dr
        diag[:, 1:] += trans
        rows.append(ID[:, 1:].ravel()); cols.append(ID[:, :-1].ravel()); vals.append(-trans.ravel())
        diag[:, :-1] += trans
        rows.append(ID[:, :-1].ravel()); cols.append(ID[:, 1:].ravel()); vals.append(-trans.ravel())
        for edge in (0, -1):
            tb = k_eff[:, edge] * face_area_z / (dr / 2)
            diag[:, edge] += tb
            rhs[:, edge] += tb * tissue.T_body

        rows = np.concatenate(rows); cols = np.concatenate(cols); vals = np.concatenate(vals)
        keep = ~probe_mask.ravel()[rows]                 # probe cells: identity rows (Dirichlet)
        rows, cols, vals = rows[keep], cols[keep], vals[keep]
        d = np.where(probe_mask, 1.0, diag).ravel()
        bvec = np.where(probe_mask, probe_temp_grid, rhs).ravel()
        rows = np.concatenate([rows, np.arange(N)]); cols = np.concatenate([cols, np.arange(N)])
        vals = np.concatenate([vals, d])
        A = csr_matrix(coo_matrix((vals, (rows, cols)), shape=(N, N)))

        T = spsolve(A, bvec).reshape(nr, nz)

    return SimulationResult(
        times_s=np.array(times_record),
        iceball_radius_mm=np.array(iceball_record),
        lethal_radius_mm=np.array(lethal_record),
        final_T=T,
        r_centers_mm=r * 1000.0,
        z_centers_mm=z * 1000.0,
        extra_radii={x: np.array(v) for x, v in extra_record.items()},
    )


def _isotherm_radius(r, profile_1d, target_temp):
    below = profile_1d <= target_temp
    if not below.any():
        return 0.0
    last_idx = np.where(below)[0][-1]
    if last_idx == len(r) - 1:
        return r[last_idx]
    r0, r1 = r[last_idx], r[last_idx + 1]
    T0, T1 = profile_1d[last_idx], profile_1d[last_idx + 1]
    if T1 == T0:
        return r0
    frac = (target_temp - T0) / (T1 - T0)
    return r0 + frac * (r1 - r0)
