"""Physics invariant integration tests.

Each test asserts a universally required property of any valid CARMA run —
no NaN/Inf, non-negative number density, supersaturated initial conditions
must produce particles, and growth/evaporation must not fire simultaneously
(they are opposing processes). Comparative or configuration-dependent
"physics" tests were removed because they don't cleanly isolate a single
process.

All tests require CARMAPY_RUN_INTEGRATION=1.
"""
import pytest
import numpy as np


def _build_tio2(path, example_levels, *, n_tstep=50, output_gap=5, dt=250):
    from carmapy.carmapy import Carma
    from carmapy.chemistry import populate_abundances_at_cloud_base
    P, T, kzz, mu = example_levels
    c = Carma(path)
    c.set_physical_params(surface_grav=31600, wt_mol=float(np.mean(mu)))
    c.set_atmospheric_parameters_from_defaults("Pure H2")
    c.set_stepping(dt=dt, output_gap=output_gap, n_tstep=n_tstep)
    c.add_gas("TiO2")
    c.add_hom_group("TiO2", 1e-8)
    c.add_P(P)
    c.add_T(T)
    c.add_kzz(kzz)
    c.calculate_z(mu)
    populate_abundances_at_cloud_base(c)
    return c


@pytest.mark.integration
def test_no_nans_or_infs(tmp_path, example_levels):
    """Numerical sanity: no NaN or Inf in any field."""
    from carmapy.results import Results
    c = _build_tio2(str(tmp_path / "nan_check"), example_levels)
    c.run(suppress_output=True)
    r = Results(c, read_diag=False)
    assert np.all(np.isfinite(r.numden)), "NaN or Inf in numden"
    assert np.all(np.isfinite(r.gas_abund)), "NaN or Inf in gas_abund"


@pytest.mark.integration
def test_numden_nonneg_all_timesteps(tmp_path, example_levels):
    """Physical invariant: number density cannot be negative at any timestep."""
    from carmapy.results import Results
    c = _build_tio2(str(tmp_path / "nonneg"), example_levels)
    c.run(suppress_output=True)
    r = Results(c, read_diag=False)
    assert np.all(r.numden >= 0), "Negative numden detected"


@pytest.mark.integration
def test_nucleation_creates_particles(tmp_path, example_levels):
    """Given a supersaturated initial condition (populate_abundances_at_cloud_base
    fills gas to the saturation limit), at least one timestep must show particles
    forming. If total numden is identical at all times, nucleation is dead."""
    from carmapy.results import Results
    c = _build_tio2(str(tmp_path / "nuc_parts"), example_levels)
    c.run(suppress_output=True)
    r = Results(c, read_diag=False)
    total_initial = np.sum(r.numden[:, :, :, 0])
    total_final = np.sum(r.numden[:, :, :, -1])
    assert total_final > total_initial, (
        "Total numden did not increase from initial state — nucleation is not firing"
    )


@pytest.mark.integration
def test_grow_and_evap_not_simultaneously_large(tmp_path, example_levels):
    """Growth and evaporation are opposing processes at any single point:
    a level/bin/timestep is either net-condensing or net-evaporating, not
    both at large magnitudes. Both rates being large simultaneously would
    indicate a logic bug in CARMA's grow/evap dispatcher."""
    from carmapy.results import Results
    c = _build_tio2(str(tmp_path / "grow_evap"), example_levels)
    c.run(suppress_output=True)
    r = Results(c, read_diag=True)
    grow = r.clouds["Pure TiO2"]["grow_gain_rate"]  # (NZ, NBIN, NT)
    evap = r.clouds["Pure TiO2"]["evap_gain_rate"]
    threshold = 0.1 * max(grow.max(), evap.max())
    both_large = np.sum((grow > threshold) & (evap > threshold))
    assert both_large == 0, (
        f"{both_large} level/bin/time points have both large growth and evaporation rates"
    )


def _build_silicate_closed_box(path, example_levels):
    """Mg2SiO4 and MgSiO3 on TiO2 and homogeneous SiO sharing the Mg and SiO
    reservoirs, with zero-flux boundaries for gas and particles."""
    from carmapy.carmapy import Carma
    from carmapy.chemistry import populate_abundances_at_cloud_base
    P, T, kzz, mu = example_levels
    c = Carma(path)
    c.set_physical_params(surface_grav=31600, wt_mol=float(np.mean(mu)))
    c.set_atmospheric_parameters_from_defaults("Pure H2")
    c.set_stepping(dt=250, output_gap=5, n_tstep=50)
    c.add_P(P)
    c.add_T(T)
    c.add_kzz(kzz)
    c.calculate_z(mu)
    c.add_hom_group("TiO2", 1e-8)
    r_het = 1e-8 * 2 ** (1 / 3)
    c.add_het_group("Mg2SiO4", "TiO2", r_het)
    c.add_het_group("MgSiO3", "TiO2", r_het)
    c.add_hom_group("SiO", 1e-8)
    c.set_gas_boundary_type("fixed_flux", "fixed_flux")
    c.set_cloud_boundary_type("fixed_flux", "fixed_flux")
    populate_abundances_at_cloud_base(c)
    return c


@pytest.mark.integration
def test_silicate_elemental_conservation(tmp_path, example_levels):
    """In the closed silicate box, column Mg and Si (gas plus
    stoichiometry-weighted condensate) stay constant."""
    from carmapy.carmapy import Carma
    from carmapy.chemistry import populate_abundances_at_cloud_base
    from carmapy.results import Results
    import scipy.constants as sc

    BK_CGS = 1.380649e-16  # Boltzmann constant [erg/K]
    AVOG = sc.N_A          # Avogadro's number [1/mol]

    c = _build_silicate_closed_box(str(tmp_path / "elemcons"), example_levels)
    c.run(suppress_output=True)

    r = Results(c, read_diag=True)

    n_air = r.P / (BK_CGS * r.T)                       # [#/cm3], (NZ,)
    dz = np.abs(np.gradient(np.asarray(r.Z)))          # [cm], (NZ,)

    igas = {name: i for i, name in enumerate(r.gas_names)}
    igrp = {name: i for i, name in enumerate(r.group_names)}

    # Moles of formula units per cm3 in each condensate: condensate mass /
    # condensate molar mass.  For a heterogeneous group the condensate is the
    # mantle, so the seed core's mass is excluded.
    groups = c.groups
    wtmol = {name: groups[name].wtmol for name in r.group_names}
    NT = r.numden.shape[-1]

    def condensate_moles(name, t):
        g = igrp[name]
        mass = r.numden[:, g, :, t] * r.rmass[:, g][None, :]
        if groups[name].mantle is not None:
            mass = mass * (1.0 - r.clouds[name]["coremass_frac"][:, :, t])
        return mass.sum(axis=1) / wtmol[name]

    totals = {"Mg": [], "Si": []}
    for t in range(NT):
        n_mg_gas = r.gas_abund[:, igas["Mg"], t] * 1e-6 * n_air
        n_sio_gas = r.gas_abund[:, igas["SiO"], t] * 1e-6 * n_air

        m_fo = condensate_moles("Mg2SiO4 on TiO2", t)
        m_en = condensate_moles("MgSiO3 on TiO2", t)
        m_sio = condensate_moles("Pure SiO", t)

        mg = (n_mg_gas / AVOG + 2 * m_fo + 1 * m_en)
        si = (n_sio_gas / AVOG + 1 * m_fo + 1 * m_en + 1 * m_sio)
        totals["Mg"].append(float((mg * dz).sum()))
        totals["Si"].append(float((si * dz).sum()))

    for elem in ("Mg", "Si"):
        series = np.array(totals[elem])
        assert series[0] > 0
        # The t=0 row precedes the engine's first-step gas adjustment, so
        # drift is measured from the first post-step output.  numden is
        # written to 3 significant figures, which limits resolution to ~1e-3.
        drift = np.abs(series[1:] - series[1]).max() / series[1]
        assert drift < 1e-3, (
            f"Column elemental {elem} drifts by {drift:.2e} "
            f"(> 0.1%) in a closed box after the first step: {series}"
        )
        total_drift = np.abs(series - series[0]).max() / series[0]
        assert total_drift < 1e-2, (
            f"Column elemental {elem} moves by {total_drift:.2e} "
            f"(> 1%) including the initial-output adjustment: {series}"
        )


# LLM-GENERATED — UNVERIFIED (tier: simulation)
# Remove this banner only after completing the checks in Assistant/POLICY.md §2.
def _ffloat(tok):
    """Parse a Fortran real, including the E-less form ('1.2-139') that an
    E-format edit descriptor writes for three-digit exponents."""
    import re
    return float(re.sub(r"(?<=\d)([+-]\d{3})$", r"E\1", tok))


def _reservoir_columns(c):
    """Column totals [mol/cm2 x const] of every non-water gas reservoir at each
    output, from the idiag = 1 full-precision dump: the gas, plus nu times
    the mantle of every condensate depleting it, plus seed mass for the seed's
    own reservoir.  The dump's air density carries a constant factor, which
    cancels in relative drift."""
    import os
    from carmapy.constants import gas_dict
    name = os.path.basename(c.name.rstrip("/"))
    path = os.path.join(c.name, f"fullprec_{name}.txt")
    NZ, NBIN = c.NZ, c.NBIN
    NELEM = len(c.elems)
    gnames = list(c.gases.keys())
    blocks, cur = [], None
    for line in open(path):
        if line.startswith("T"):
            if cur is not None:
                blocks.append(cur)
            cur = []
        else:
            cur.append([_ffloat(x) for x in line.split()])
    if cur is not None:
        blocks.append(cur)
    blocks = [b for b in blocks if len(b) == NZ + NZ * NBIN]
    out = {g: [] for g in gnames if g != "H2O"}
    for b in blocks:
        rho = np.array([row[0] for row in b[:NZ]])
        dz = np.array([row[1] for row in b[:NZ]])
        m = np.array(b[NZ:]).reshape(NBIN, NZ, NELEM + len(gnames))
        me, mg = m[..., :NELEM], m[0, :, NELEM:]
        for gas in out:
            tot = mg[:, gnames.index(gas)] * rho / gas_dict[gas]["wtmol_dif"]
            for gname, grp in c.groups.items():
                seed = gname.split(" on ")[-1].replace("Pure ", "")
                ce = grp.core.ielem - 1
                if grp.mantle is None:
                    if seed == gas:
                        tot = tot + me[:, :, ce].sum(axis=0) * rho / grp.wtmol
                    continue
                if seed == gas:
                    tot = tot + me[:, :, ce].sum(axis=0) * rho / c.groups["Pure " + seed].wtmol
                stofact = grp.stofact
                nu = (stofact.get(gas, 0) if isinstance(stofact, dict)
                      else stofact * (grp.gas.name == gas))
                if nu:
                    mantle = (me[:, :, grp.mantle.ielem - 1] - me[:, :, ce]).sum(axis=0)
                    tot = tot + nu * mantle * rho / grp.wtmol
            out[gas].append(float((tot * dz).sum()))
    return {g: np.array(v) for g, v in out.items()}


def _floor_clamps(c):
    """Number of gsolve floor clamps of gas that microphysics was consuming,
    from the idiag = 1 failure log.  Clamps of round-off negatives
    transported into empty layers are not counted."""
    import os
    name = os.path.basename(c.name.rstrip("/"))
    path = os.path.join(c.name, f"failures_{name}.txt")
    n = 0
    for line in open(path):
        p = line.split()
        if p and p[0] == "F" and p[8] == "1" and _ffloat(p[6]) < 0.0:
            n += 1
    return n


@pytest.mark.integration
def test_silicate_reservoirs_conserved_full_precision(tmp_path, example_levels):
    """With implicit_gas, every reservoir in the closed
    silicate box (gas + stoichiometric condensate + seed) holds to round-off
    at full precision, with no floor clamp of consumed gas."""
    c = _build_silicate_closed_box(str(tmp_path / "fpcons"), example_levels)
    c.implicit_gas = True
    c.idiag = 1
    c.run(suppress_output=True)

    totals = _reservoir_columns(c)
    for gas, series in totals.items():
        assert series[0] > 0
        drift = np.abs(series - series[0]).max() / series[0]
        assert drift < 1e-10, f"{gas} reservoir drifts by {drift:.2e} in a closed box: {series}"
    assert _floor_clamps(c) == 0, "consumed gas was floor-clamped (non-conservative fallback)"
