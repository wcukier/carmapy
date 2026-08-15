"""The Python side of the radiative coupling, end to end.

``test_rce.py`` pins the physics of ``carma_rce`` against numpy references;
these tests pin the plumbing around it -- that ``set_radiation`` reaches the
Fortran, that the temperature output parses back into the shapes it claims, and
that a configured run survives ``load_carma``.

They run the real binary, so they need ``CARMAPY_RUN_INTEGRATION=1``. The ck
table is the grey one from ``test_rce``: it needs no Sonora download, and with
a wavelength-independent opacity the emergent flux can be compared against
``sigma*Teff^4`` directly.
"""
import os

import numpy as np
import pytest

from .test_rce import write_grey_ck

pytestmark = pytest.mark.integration

T_INT = 1000.0


@pytest.fixture(scope="module")
def ck_table(tmp_path_factory):
    path = tmp_path_factory.mktemp("ck") / "grey.ck"
    write_grey_ck(str(path), kappa_ln=-57.0)
    return str(path)


@pytest.fixture
def coupled_carma(tmp_path, example_levels, ck_table):
    """A minimal TiO2 run with the temperature profile radiatively coupled."""
    from carmapy.carmapy import Carma
    from carmapy.chemistry import populate_abundances_at_cloud_base

    P, T, kzz, mu = example_levels
    c = Carma(str(tmp_path / "coupled"))
    c.set_physical_params(surface_grav=31600, wt_mol=float(np.mean(mu)))
    c.set_atmospheric_parameters_from_defaults("Pure H2")
    c.set_stepping(dt=250, output_gap=5, n_tstep=20)
    c.add_gas("TiO2")
    c.add_hom_group("TiO2", 1e-8)
    c.add_P(P)
    c.add_T(T)
    c.add_kzz(kzz)
    c.calculate_z(mu)
    populate_abundances_at_cloud_base(c)
    c.set_radiation(t_int=T_INT, ck_table_path=ck_table,
                    accel=20.0, dT_max=5.0)
    return c


def test_set_radiation_forces_t_evolves(coupled_carma):
    """Without t_evolves the growth kernels are cached from step 1 and the
    coupling is silently inert, so enabling radiation must set it."""
    assert coupled_carma.t_evolves is True


def test_set_radiation_rejects_accel_in_physical_mode(coupled_carma, ck_table):
    with pytest.raises(ValueError, match="only meaningful in equilibrium"):
        coupled_carma.set_radiation(t_int=T_INT, ck_table_path=ck_table,
                                    mode="physical", accel=20.0)


def test_namelist_has_no_radiative_convective_boundary(coupled_carma):
    """The boundary is diagnosed, not configured, so nothing may carry it in.

    A stale ``p_rcb`` in the group would make the Fortran namelist read fail,
    which falls back to ``do_radiation = 0`` -- radiation silently off rather
    than an error.
    """
    import f90nml

    coupled_carma.run(suppress_output=True)
    nml = f90nml.read(os.path.join(coupled_carma.name, "inputs", "input.nml"))

    assert nml["radiation"]["do_radiation"] == 1
    assert nml["radiation"]["t_int"] == pytest.approx(T_INT)
    assert "p_rcb" not in nml["radiation"]


def test_radiation_survives_load_carma(coupled_carma):
    """A restarted run must reproduce the same coupling, not a default one."""
    from carmapy.base import load_carma

    coupled_carma.run(suppress_output=True)
    reloaded = load_carma(coupled_carma.name, restart=1)

    assert reloaded.t_int == pytest.approx(T_INT)
    assert not hasattr(reloaded, "p_rcb")
    assert reloaded.rad_mode == "equilibrium"
    assert reloaded.rad_accel == pytest.approx(20.0)
    assert reloaded.cloud_rad is True
    assert reloaded.t_evolves is True
    assert os.path.exists(reloaded.ck_table_path)


def test_temperature_output_parses(coupled_carma):
    from carmapy.results import Results

    coupled_carma.run(suppress_output=True)
    r = Results(coupled_carma)

    nt = r.T_history.shape[1]
    assert nt == r.ts.size
    assert r.T_history.shape == (coupled_carma.NZ, nt)
    assert r.dTdt.shape == (coupled_carma.NZ, nt)
    assert r.flux_net.shape == (coupled_carma.NZ + 1, nt)
    assert r.F_TOA.shape == (nt,)
    assert r.t_int == pytest.approx(T_INT)
    assert r.convective.shape == (coupled_carma.NZ, nt)
    assert r.convective.dtype == bool
    assert r.nzone.shape == (nt,)
    assert r.conv_resid.shape == (nt,)
    assert r.P_levels_rad.shape == (coupled_carma.NZ + 1,)

    assert np.all(r.T_history > 0)
    assert np.all(np.isfinite(r.flux_net))
    # F_TOA is the top level of flux_net, so the two must agree exactly.
    assert np.allclose(r.F_TOA, r.flux_net[-1, :], rtol=0, atol=0)


def test_temperature_actually_evolves(coupled_carma):
    """The whole point of the coupling: t(:) must not be the input profile."""
    from carmapy.results import Results

    coupled_carma.run(suppress_output=True)
    r = Results(coupled_carma)

    assert not np.allclose(r.T_history[:, 0], r.T_history[:, -1])
    # nzone counts the mask's contiguous runs, so it has to agree with the mask
    # it is derived from at every output.
    for it in range(r.nzone.size):
        assert len(r.convective_zones(it)) == r.nzone[it]


def test_tabulated_adiabat_end_to_end(coupled_carma):
    """A coupled run with ``adiabat="table"`` must put its convective interior
    on the tabulated adiabat, not the analytic one.

    This is the whole Fortran path in one assertion: the table is found by the
    path in the namelist, read and interpolated, and the temperature ratio
    across each convective interface is the one the *tabulated* gradient
    demands rather than the analytic fit's.
    """
    from carmapy.adiabat import grad_ad
    from carmapy.constants import (PARMENTIER_A_COEFF as PARM_A,
                                   PARMENTIER_B_COEFF as PARM_B)
    from carmapy.results import Results

    coupled_carma.adiabat = "table"
    coupled_carma.run(suppress_output=True)
    r = Results(coupled_carma)

    mask = r.convective[:, -1]
    inside = mask[:-1] & mask[1:]
    assert inside.sum() > 2, "no convective interfaces to check"

    t = r.T_history[:, -1]
    p, pl = r.P, r.P_levels_rad
    dp = pl[:-1] - pl[1:]
    t_bar = (dp[:-1] * t[:-1] + dp[1:] * t[1:]) / (dp[:-1] + dp[1:])

    # The midpoint rule the Fortran's rce_pfact uses, both ways.
    ratio = t[:-1] / t[1:]
    tab = (p[:-1] / p[1:]) ** grad_ad(t_bar, np.sqrt(p[:-1] * p[1:]))
    fit = (p[:-1] / p[1:]) ** (PARM_A - PARM_B * t_bar)

    # The adjustment stops after a fixed number of sweeps rather than at an
    # exactly neutral column, so the floor here is its own residual -- a few
    # parts in 1e6 -- not machine precision.
    np.testing.assert_allclose(ratio[inside], tab[inside], rtol=1e-4)

    # And the two gradients genuinely differ here, by a good margin over that
    # tolerance, or the test would pass just as well with the flag ignored.
    # The margin is only ~45x because this run's convective layers sit near the
    # top of the column, where the table is clamped at its 0.01 bar floor and
    # the dissociation flattening that separates the two has not switched on.
    separation = np.abs(tab[inside] / fit[inside] - 1).max()
    assert separation > 2e-3, f"table and fit differ by only {separation:.2e}"


def test_cadence_is_reported(coupled_carma):
    from carmapy.results import Results

    coupled_carma.run(suppress_output=True)
    r = Results(coupled_carma)

    assert r.rad_interval.shape == r.F_TOA.shape
    assert np.all(r.rad_interval >= 1)
    assert np.all(r.rad_interval <= coupled_carma.rad_gap_max)


def test_final_profiles_track_the_coupled_run(coupled_carma):
    """The PICASO writers must not be handed the input profile and the grid
    calculate_z built from it -- the Fortran moved both."""
    from carmapy.results import Results

    coupled_carma.run(suppress_output=True)
    r = Results(coupled_carma)

    tl = r.final_T_levels()
    assert tl.shape == (coupled_carma.NZ + 1,)
    assert np.all(tl > 0)
    assert not np.allclose(tl, coupled_carma.T_levels)

    zl = r.final_z_levels()
    assert zl.shape == (coupled_carma.NZ + 1,)
    assert zl[0] == 0
    assert np.all(np.diff(zl) > 0)
    assert not np.allclose(zl, coupled_carma.z_levels)


def test_final_profiles_are_the_inputs_when_uncoupled(tmp_path, example_levels):
    from carmapy.carmapy import Carma
    from carmapy.chemistry import populate_abundances_at_cloud_base
    from carmapy.results import Results

    P, T, kzz, mu = example_levels
    c = Carma(str(tmp_path / "plain"))
    c.set_physical_params(surface_grav=31600, wt_mol=float(np.mean(mu)))
    c.set_atmospheric_parameters_from_defaults("Pure H2")
    c.set_stepping(dt=250, output_gap=5, n_tstep=10)
    c.add_gas("TiO2")
    c.add_hom_group("TiO2", 1e-8)
    c.add_P(P)
    c.add_T(T)
    c.add_kzz(kzz)
    c.calculate_z(mu)
    populate_abundances_at_cloud_base(c)
    c.run(suppress_output=True)

    r = Results(c)
    assert r.final_T_levels() is c.T_levels
    assert r.final_z_levels() is c.z_levels


def test_uncoupled_run_has_no_temperature_output(tmp_path, example_levels):
    """Radiation off must leave the outputs exactly as they were."""
    from carmapy.carmapy import Carma
    from carmapy.chemistry import populate_abundances_at_cloud_base
    from carmapy.results import Results

    P, T, kzz, mu = example_levels
    c = Carma(str(tmp_path / "uncoupled"))
    c.set_physical_params(surface_grav=31600, wt_mol=float(np.mean(mu)))
    c.set_atmospheric_parameters_from_defaults("Pure H2")
    c.set_stepping(dt=250, output_gap=5, n_tstep=10)
    c.add_gas("TiO2")
    c.add_hom_group("TiO2", 1e-8)
    c.add_P(P)
    c.add_T(T)
    c.add_kzz(kzz)
    c.calculate_z(mu)
    populate_abundances_at_cloud_base(c)
    c.run(suppress_output=True)

    name = os.path.basename(c.name)
    assert not os.path.exists(
        os.path.join(c.name, f"temperature_{name}.txt"))

    r = Results(c)
    assert r.T_history is None
    assert r.F_TOA is None
