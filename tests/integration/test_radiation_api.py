"""The Python side of the radiative coupling, end to end: ``set_radiation``
reaches the Fortran, the temperature output parses, and a configured run
survives ``load_carma``.  Uses the grey ck table from ``test_rce``.
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
    c.set_radiation(t_int=T_INT, ck_table_path=ck_table, accel=20.0)
    return c


def test_set_radiation_forces_t_evolves(coupled_carma):
    """Enabling radiation sets t_evolves, without which growth kernels are
    cached from step 1."""
    assert coupled_carma.t_evolves is True


def test_namelist_has_no_radiative_convective_boundary(coupled_carma):
    """The radiation namelist has radiation on and carries no ``p_rcb``; an
    unknown key would fail the Fortran namelist read and turn radiation off."""
    import f90nml

    coupled_carma.run(suppress_output=True)
    nml = f90nml.read(os.path.join(coupled_carma.name, "inputs", "input.nml"))

    assert nml["radiation"]["do_radiation"] == 1
    assert nml["radiation"]["t_int"] == pytest.approx(T_INT)
    assert "p_rcb" not in nml["radiation"]


def test_radiation_survives_load_carma(coupled_carma):
    """A restarted run reproduces the configured coupling."""
    from carmapy.base import load_carma

    coupled_carma.run(suppress_output=True)
    reloaded = load_carma(coupled_carma.name, restart=1)

    rad, configured = reloaded.radiation, coupled_carma.radiation
    assert rad.t_int == pytest.approx(T_INT)
    assert not hasattr(rad, "p_rcb")
    assert rad.accel == pytest.approx(20.0)
    assert rad.cloud_rad is True
    assert reloaded.t_evolves is True
    assert os.path.exists(rad.ck_table_path)

    assert rad.self_consistent_kzz == configured.self_consistent_kzz
    assert rad.kzz_mixl_scale == pytest.approx(configured.kzz_mixl_scale)
    assert reloaded.adiabat == coupled_carma.adiabat
    assert rad.gap_max == configured.gap_max


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
    """The temperature profile moves away from the input profile."""
    from carmapy.results import Results

    coupled_carma.run(suppress_output=True)
    r = Results(coupled_carma)

    assert not np.allclose(r.T_history[:, 0], r.T_history[:, -1])
    # nzone counts the mask's contiguous runs, so it has to agree with the mask
    # it is derived from at every output.
    for it in range(r.nzone.size):
        assert len(r.convective_zones(it)) == r.nzone[it]


def test_tabulated_adiabat_end_to_end(coupled_carma):
    """A coupled run with ``adiabat="table"`` puts its convective interior on
    the tabulated adiabat, not the analytic one."""
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

    # The adjustment stops after a fixed number of sweeps, so the tolerance
    # is its residual (a few parts in 1e6), not machine precision.
    np.testing.assert_allclose(ratio[inside], tab[inside], rtol=1e-4)

    # The two gradients must differ by well over that tolerance, or the check
    # above would also pass with the flag ignored.
    separation = np.abs(tab[inside] / fit[inside] - 1).max()
    assert separation > 2e-3, f"table and fit differ by only {separation:.2e}"


def test_cadence_is_reported(coupled_carma):
    from carmapy.results import Results

    coupled_carma.run(suppress_output=True)
    r = Results(coupled_carma)

    assert r.rad_interval.shape == r.F_TOA.shape
    assert np.all(r.rad_interval >= 1)
    assert np.all(r.rad_interval <= coupled_carma.radiation.gap_max)


def test_kzz_follows_the_profile(coupled_carma):
    """The recorded Kzz is the one ``rce_kzz`` recomputes each output, not the
    profile ``add_kzz`` supplied."""
    from carmapy.results import Results

    coupled_carma.run(suppress_output=True)
    r = Results(coupled_carma)

    nt = r.T_history.shape[1]
    assert r.kzz_history is not None
    assert r.kzz_history.shape == (coupled_carma.NZ + 1, nt)
    assert np.all(np.isfinite(r.kzz_history))
    assert np.all(r.kzz_history > 0)

    assert not np.allclose(r.kzz_history[:, -1], coupled_carma.kzz_levels)

    # Recomputed as the profile moves, not once.
    assert not np.allclose(r.kzz_history[:, 0], r.kzz_history[:, -1])


def test_static_kzz_holds_the_input_profile(tmp_path, example_levels, ck_table):
    """``self_consistent_kzz=False`` keeps ``add_kzz``'s profile at every
    output."""
    from carmapy.carmapy import Carma
    from carmapy.chemistry import populate_abundances_at_cloud_base
    from carmapy.results import Results

    P, T, kzz, mu = example_levels
    c = Carma(str(tmp_path / "static_kzz"))
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
    c.set_radiation(t_int=T_INT, ck_table_path=ck_table, accel=20.0,
                    self_consistent_kzz=False)
    c.run(suppress_output=True)

    r = Results(c)
    for it in range(r.kzz_history.shape[1]):
        np.testing.assert_allclose(r.kzz_history[:, it], c.kzz_levels,
                                   rtol=1e-12)


def test_self_consistent_kzz_survives_load_carma(coupled_carma):
    from carmapy.base import load_carma

    coupled_carma.run(suppress_output=True)
    reloaded = load_carma(coupled_carma.name, restart=1)
    assert reloaded.radiation.self_consistent_kzz is True
    assert reloaded.radiation.kzz_mixl_scale == pytest.approx(1.0)


def test_kzz_mixl_scale_reaches_the_fortran(tmp_path, example_levels, ck_table):
    """At the first output, before the columns diverge, scaling
    ``kzz_mixl_scale`` by ``f`` scales the recorded Kzz by ``f**(4/3)``."""
    from carmapy.base import load_carma
    from carmapy.carmapy import Carma
    from carmapy.chemistry import populate_abundances_at_cloud_base
    from carmapy.results import Results

    P, T, kzz, mu = example_levels
    f = 2.0

    def run(scale, tag):
        c = Carma(str(tmp_path / tag))
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
        c.set_radiation(t_int=T_INT, ck_table_path=ck_table, accel=20.0,
                        kzz_mixl_scale=scale)
        c.run(suppress_output=True)
        return Results(c)

    base = run(1.0, "mixl_1")
    scaled = run(f, "mixl_2")

    np.testing.assert_allclose(scaled.kzz_history[:, 0] / base.kzz_history[:, 0],
                               f ** (4 / 3), rtol=1e-9)
    assert load_carma(str(tmp_path / "mixl_2")).radiation.kzz_mixl_scale == pytest.approx(f)


def test_final_profiles_track_the_coupled_run(coupled_carma):
    """final_T_levels and final_z_levels return the coupled run's evolved
    profile, not the inputs."""
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
    """An uncoupled run writes no temperature output."""
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
