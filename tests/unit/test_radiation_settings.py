"""The Radiation settings object: validation, its &radiation namelist round
trip, and how set_radiation builds it.  No Fortran binary runs.
"""
import os

import f90nml
import pytest

from carmapy import Carma, Radiation


SETTINGS = dict(t_irr=1500.0, t_star=5000.0, mu0=0.25, w_surf=0.1,
                cloud_rad=False, self_consistent_kzz=False,
                kzz_mixl_scale=0.5, accel=20.0, gap_max=7, dT_tol=0.5,
                dtau_tol=0.05)


@pytest.mark.unit
def test_namelist_round_trip(tmp_path):
    """to_nml -> f90nml file -> from_nml rebuilds every setting; f90nml
    lower-cases the keys on the way."""
    ck = str(tmp_path / "table.ck")
    rad = Radiation(ck, 400.0, **SETTINGS)

    path = str(tmp_path / "input.nml")
    f90nml.Namelist({"radiation": rad.to_nml()}).write(path)
    back = Radiation.from_nml(f90nml.read(path)["radiation"])

    assert vars(back) == vars(rad)


@pytest.mark.unit
def test_uncoupled_namelist_reads_as_none():
    assert Radiation.from_nml(Radiation.UNCOUPLED_NML) is None
    assert Radiation.from_nml({}) is None


@pytest.mark.unit
def test_uncoupled_namelist_carries_the_defaults():
    """Apart from switching radiation off, the uncoupled entries are a
    default coupling's."""
    on = Radiation("table.ck", 1000.0).to_nml()
    off = Radiation.UNCOUPLED_NML
    assert on.keys() == off.keys()
    differ = {k for k in on if on[k] != off[k]}
    assert differ == {"do_radiation", "ck_table_file", "t_int"}


@pytest.mark.unit
@pytest.mark.parametrize("bad", [
    dict(t_int=0.0),
    dict(t_irr=-1.0),
    dict(t_star=-1.0),
    dict(mu0=0.0),
    dict(mu0=1.5),
    dict(w_surf=1.0),
    dict(gap_max=0),
    dict(kzz_mixl_scale=0.0),
])
def test_out_of_range_settings_raise(bad):
    kwargs = {"t_int": 1000.0, **bad}
    with pytest.raises(ValueError):
        Radiation("table.ck", **kwargs)


@pytest.mark.unit
def test_missing_table_loads_but_cannot_be_set(tmp_path):
    """A run whose ck table has moved still rebuilds from its namelist, but
    set_radiation refuses a table that is not there."""
    missing = str(tmp_path / "moved.ck")
    nml = Radiation(missing, 1000.0).to_nml()
    assert Radiation.from_nml(nml).ck_table_path == missing

    with pytest.raises(FileNotFoundError):
        Carma("_t").set_radiation(t_int=1000.0, ck_table_path=missing)


@pytest.mark.unit
def test_set_radiation_stores_the_settings(tmp_path):
    ck = tmp_path / "table.ck"
    ck.touch()
    c = Carma("_t")
    assert c.radiation is None

    c.set_radiation(t_int=1000.0, ck_table_path=str(ck), accel=20.0,
                    rad_gap_max=7, rad_dT_tol=0.5, rad_dtau_tol=0.05)

    rad = c.radiation
    assert isinstance(rad, Radiation)
    assert (rad.accel, rad.gap_max, rad.dT_tol, rad.dtau_tol) == (20.0, 7, 0.5, 0.05)
    assert c.t_evolves is True
    assert c._citation["radiation"] and c._citation["kzz_mixing_length"]


@pytest.mark.unit
def test_rejected_settings_leave_the_run_uncoupled(tmp_path):
    ck = tmp_path / "table.ck"
    ck.touch()
    c = Carma("_t")
    with pytest.raises(ValueError):
        c.set_radiation(t_int=1000.0, ck_table_path=str(ck), mu0=2.0)
    assert c.radiation is None
    assert c.t_evolves is False
