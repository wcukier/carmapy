import pytest
import warnings
from carmapy.carmapy import Carma
from carmapy.constants import JUPITER_RADIUS


# ── set_stepping ──────────────────────────────────────────────────────────────

@pytest.mark.unit
def test_set_stepping_sets_dt():
    c = Carma("_t")
    c.set_stepping(dt=100)
    assert c.dt == 100


@pytest.mark.unit
def test_set_stepping_dt_float_raises():
    c = Carma("_t")
    with pytest.raises(TypeError):
        c.set_stepping(dt=100.5)


@pytest.mark.unit
def test_set_stepping_dt_negative_raises():
    c = Carma("_t")
    with pytest.raises(ValueError):
        c.set_stepping(dt=-1)


@pytest.mark.unit
def test_set_stepping_sets_output_gap():
    c = Carma("_t")
    c.set_stepping(dt=100, output_gap=50)
    assert c.output_gap == 50


@pytest.mark.unit
def test_set_stepping_sets_n_tstep():
    c = Carma("_t")
    c.set_stepping(n_tstep=5000)
    assert c.n_tstep == 5000


@pytest.mark.unit
def test_set_stepping_partial_update_preserves_output_gap():
    c = Carma("_t")
    c.output_gap = 999
    c.set_stepping(dt=100)
    assert c.output_gap == 999


@pytest.mark.unit
def test_set_stepping_partial_update_preserves_n_tstep():
    c = Carma("_t")
    c.n_tstep = 8880
    c.set_stepping(dt=100)
    assert c.n_tstep == 8880


@pytest.mark.unit
def test_set_stepping_output_gap_negative_raises():
    c = Carma("_t")
    with pytest.raises(ValueError):
        c.set_stepping(dt=100, output_gap=-1)


# ── set_physical_params ───────────────────────────────────────────────────────

@pytest.mark.unit
def test_set_physical_params_surface_grav():
    c = Carma("_t")
    c.set_physical_params(surface_grav=31600.0)
    assert c.surface_grav == 31600.0


@pytest.mark.unit
def test_surface_grav_negative_raises():
    c = Carma("_t")
    with pytest.raises(ValueError):
        c.set_physical_params(surface_grav=-980.0)


@pytest.mark.unit
def test_set_wt_mol():
    c = Carma("_t")
    c.set_physical_params(wt_mol=2.3)
    assert c.wt_mol == pytest.approx(2.3)


@pytest.mark.unit
def test_wt_mol_negative_raises():
    c = Carma("_t")
    with pytest.raises(ValueError):
        c.set_physical_params(wt_mol=-2.0)


@pytest.mark.unit
def test_r_planet_jovian_radius_conversion():
    c = Carma("_t")
    c.set_physical_params(r_planet=1.0, use_jovian_radius=True)
    assert abs(c.r_planet - JUPITER_RADIUS) < 1e6


@pytest.mark.unit
def test_r_planet_small_auto_assumes_jovian():
    c = Carma("_t")
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        c.set_physical_params(r_planet=1.5)
    assert len(w) == 1
    assert "jovian" in str(w[0].message).lower() or "jovian radius" in str(w[0].message).lower()
    assert abs(c.r_planet - 1.5 * JUPITER_RADIUS) < 1e6


@pytest.mark.unit
def test_r_planet_too_large_raises():
    c = Carma("_t")
    with pytest.raises(ValueError):
        c.set_physical_params(r_planet=2000, use_jovian_radius=True)


@pytest.mark.unit
def test_r_planet_negative_raises():
    c = Carma("_t")
    with pytest.raises(ValueError):
        c.set_physical_params(r_planet=-1.0)


@pytest.mark.unit
def test_set_physical_params_partial_update():
    c = Carma("_t")
    c.set_physical_params(surface_grav=980.0)
    c.set_physical_params(wt_mol=2.3)
    assert c.surface_grav == 980.0
    assert c.wt_mol == pytest.approx(2.3)


# ── set_cloud_boundary_type ───────────────────────────────────────────────────

@pytest.mark.unit
@pytest.mark.parametrize("bc", ["fixed_conc", "fixed_flux", "zero_grad"])
def test_valid_cloud_top_boundary(bc):
    c = Carma("_t")
    c.set_cloud_boundary_type(bc, "fixed_conc")


@pytest.mark.unit
@pytest.mark.parametrize("bc", ["fixed_conc", "fixed_flux", "zero_grad"])
def test_valid_cloud_bot_boundary(bc):
    c = Carma("_t")
    c.set_cloud_boundary_type("fixed_conc", bc)


@pytest.mark.unit
def test_invalid_cloud_top_boundary_raises():
    c = Carma("_t")
    with pytest.raises(ValueError):
        c.set_cloud_boundary_type("bogus", "fixed_conc")


@pytest.mark.unit
def test_invalid_cloud_bot_boundary_raises():
    c = Carma("_t")
    with pytest.raises(ValueError):
        c.set_cloud_boundary_type("fixed_conc", "bogus")


@pytest.mark.unit
def test_cloud_boundary_type_stored():
    c = Carma("_t")
    c.set_cloud_boundary_type("fixed_flux", "zero_grad")
    assert c.top_bound_type_cloud == "fixed_flux"
    assert c.bot_bound_type_cloud == "zero_grad"


# ── set_gas_boundary_type ─────────────────────────────────────────────────────

@pytest.mark.unit
def test_set_gas_boundary_type_valid():
    c = Carma("_t")
    c.set_gas_boundary_type("fixed_flux", "fixed_conc")
