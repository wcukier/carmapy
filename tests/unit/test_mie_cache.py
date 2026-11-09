"""The shared Mie-table cache in gen_band_mie_table.  The Mie calculation and
the ck table are replaced by stand-ins, so neither PyMieScatt work nor a ck
file is needed.
"""
import filecmp

import numpy as np
import pytest

from carmapy import mie, radiation


LAMBDAS = np.array([3e-3, 1e-3, 3e-4, 1e-4])


@pytest.fixture
def fake_mie(tmp_path, monkeypatch):
    """Point the cache at tmp_path, fix the band set, and count Mie runs."""
    monkeypatch.setenv("CARMAPY_MIE_CACHE", str(tmp_path / "cache"))
    monkeypatch.setattr(radiation, "band_wavelengths", lambda path: LAMBDAS)
    calls = []

    def gen_mie_table(species, min_r, ratio, n_r, lambdas, indices_path,
                      out_path, clamp_indices):
        calls.append(species)
        with open(out_path, "w") as f:
            f.write("r[cm]\tλ[cm]\tQ_ext\tQ_sca\tg\n")
            for ir in range(n_r):
                for lam in lambdas:
                    f.write(f"{min_r * ratio ** ir}\t{lam}\t2.0\t1.0\t0.5\n")
        return out_path

    monkeypatch.setattr(mie, "gen_mie_table", gen_mie_table)
    return calls


def _table(tmp_path, run, rmin=1e-8, **kwargs):
    out = tmp_path / run / "mie.dat"
    out.parent.mkdir(exist_ok=True)
    table = radiation.gen_band_mie_table("Mg2SiO4 on TiO2", rmin, 6,
                                         "unused.ck", str(out), **kwargs)
    return out, table


@pytest.mark.unit
def test_new_run_copies_the_cached_table(tmp_path, fake_mie):
    out1, t1 = _table(tmp_path, "run1")
    out2, t2 = _table(tmp_path, "run2")
    assert fake_mie == ["Mg2SiO4"]
    assert filecmp.cmp(out1, out2, shallow=False)
    np.testing.assert_array_equal(t1, t2)
    assert len(list((tmp_path / "cache").iterdir())) == 1


@pytest.mark.unit
def test_different_bins_miss_the_cache(tmp_path, fake_mie):
    _table(tmp_path, "run1")
    _table(tmp_path, "run2", rmin=2e-8)
    assert len(fake_mie) == 2


@pytest.mark.unit
def test_overwrite_recomputes(tmp_path, fake_mie):
    _table(tmp_path, "run1")
    _table(tmp_path, "run2", overwrite=True)
    assert len(fake_mie) == 2


@pytest.mark.unit
def test_cache_can_be_disabled(tmp_path, fake_mie):
    _table(tmp_path, "run1", use_cache=False)
    _table(tmp_path, "run2", use_cache=False)
    assert len(fake_mie) == 2
    assert not (tmp_path / "cache").exists()
