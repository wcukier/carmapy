import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os
import pickle
import dill
from types import SimpleNamespace
from carmapy.constants import gas_dict

# petroff 10 color cycle
petroff10 = ["#3f90da", "#ffa90e", "#bd1f01", "#94a4a2", "#832db6",
              "#a96b59", "#e76300", "#b9ac70", "#717581", "#92dadd"]


def save_results_pkl(results, file_path=None):
    """Save a CARMApy Results object to a pickle file.

    Parameters
    ----------
    results : carmapy.Results
        Results object to save.
    file_path : str, optional
        Output pickle path. If not given, writes to <results.path>/results.pkl.
    """
    # Build default output path from results object.
    if file_path is None:
        run_name = os.path.basename(results.path)
        file_path = os.path.join(results.path, f"results_obj_{run_name}.pkl")

    # Write the full Results object.
    with open(file_path, "wb") as handle:
        dill.dump(results, handle)

    return file_path


def save_carma_pkl(carma, file_path=None):
    """Save a CARMA object to a pickle file.

    Parameters
    ----------
    carma : carmapy.Carma
        CARMA object to save.
    file_path : str, optional
        Output pickle path. If not given, writes to <carma.path>/carma.pkl.
    """
    # Build default output path from the CARMA object.
    if file_path is None:
        file_path = os.path.join('./', f"carma_obj_{carma.name}.pkl")

    # Write the full CARMA object.
    with open(file_path, "wb") as handle:
        dill.dump(carma, handle)

    return file_path


def fake_carma(path, is_2d=False, nlongitude=1):
    """Build a minimal CARMA-like object from an output directory.

    Parameters
    ----------
    path : str
        Path to a CARMA output directory containing <dirname>.txt.
    is_2d : bool, optional
        Whether the run should be treated as 2D, by default False.
    nlongitude : int, optional
        Number of longitude bins for 2D runs, by default 1.
    """
    print('Oh no, you lost your carma? I\'ll rustle up a fake one quickly for you')
    path = os.path.abspath(path)
    path_end = os.path.basename(path)
    file_path = os.path.join(path, f"{path_end}.txt")

    with open(file_path) as handle:
        NZ, NGROUP, NELEM, NBIN, NGAS, nstep, iskip = np.array(
            handle.readline().split(), dtype=int
        )

    def _mk_group(i):
        elem = SimpleNamespace(ielem=min(i + 1, NELEM))
        return SimpleNamespace(igroup=i + 1, mantle=elem, core=elem)

    print('Here you go!')

    return SimpleNamespace(
        name=path,
        NZ=NZ,
        NBIN=NBIN,
        groups={f"group_{i+1}": _mk_group(i) for i in range(NGROUP)},
        elems={f"elem_{i+1}": SimpleNamespace(ielem=i + 1) for i in range(NELEM)},
        gases={f"gas_{i+1}": None for i in range(NGAS)},
        n_tstep=nstep - 1,
        output_gap=iskip,
        dt=1.0,
        is_2d=is_2d,
        NLONGITUDE=nlongitude,
        _citation={},
        log_metallicity=0.0,
    )

def plot_condensation_curves(self, ax=None, skip_gases=None, lon_idxs=None, **kwargs):
    """Plots condensation curves for each gas species in the CARMA model
    overlaid on the model's P-T profile.

    For each gas the condensation curve is the locus of (T, P) points at
    which the gas partial pressure equals its saturation vapour pressure
    (i.e. saturation ratio = 1).  The curve is computed analytically from
    the vapour-pressure coefficients stored in ``constants.gas_dict`` using
    the deep-atmosphere gas mixing ratio from the simulation.  Species whose
    vapour-pressure formula does not have a temperature coefficient (e.g. H2O,
    which uses a custom Fortran routine) are plotted using the saturation
    vapour pressures stored directly in ``results.sat_vp``.

    Parameters
    ----------
    ax : matplotlib.axes.Axes, optional
        Axes to draw on.  A new figure is created when *None*.
    skip_gases : list[int], optional
        List of gas indices to exclude.  Defaults to ``[]`` (plot all
        species).
    lon_idxs : list[int], optional
        Longitude indices to average over for the reference P-T profile.
        Defaults to ``None`` (average all longitudes).
    **kwargs
        Extra keyword arguments forwarded to ``ax.plot`` for the P-T
        profile line only.

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax : matplotlib.axes.Axes
    """
    if skip_gases is None:
        skip_gases = []

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()

    ax.set_prop_cycle(mpl.cycler(color=petroff10))

    P = self.P_centers      # [barye], shape (NZ,)
    T = self.T_centers      # [K], shape (NZ,) or (NZ, NLONG)
    log_met = getattr(self, "log_metallicity", 0.0)

    if T.ndim > 1:
        cmap = plt.get_cmap("viridis")
        for i in range(T.shape[1]):
            print("Plotting longitude index", i)
            color = cmap(i / max(T.shape[1] - 1, 1))
            ax.plot(T[:, i], P / 1e6, color=color, lw=0.5, alpha=0.5, **kwargs)
        if lon_idxs is None:
            T = np.mean(T, axis=1)
        else:
            T = np.mean(T[:, lon_idxs], axis=1)
        ax.scatter(T, P / 1e6, color="black", label="Average", **kwargs, s=3)#, lw=3)
    else:
        ax.plot(T, P / 1e6, color="black", lw=2, label="P-T profile", **kwargs)

    # Pressure grid for computing smooth condensation curves
    P_curve = np.logspace(np.log10(P.min()), np.log10(P.max()), 300)

    seen_materials = set()
    for group_obj in self.groups.values():
        material = group_obj.material
        if material in seen_materials:
            continue
        seen_materials.add(material)

        gas_obj = group_obj.gas
        if material in skip_gases or gas_obj.name in skip_gases:
            continue

        # Deep abundance: max over altitude gives uninhibited mixing ratio
        x_gas = np.max(gas_obj.nmr)

        if x_gas <= 0:
            continue

        tcoeff = group_obj.vp_tcoeff

        if tcoeff == 0:
            # No analytic inversion available (e.g. H2O Murphy 2005).
            # Cannot plot without sat_vp from simulation results; skip.
            continue
        else:
            offset    = group_obj.vp_offset
            metcoeff  = group_obj.vp_metcoeff
            logpcoeff = group_obj.vp_logpcoeff

            # Condensation condition: sat_vp(T, P) = x_gas * P
            #   1e6 * 10^(offset - tcoeff/T - metcoeff*log_met
            #              - logpcoeff*log10(P*1e-6)) = x_gas * P
            # Solving for T:
            #   T = tcoeff / (6 + offset - metcoeff*log_met
            #                 - logpcoeff*log10(P*1e-6)
            #                 - log10(x_gas * P))
            denom = (
                6.0
                + offset
                - metcoeff  * log_met
                - logpcoeff * np.log10(P_curve * 1e-6)
                - np.log10(x_gas * P_curve)
            )
            with np.errstate(invalid="ignore", divide="ignore"):
                T_cond = np.where(denom > 0, tcoeff / denom, np.nan)

            ax.plot(T_cond, P_curve / 1e6, label=material)

    ax.set_yscale("log")
    ax.invert_yaxis()
    ax.set_xlabel("Temperature [K]")
    ax.set_ylabel("Pressure [bar]")
    ax.legend(bbox_to_anchor=(1.0, 1.0), loc="upper left")
    fig.tight_layout()

    return fig, ax