[![Tests](https://github.com/wcukier/carmapy/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/wcukier/carmapy/actions/workflows/tests.yml)
[![Documentation Status](https://readthedocs.org/projects/carmapy/badge/?version=latest)](https://carmapy.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://img.shields.io/pypi/v/carmapy.svg)](https://pypi.org/project/carmapy/)
[![Python versions](https://img.shields.io/pypi/pyversions/carmapy.svg)](https://pypi.org/project/carmapy/)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)

# CARMApy

![CARMApy logo](https://raw.githubusercontent.com/wcukier/carmapy/main/docs/source/logo.png)


CARMApy is a Python wrapper of the Community Aerosol and Radiation Model for Atmospheres (CARMA), originally developed by Turco et al. (1979), Toon et al. (1988), and Bardeen et al. (2008). CARMApy specifically wraps ExoCARMA, a version of CARMA developed by Gao et al. (2018) and Powell et al. (2018) to model clouds on exoplanets.

This code is still under development: while we believe there are no significant errors in the code, the CARMApy wrapper is still young. If you encounter any issues or bugs, please open an issue or email Wolf Cukier (wcukier@uchicago.edu).

The documentation, including detailed tutorials and installation instructions, are available here: [https://carmapy.readthedocs.io](https://carmapy.readthedocs.io)


## Features

- Detailed microphysical modeling including 
    - Homogeneous and heterogeneous nucleation
    - Particle growth and evaporation
    - Vertical transport including eddy diffusion, falling, and bulk winds
    - Coagulation
- The ability to use user-defined condensates by specifying their physical properties
- A beginner-friendly, opinionated API that makes the code easy to learn
- Easy reading of output data and plotting of results
- Helper functions with reasonable physical assumptions to supplement input data


## Installation

CARMApy is available on PyPI:

```bash
pip install carmapy
```

CARMApy requires Python 3.10 or newer. Pre-built wheels bundle a standalone
Fortran binary, so no compiler is needed for a standard install. See the
[documentation](https://carmapy.readthedocs.io) for source builds and details.


## Quick Start

```python
import carmapy

# Create a simulation (the name becomes the output directory)
carma = carmapy.Carma("my_first_run")

# Set the timestep [s], output cadence, and number of steps
carma.set_stepping(dt=100, output_gap=100, n_tstep=24000)

# Load the example atmosphere (a 2000 K Sonora Diamondback profile).
# All inputs are cgs: P in barye, T in K, kzz in cm^2/s.
P_levs, T_levs, kzz_levs, mu_levs = carmapy.example.example_levels()
carma.add_P(P_levs)
carma.add_T(T_levs)
carma.add_kzz(kzz_levs)

# Surface gravity [cm/s^2] and mean molecular weight 
carma.set_physical_params(surface_grav=31600, wt_mol=mu_levs[0])
carma.set_atmospheric_parameters_from_defaults("Pure H2")  # H2-dominated
carma.calculate_z(mu_levs)        # derive altitudes from P/T/mu

# Add cloud species: TiO2 nucleates homogeneously, Mg2SiO4 grows on it.
# The second argument is the minimum particle radius [cm].
carma.add_hom_group("TiO2", 1e-8)
carma.add_het_group("Mg2SiO4", "TiO2", 1e-8 * 2**(1/3))

# Initialize gas abundances (via FastChem) and run
carmapy.chemistry.populate_abundances_at_cloud_base(carma)
carma.run()

results = carma.read_results()
```

See the [tutorials](https://carmapy.readthedocs.io) for complete, runnable examples.


## Citation

If you use CARMApy in your work, please cite CARMApy itself along with the
underlying models it is built on. At a minimum:

> This work uses CARMApy (Cukier et al. 2026), which is built on the
> ExoCARMA (Gao et al. 2018, Powell et al. 2018) version of the CARMA model
> (Turco et al. 1979, Toon et al. 1988, Bardeen et al. 2008).

If you use additional features such as the 2D mode or the FastChem chemistry
interface, please also cite the relevant papers. See the
[citation page](https://carmapy.readthedocs.io/en/latest/citation.html) for the
complete list of references and BibTeX entries.


## License

CARMApy is released under the [Apache License 2.0](LICENSE). 