[![Documentation Status](https://readthedocs.org/projects/carmapy/badge/?version=latest)](https://carmapy.readthedocs.io/en/latest/?badge=latest)

# CARMApy

![https://github.com/wcukier/carmapy/blob/main/docs/source/logo.png](https://github.com/wcukier/carmapy/blob/main/docs/source/logo.png)


CARMApy is a python wrapper of the Community Aerosol and Radiation Model for Atmospheres (CARMA), originally developed by Turco et al. (1979) and Toon et al. (1979).  CARMApy specifically wraps ExoCARMA, a version of CARMA developed by Gao et al. (2018) and Powell et al. (2018) to model clouds on exoplanets.  This code is still under development: while we believe there are no significant errors in the code, the CARMApy wrapper is still young.  If you encounter any issues or bugs, please open an issue or email Wolf Cukier (wcukier@uchicago.edu)

The documentation, including detailed tutorials and installation instructions, are available here: [https://carmapy.readthedocs.io](https://carmapy.readthedocs.io)


## Features

- Detailed microphysical modeling including 
    - Homogeneous and heterogeneous nucleation
    - Particle growth and evaporation
    - Vertical transport including eddy diffusion, falling, and bulk winds
    - Coagulation
- The ability to use user-defined condensates by specifying their physical properties
- An opinionated API that allows for easy learning of the code
- Easy reading of output data and plotting of results
- Helper functions with reasonable physical assumptions to supplement input data 