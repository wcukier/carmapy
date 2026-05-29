# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: base
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Other Features
#
# CARMApy offers a number of other ways to customize your models. Below we will
# initialize a carma object for use in these examples

# %%
import carmapy
import numpy as np

carma = carmapy.example.example_carma("other_features")

# %% [markdown]
# ## Winds
#
# In addition to eddy diffusion and settling, CARMApy can take into account
# vertical advection of the atmosphere by specifying vertical winds.  For
# example, to specify an updraft of 100 cm/s at every level of the atmosphere
# we can do the following:

# %%
winds = np.ones(carma.NZ) * 100
carma.add_vertical_winds(winds)

# %% [markdown]
# Note that the length of the `winds` array is `carma.NZ` because winds are
# specified at the center of each atmosphere layer, not the boundary.

# %% [markdown]
# ## Boundary Conditions

# %% [markdown]
# CARMA requires that boundary conditions are specified for both the gas and the
# cloud particles at the top and bottom of the atmosphere—if no boundary
# conditions are specified CARMApy will supply CARMA with default boundary
# conditions (discussed below).  There are three types of boundary conditions
# allowed in CARMApy:
#
# - `"fixed_conc"`: The concentration of the species is fixed at the boundary
# - `"fixed_flux"`: The flux of the species is fixed at the boundary
# - `"zero_grad"`: The concentration gradient of the species is zero at the
#   boundary
#
# When selecting boundary condition types, note that the same type of boundary
# condition must be used for all gas species at a given boundary, and similarly
# for cloud species.
#
# The default boundary condition behavior is as follows:
# - At the bottom of the atmosphere all cloud particles are assumed to have a
#   concentration of 0 and the gas concentration is set to that specified in
#   `gas.nmr`
# - At the top of the atmosphere the flux of both gas and cloud particles is
#   assumed to be 0
#
# Additionally if the boundary condition type is changed, CARMApy will default
# to all fluxes being 0, all cloud particle concentrations being 0, and throwing
# an error for unspecified gas concentrations.
#
# As an example, to specify a cloud particle flux at the top of the atmosphere
# (to simulate haze generation or infalling dust for example) you can specify
# the boundary conditions as follows:

# %%
carma.set_cloud_boundary_type("fixed_conc", "fixed_flux")

top_flux = np.zeros(carma.NBIN)
top_flux[0] = 1000

carma.set_cloud_boundary("Pure TiO2", top_flux=top_flux)

# %% [markdown]
# Here we set there to be a flux of 1000 TiO2 particles/cm²/s at the top of
# the atmosphere.  Specifically we set all of these particles to be at the
# nucleation size (ie the smallest size bin).  By changing other elements of
# the `top_flux` array to be non-zero we could add a flux of larger particles.
#
# Note that non-zero boundary conditions are only supported for
# single-composition cloud species.  For example, the above code works because
# the `"Pure TiO2"` group consists only of TiO2, but a non-zero
# `"Mg2SiO4 on TiO2"` group would not work as it consists of both Mg2SiO4 and
# TiO2.
#
# **Warning**: using `carmapy.load_carma()` to load a simulation which has
# boundary conditions set will only load the boundary conditions actually used
# in the simulation.  The unused boundary conditions have undefined behavior
# upon `carmapy.load_carma()`
