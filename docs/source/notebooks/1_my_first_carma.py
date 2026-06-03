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

# %% {"nbsphinx": "hidden"}
import os
import carmapy as _cm
if os.environ.get("RUN_CARMA") != "True":
    _orig = _cm.Carma.run
    def _rtd_run(self, *a, suppress_output=False, **kw):
        if suppress_output:
            return
        return _orig(self, *a, suppress_output=True, **kw)
    _cm.Carma.run = _rtd_run

# %% [markdown]
# # My First Carma Object Model
#
# We are going to start by importing carmapy and creating a carma object.  Each
# carma object has a name that corresponds with a path to a directory that will be
# created upon running the model to house the model and its output—we will name
# this carma object "my_first_carma"

# %%
import carmapy

carma = carmapy.Carma("my_first_carma")

# %% [markdown]
# ## Setting Time Steps
#
# We will now tell carmapy its timestep, along with how often to output data, and
# how many timesteps to run

# %%
carma.set_stepping(dt=100, output_gap=1000, n_tstep=24000)

# %% [markdown]
# The parameters tell carmapy the following
#  - `dt` is the time step, in seconds
#  - `output_gap` is the number of time steps between each time the CARMA
#    simulation writes to file
#  - `n_tstep` is the total number of timesteps to run CARMA for

# %% [markdown]
# ## Setting Atmospheric Structure and Physical Parameters
#
# For this example we will use the example atmospheric structure that is included
# in carmapy but you of course can provide your own.  This atmospheric structure
# comes from a sonora diamondback
# ([Morley et al. 2024](https://doi.org/10.48550/arXiv.2402.00758)) run of a
# 2000 K solar-metallicity brown dwarf with a surface gravity of 31,600 cm/s²
# ($f_\text{sed} = 4$)

# %%
P_levs, T_levs, kzz_levs, mu_levs = carmapy.example.example_levels()

carma.add_P(P_levs)
carma.add_T(T_levs)
carma.add_kzz(kzz_levs)

# %% [markdown]
# Note that all of the inputs to carmapy are in cgs units so the units of
# `P_levs` is baryes, `T_levs` is in K, and `kzz_levs` is in cm²/s.  Each of
# these arrays are "levels" arrays which mean they describe the atmosphere at the
# boundary between layers—carmapy will model the formation of clouds at the
# "centers" location, between each of the provided "levels" locations.  Thus, if
# you want $N$ locations vertically where clouds might form, your "levels" arrays
# must be $N+1$ elements long.
#
# Also note that the first element of each array corresponds to the bottom of the
# atmosphere and the last element in each array corresponds to the top of the
# atmosphere.
#
# We also need to specify the surface gravity and the mean molecular weight of the
# atmosphere.  While our mean molecular weight varies throughout the atmosphere,
# carmapy requires a single value of the mean molecular weight to represent the
# entire atmosphere—in this case we will use the mean molecular weight at the
# bottom of the atmosphere

# %%
carma.set_physical_params(surface_grav=31600,
                            wt_mol=mu_levs[0])


# %% [markdown]
# Here we set the surface gravity to 31,600 cm/s² using the `surface_grav`
# argument and the mean molecular weight to the mean molecular weight at the
# bottom of our atmosphere using the `wt_mol` argument.
#
# We also need to tell CARMApy that the atmosphere is dominated by hydrogen so it
# can set the viscosity, thermal conductivity, and specific heat of the atmosphere
# accordingly

# %%
carma.set_atmospheric_parameters_from_defaults("Pure H2")

# %% [markdown]
# Note that you could instead set the viscosity, thermal conductivity, and
# specific heat manually with `carma.set_atmospheric_parameters()`
#
# Since we do not have altitude levels, we can tell carmapy to calculate the
# altitudes for us.

# %%
carma.calculate_z()

# %% [markdown]
# Since we have them, we can instead optionally use the complete mean molecular
# weight profile of our atmosphere for added accuracy:

# %%
carma.calculate_z(mu_levs)

# %% [markdown]
# For more accurate cloud calculations, sometimes we need deeper atmospheres than
# those commonly provided by evolutionary codes.  In those cases, we can tell
# carmapy to extend the bottom of the atmosphere adiabatically to a deeper
# pressure.

# %%
carma.extend_atmosphere(1e9)

# %% [markdown]
# Above we extended the atmosphere to a pressure of $10^{9}$ baryes.

# %% [markdown]
# ## Setting Cloud Physics
#
# carmapy requires that we tell it which species we should let form.  We do so by
# adding groups—types of particles that can form.  There are two types of groups:
#
# 1. Homogeneous Groups: These groups are particles which form via homogeneous
#    nucleation and thus do not need a seed particle to form
# 2. Heterogeneous Groups: These groups are particles which consist of one species
#    which heterogeneously nucleated onto another species, the latter of which
#    serves as a seed particle.
#
# We can add a homogeneous group and then a heterogeneous group that nucleates on
# the homogeneous group as follows:

# %%
carma.add_hom_group("TiO2", 1e-8)
carma.add_het_group("Mg2SiO4", "TiO2", 1e-8 * 2**(1/3))

# %% [markdown]
# `carma.add_hom_group("TiO2", 1e-8)` tells carmapy to allow TiO₂ to
# homogeneously condense with a minimum radius of $10^{-8}$ cm using the physical
# properties of TiO₂ preprogrammed into carmapy (some of these properties can be
# changed: see the documentation for Elements and Gases once they are written).
#
# `carma.add_het_group("Mg2SiO4", "TiO2", 1e-8 * 2**(1/3))` similarly tells
# carmapy to allow Mg₂SiO₄ to heterogeneously condense upon the TiO₂ seed
# particles in the group above.  This group has a minimum particle radius of
# $1.41 \times 10^{-8}$ cm.  Note that a homogeneous group must be created before
# a heterogeneous group can be defined to condense upon it.
#
# Any number of homogeneous and heterogeneous groups can be created, we created
# only 2 here for simplicity's sake.
#
# carmapy has a number of included default species. You can see the supported 
# species as with the following function:

# %%
carmapy.available_species()

# %% [markdown]
# If you wish to add a non-default species, see
# [Tutorial 4: Custom Condensates](4_custom_condensates.py).
#
# Additionally, the interaction term, `mucos`, for heterogeneously nucleating
# particles is only preprogrammed for certain combinations of species.  To see
# which species a given species is able to condense upon without manually
# specifying `mucos` you can run the following code:

# %%
carmapy.included_mucos("Mg2SiO4")

# %% [markdown]
# This indicates that Mg₂SiO₄ has an included interaction term for when it
# heterogeneously nucleates onto TiO₂.  If this term was not included in carmapy
# or we wanted to explore how the model changes if this interaction term is weaker
# we could, instead of running the above code to add groups, uncomment and run the
# below code:

# %%
# carma.add_hom_group("TiO2", 1e-8)
# carma.add_het_group("Mg2SiO4", "TiO2", 1e-8 * 2**(1/3), mucos=0.9)

# %% [markdown]
# ## Atmospheric Composition

# %% [markdown]
# We finally need to tell carma the gas abundances for each gas reservoir that is
# limiting for each species.  (The limiting gas for each species has defaults
# hard-coded—see the defaults documentation page once it is written).  The most
# common way to do that is to calculate the gas abundance at the base of each
# cloud—the place where the saturation vapor pressure for the species reaches the
# atmospheric pressure—and consider that the gas mixing ratio of the species at the
# bottom of the atmosphere.  Included in carmapy is a helper function which uses
# [fastchem](https://github.com/newstrangeworlds/fastchem) to calculate these gas
# phase abundances.

# %%
carmapy.chemistry.populate_abundances_at_cloud_base(carma)

# %% [markdown]
# If instead you want to set the gas abundances manually, see
# [Carma.set_nmr()](https://carmapy.readthedocs.io/en/latest/_autosummary/carmapy.Carma.html#carmapy.Carma.set_nmr)

# %% [markdown]
# ## Running the Model

# %% [markdown]
# To run the model all you need to do now is call `carma.run()`.  Note that the
# simulation should take 15-25 minutes to run, depending on your computer, and
# will give you multiple warnings/errors about negative gas concentrations and
# conditions at the beginning of your timestep—this is normal and expected.

# %%
carma.run(suppress_output=True)

# %% [markdown]
# Note that when you call `carma.run()` a directory with the name
# `"my_first_carma"` (or whatever you set the name of your carma object as)
# should have been created.  The directory stores the state of the model along
# with the model output—do not move or delete it while the model is running.  The
# next tutorial will go over how to read the output from that directory.
