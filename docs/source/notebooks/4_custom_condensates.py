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
# # Custom Condensates
#
# CARMApy has a number of default condensates built in but also allows the user
# to define custom condensates.  As a reminder, you can check what default
# condensates are included in CARMApy by running the following command:

# %%
import carmapy

carmapy.available_species()

# %% [markdown]
# ## Relevant Properties
# For reasons related to how CARMA was originally designed, condensate properties
# are tied to the limiting gas reservoir the condensate is formed from.  Thus,
# defining a custom condensate means creating a custom `carmapy.Gas` object.
# This unfortunately has a side effect that each gas can only form condensates of
# a single chemical formula (ie a gas can form Pure Fe and Fe on TiO₂ condensates
# but cannot also form Fe₂SiO₄ on TiO₂ as the latter has a different chemical
# formula).
#
# A large number of properties need to be defined for a custom condensate; we
# will go through them below, but you can also see a reference page for the
# `carmapy.Gas` object
# [here](https://carmapy.readthedocs.io/en/latest/_autosummary/carmapy.Gas.html)
#
# For purposes of this tutorial, we will be recreating Al₂O₃
#
# To do so we need to know the following properties:
# - The density of the condensate
# - The molar mass of the condensate
# - The molar mass of the limiting gas species
# - The collisional diameter
# - The cosine of the contact angle between the gas species and any seed
#   particles it might be condensing on
# - The [Hill Notation](https://en.wikipedia.org/wiki/Chemical_formula#Hill_system)
#   formula for the gas species
# - The formula for the saturation vapor pressure of the condensate (explained
#   below)
# - The formula for the surface tension of the condensate (explained below)
# - The stoichiometry factor of the condensation reaction (explained below)
# - Whether or not the reaction is considered a type III reaction (explained
#   below)
# - (Optionally) the latent heat of evaporation for the condensate
#
#

# %% [markdown]
# ### Simple Properties
# - Al₂O₃ has a density of 3.99 g/cm³
# - Al₂O₃ has a molar mass of 101.926 g/mol
# - The limiting gas species, Al, has a molar mass of 26.98 g/mol
# - Al₂O₃'s collisional diameter is approximately $3.825\times10^{-8}$ cm
#   (value estimated from Dobrovinskaya et al. 2009)
# - The cosine of the contact angle between Al₂O₃ and TiO₂ is 0.724172
#   (value taken from P. Gao et al. 2020)
#
# ### Saturation Vapor Pressure
# The saturation vapor pressure of a species can be a function of both pressure
# and temperature.  carmapy assumes that the formula for vapor pressure can be
# parameterized as follows:
#
# $$
# \begin{align*}
# \log_{10} \frac{p'}{10^6 \text{ barye}} &=  \texttt{vp\_offset} \\
# &- \frac{\texttt{vp\_tcoeff}}{T} \\
# &- \texttt{vp\_metcoeff} \cdot [Fe/H] \\
# &- \texttt{vp\_logpcoeff} \cdot \log_{10} \frac{P}{10^6 \text{ barye}}
# \end{align*}$$
#
# where $p'$ is the saturation vapor pressure, $T$ is the temperature, and $P$
# is the atmospheric pressure
#
# Using the values from Wakeford et al. 2017 for Al₂O₃ we have:
# - `vp_offset` = 17.7
# - `vp_tcoeff` = 45892.6 K
# - `vp_metcoeff` = 1.66
# - `vp_logpcoeff` = 0
#
# ### Surface Tension
# CARMApy also parameterizes the surface tension of condensates as follows:
#
# $$\texttt{surface\_tension} = \texttt{surften\_0}  - \texttt{surften\_slope} * T$$
#
# where $T$ is the temperature in K.
#
# From Kozasa et al. 1989 we estimate that for Al₂O₃
# - `surften_0` = 690 dyne/cm
# - `surften_slope` = 0
#
# ### Stoichiometry Factor
# The stoichiometry factor is determined by the reaction that forms the
# condensate, ie for a reaction which is
# $$X (\text{limiting gas}) + (\text{other stuff}) \to Y (\text{condensate})
# + (\text{some other stuff})$$
#
# the stoichiometry factor is $X/Y$
#
# As we are considering Al as our limiting gas and it takes 2 monatomic Al to
# form 1 molecule Al₂O₃, the stoichiometry factor will be 2
#
# ### Type III reaction
#
# Condensates which form via type III reactions, as defined by Helling and
# Woitke 2012, have a modified supersaturation ratio.  Type III reactions are
# reactions that, besides the condensate, have 2 or more products or 2 or more
# reactants.
#
# Since the main chemical reaction which forms Al₂O₃ is
# $$2 \text{AlOH} + \text{H}_2\text{O} \to \text{Al}_2\text{O}_3 [\text{s}]
# + 2 \text{H}_2$$
#
# Al₂O₃ is a type III reaction
#
# ### Latent heat of evaporation
#
# You can choose to provide a latent heat of evaporation; if you don't, CARMApy
# will calculate one for you following the assumptions of Charnay et al. 2015
# (TODO: check this citation):
#
# $$\texttt{lat\_heat\_e} = \texttt{vp\_tcoeff} \cdot \ln(10) *
# \frac{R}{\texttt{wtmol\_dif}}$$
#
# where $R$ is the ideal gas constant
#
# ## Defining the Gas Object
# First we will set up the carma simulation as we did previously, without adding
# the Mg₂SiO₄

# %%
carma = carmapy.Carma("custom_gases")

carma.set_stepping(dt=100, output_gap=100, n_tstep=3000)


P_levs, T_levs, kzz_levs, mu_levs = carmapy.example.example_levels()

carma.add_P(P_levs)
carma.add_T(T_levs)
carma.add_kzz(kzz_levs)

carma.set_physical_params(surface_grav=31600,
                            wt_mol=mu_levs[0])
carma.set_atmospheric_parameters_from_defaults("Pure H2")

carma.calculate_z(mu_levs)
carma.extend_atmosphere(1e10)

carma.add_hom_group("TiO2", 1e-8)

# %% [markdown]
# Then to add the user-defined condensate, first we add a gas to our carma
# object with all the relevant properties, and then we add the nucleation
# pathways we want to consider:

# %%
carma.add_gas("Al_oxide",
    rho_cond = 3.99,
    wtmol = 101.926,
    wtmol_dif = 26.98,
    coldia = 3.825e-8,
    hill_formula = "Al",
    vp_offset = 17.7,
    vp_tcoeff = 45892.6,
    vp_metcoeff = 1.66,
    vp_logpcoeff = 0,
    surften_0 = 690,
    surften_slope = 0,
    stofact = 2,
    is_typeIII = True
)

carma.add_het_group("Al_oxide",
                    "TiO2",
                    1e-8 * 2**(1/3),
                    mucos=0.724172)

# %% [markdown]
# Note that we can add the nucleation pathway by using the name of the gas we
# defined.
#
# **Note: use a short name for the gas — long names (around 12 or more
# characters) have a chance of breaking CARMA**
#
# Then, as before, we can populate the gas abundances as follows:

# %%
carmapy.chemistry.populate_abundances_at_cloud_base(carma)

# %% [markdown]
# As with other carma objects, you can run the simulation by running the
# following cell:

# %%
# carma.run()
