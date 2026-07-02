---
title: '`CARMApy`: A User-Friendly Open-Source Cloud Microphysics Code in Python'
tags:
  - Python
  - astronomy
  - exoplanets
  - atmospheres
authors:
  - name: Wolf Cukier
    orcid: 0000-0002-8658-3811
    affiliation: 1
  - name: Diana Powell
    orcid: 0000-0002-4250-0957
    affiliation: 1
  - name: Xi Zhang
    orcid: 0000-0002-8706-6963
    affiliation: 2
  - name: Peter Gao
    orcid: 0000-0002-8518-9601
    affiliation: 3
  - name: Dominic Samra
    orcid: 0000-0002-8956-2047
    affiliation: 1
  - name: Vighnesh Nagpal
    orcid: 0000-0001-5909-4433
    affiliation: 1
affiliations:
 - name: Department of Astronomy & Astrophysics, the University of Chicago, Chicago, IL, 60637, USA
   index: 1
   ror: 024mw5h28
 - name: Department of Earth and Planetary Sciences, University of California Santa Cruz, Santa Cruz, CA, 95064, USA
   index: 2
   ror: 03s65by71
 - name: Earth and Planets Laboratory, Carnegie Institution for Science, Washington, DC, USA
   index: 3
   ror: 04jr01610

date: 11 June 2026?
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: ??? <- update this with the DOI from AAS once you know it.
aas-journal: The Astrophysical Journal <- The name of the AAS journal.
---
# Summary

Clouds are poorly understood yet incredibly important for understanding our
observations of exoplanets.  The complexity of clouds is due to the
fact that macroscopic observables of clouds are dictated by the 
properties, and in particular size-distributions, of the tiny particles that
make up the cloud.  These properties are in turn determined by the physics
that happens at the molecular scale.  The `carmapy` package allows for the 
easy simulation of these microscopic processes and interpretation of these
results for comparisons with observable properties.

# Statement of need

`CARMApy` sits towards the high complexity end of the efficiency-complexity
trade-off space.  As such, it is best suited for tasks where getting an accurate result
or a result with minimal physical assumptions is desired over getting a result 
fast.  `CARMApy` is best suited for tasks such as:

 - Forward modeling the aerosol properties of individual exoplanets
 - Post-processing General Circulation Model (GCM) outputs to include more 
   accurate cloud treatments
 - Running a grid of models to understand demographic trends in exoplanets and/or
   brown dwarfs
 - Modeling morning/evening limb asymmetries caused by cloud formation and 
   transport
 - Validating the result of a retrieval with a detailed forward cloud model
 - Understanding which microphysical processes are important to cloud formation
   in different atmospheric regimes

`CARMApy` builds upon the previous papers that have used the underlying `CARMA`
model as applied to exoplanets [e.g. @Gao2018MicrophysicsKClZnS; @Powell2018FormationSilicateTitanium; @Powell2024TwodimensionalModelsMicrophysical].  This work generalizes the work done by these papers
while also open-sourcing and documenting the code.  The code is designed such
that it can be used by non-experts in microphysical cloud modeling, such as a
motivated undergraduate student or an observer validating a retrieval model.

# State of the field

The community uses a variety of tools and techniques across the 
efficiency-complexity trade-off spectrum to model clouds on exoplanets.  On the 
extreme efficiency end of the scale, there are techniques such as gray cloud
models that assume constant opacity across all wavelengths.  These models are
thus incredibly fast to run and are often incorporated into retrieval models
[e.g. @Xue2024JWSTTransmissionSpectroscopy].  Slightly more complex on this scale
are equilibrium cloud condensation
models [@Ackerman2001PrecipitatingCondensationClouds].  These models incorporate some physics, assuming
that clouds only form in regions of the atmosphere which are supersaturated.  
These models require a tunable parameter, $f_\text{sed}$, that accounts for the
removal of particles from the atmosphere via sedimentation.  This framework has
been incorporated into open source python models such as `Virga` [@Batalha2025CondensationCloudsSubstellar].
A bit more complex, and thus slower to run, than those models are the 
moment-method kinetic cloud models which can incorporate detailed microphysics
and chemistry calculations, but require a prescribed fuctional form for the 
particle size distributions.  Open-source versions of these models such as
`mini-cloud` [e.g. @Lee2025MonoculturePolydisperseMoment] exist, but those models are not designed and 
packaged for widescale use.

One more step up this complexity ladder, which is where `CARMApy` sits, are 
bin-scheme microphysical models.  These models are able to incorporate detailed
microphysical modeling while resolving arbitrary particle size distributions.
The Community Aerosol and Radiation Model for Atmospheres, `CARMA` 
[@Turco1979OneDimensionalModelDescribing; @Toon1988MultidimensionalModelAerosols; @Bardeen2008NumericalSimulationsThreedimensional], originally designed for Earth clouds,
is the most widely used of this class of models, and is open source.  There have
been a number of forks of the `CARMA` code to apply the code to non-Earth planets
such as `PlanetCARMA` [@Barth2020PlanetCARMANewFramework], which is currently non-public, and 
`ExoCAM-CARMA`, which only supports hazes and not clouds and while open-source 
is currently minimally documented.  `CARMApy` builds on the `ExoCARMA` 
[@Powell2018FormationSilicateTitanium; @Gao2018SedimentationEfficiencyCondensation] version of `CARMA` which has extended the 
functionality of `CARMA` to include condensates and background atmospheres not
otherwise modelable with `CARMA`.  Until the publication of this work, however,
the source code of `ExoCARMA` has not been public.  Additionally, all of these
codes are written in `Fortran` over a almost half-century long development 
history, leading to code that only a handful of expert users are able to 
effectively use.  `CARMApy` is designed to make the scientific power of the
`CARMA` family of codes accessible to be used by the general exoplanet community.


# Software design

The fundamental design philosophy of `CARMApy` is to expose the core functionality
of `CARMA` to the user in an easy to understand manner.  This means that we 
designed `CARMApy` to be opinionated, having default values and implementations
for the common ways of using the code (e.g. automatically allowing particles to 
grow or evaporate upon specification of a nucleation pathway 
that forms these particles).  While these defaults are strong, and often will
occur without requiring the user explicitly invoke them, we allow the user to 
override these defaults to maintain feature parity with the underlying `CARMA`
code.

Due to the reliance of `CARMA` on custom types and on `selected_real_kind` for
double precision calculations, we found it difficult to directly bind the `CARMA`
types and functions from fortran to python directly.  Instead, upon installation `CARMApy`
compiles our version `CARMA` into an executable that reads all relevant input
and configuration data from files which the python layer writes.  The 
`Carma.run()` function which runs the simulation will write these input files 
and run the executable.  Once the simulation finishes running, the user can
call `Carma.read_results()` to load the results written to file by the `CARMA`
simulation into a python `Results` object for easy analysis by the user.  This
translation between layers happens once in each direction per run and thus does
not significantly contribute to computational overhead.

A complete description of the physical processes and numerical methods is 
provided in the associated paper [@Cukier2026CARMApyOpenSourcePython]

# Research impact statement

`CARMApy` has a high potential for research impact.  The `ExoCARMA` code which
`CARMApy` upgrades and wraps has already been used in a large number of papers
[e.g., @Powell2018FormationSilicateTitanium; @Gao2018MicrophysicsKClZnS; @Mang2024MicrophysicalPrescriptionsParameterized; @Rooney2022NewSedimentationModel; @Wong2020OpticalNearinfraredTransmission], although the vast majority of these papers have relied on a small
handful of coauthors to run the `CARMA` model.  By bringing this code to python,
along with providing detailed documentation and tutorials, the number of 
researchers who would be able to run `CARMApy` is significantly greater.  
`CARMApy` is already being used in a number of in prep papers [Cukier et al. in prep; Kennedy et al. in prep; Steinrueck et al. in prep; Samra et al. in prep] with the users spanning
career stages from undergraduate to post doc, and these users providing 
feedback on documentation, reporting bugs, and suggesting new features.

We additionally have performed benchmarking our code to the versions of 
`ExoCARMA` previously used in publications.  For a sample brown dwarf case run
on a M4 Macbook Pro we find that our code runs 3.8 times faster than previous 
versions when taking advantage of the newly implemented multithreading 
capabilities of `CARMApy` and 1.9 times faster than the previous versions while
running in single threaded mode.  These speedups occur while reproducing 
the observables, in this case a thermal emission spectum, to 3 parts in 10,000,
a number which is consistent with floating point error caused by our 
optimizations.

# AI Usage Disclosure

AI tools (Sonnet 4.6, Opus 4.7, Opus 4.8) were used to help generate packaging
code, generate the test suite, refactor the codebase, tweak the documentation,
optimize the codebase, and implement helper functions.  All AI generated code 
has been read through and checked by a human, with extra scrutiny given to any 
code the that touches the core simulation logic.  Additionally, manual 
benchmarking was performed to ensure that any AI generated code did not 
unexpectedly change the physics of the simulation.

The subset of the above AI tools were used for this paper to format metadata and
citations, which were then checked manually, and provide suggestions on how to 
improve the prose, which were accepted, rejected, and implemented manually.


# Acknowledgements

The authors would like to thank Thomas Kennedy for helping beta test this code.  
This work benefited from the 2025 Exoplanet Summer Program in the Other Worlds 
Laboratory (OWL) at the University of California, Santa Cruz, a program funded 
by the Heising-Simons Foundation and NASA. X.Z. acknowledges support from the 
NSF grant (AST2307463), NASA Exoplanet Research grant (80NSSC22K0236), and the 
NASA Interdisciplinary Consortia for Astrobiology Research grant (80NSSC21K0597).
VN acknowledges support from the National Science Foundation Graduate Research 
Fellowship under Grant No. DGE 2140001. D.S. acknowledges funding as part of 
JWST GO program 3969 (PIs: Espinoza, Powell).

# References
