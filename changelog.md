
### Version 1.0.1
- Fixed results reading such that the final timestep(s) should no longer be truncated
- Fixed the gas input such that the full gas input file, not just the bottom layer, gets passed to the model


### Version 1.0.2
- Added enstatite (MgSiO3) condensate
- Small fixes to coagulation
- Fixes to upwards wind velocity


### Unreleased
- Engine: rate-informed adaptive substepping. When a microphysics substep fails, the retry jumps straight to the substep count the failing rate calls for (gas depletion, gas-change threshold, or supersaturation swing) instead of doubling blindly, bounded by a new `maxntsubsteps_cap` (namelist `imaxntsubsteps_cap`, `Carma.maxntsubsteps_cap`, default 16384). The `nsubsteps` predictor is re-enabled with a gas-exhaustion criterion over each group's full depletion list, active only for condensing layers and bounded by the legacy `maxsubsteps` (truly stiff layers are handled by the informed retry). Raises the reachable time resolution from 1024 to 16384 substeps; on the 5-condensate t1100 deck this removes the substep-ceiling saturation that was manufacturing condensate non-conservatively (Fe column mass was inflated 25x at the old cap).
- Added quartz (SiO2) condensate (vp: Grant et al. 2023; surften/density: Kiefer et al. 2026)
- Added full stoichiometric depletion vectors (`stoich`) to all condensates in preparation for competitive condensation
- growth.txt v2 format with per-pathway stoichiometric coefficient and key-gas flag (legacy 2-column files still load)
- Fixed SiO gas composition flag (was copy-pasted from Al2O3)
- Fixed `plot_condensation_curves` to read vapor-pressure parameters from the group (was silently broken since they moved off the gas); one curve per condensate
- Engine: condensate properties (vapor pressure, surface tension, stoichiometry, ...) now live on the group instead of the gas, so multiple condensates can share one gas reservoir. Supersaturation, vapor pressure, Kelvin factor, latent heat, and desorption energy are tracked per condensate. Verified bitwise-identical to the previous engine on the tio2_hom, tio2_coag, and mg2sio4_het pathways.
- Engine: gas depletion applies each condensate's own stoichiometric factor (previously, condensates sharing a gas would all silently use one condensate's factor); latent-heat release likewise uses per-condensate latent heats weighted by gas mass
- Engine: contact angles are keyed by (seed element, condensate element), so two condensates nucleating on the same seed keep distinct contact angles
- Engine: fixed an uninitialized air density in the diffusivity of non-H2O/H2SO4 gases when water is not the first gas processed
- Engine driver: `CARMASTATE_PrepareStep`/`CARMASTATE_GetGroup` replace the per-gas saturation plumbing in `SetGas`/`GetGas`; the growth-file read loop now uses NGROWTH (was NNUC)
- The saturation vapor pressure column in the main output now reports, for each gas, the lowest saturation vapor pressure among its condensates (output format unchanged)
- Binary restart files from earlier versions are incompatible when the number of groups differs from the number of gases (the restart snapshot now stores per-group saturation arrays)
- Competitive condensation: a condensate can now deplete multiple gas reservoirs (many-to-many depletion). Each condensate's `stoich` vector generates secondary depletion pathways automatically (e.g. Mg2SiO4 consumes 2 Mg + 1 SiO; MgSiO3 consumes 1 SiO + 1 Mg), disable with `carma.do_secondary_depletion = False`. Growth kinetics remain driven by the single key gas (Helling & Woitke S_r = S^(1/nu_key)); secondary gases enter only the mass budget
- `CARMA_AddGasDepletion` engine API registers secondary stoichiometric depletion; gas reservoir abundances are initialized at the deepest cloud base among all condensates depleting that gas, key or secondary
- Fixed all gases sharing one boundary-condition dictionary: setting one gas's boundary condition silently applied it to every gas
- New closed-box integration test verifies column-integrated elemental Mg and Si conservation to ~5e-5 across the three-silicate competition; new `silicate_comp` pathway regression fixture