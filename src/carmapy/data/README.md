# Bundled data

## `adiabat_grad.txt`

Adiabatic gradient `∇_ad = d ln T / d ln P |_S` and specific heat `log Cp`
for a hydrogen/helium mixture with helium mass fraction `Y = 0.28`, on a
53 × 26 grid of `log10 T [K]` × `log10 P [bar]`. The gradient includes
H₂ ↔ 2H dissociation and a detailed accounting of the molecular vibrational
and rotational levels, which is what makes it differ sharply from an analytic
fit above ~2000 K.

Computed by Didier Saumon (`dsaumon@lanl.gov`), and obtained here from the
PICASO reference data (`reference/climate_INPUTS/specific_heat_p_adiabat_grad.json`).
The original JSON was converted to plain text so that both `carmapy.adiabat`
and the Fortran engine can read a single file; the numbers round-trip exactly.
The upstream description and attribution strings are reproduced verbatim in the
file header.

Read by `carmapy.adiabat` when `Carma.adiabat` is set to `"table"`, and by
the Fortran engine (`rce_load_adiabat` in `carma_rce.F90`), which is handed
this file's path through the namelist rather than a copy of it. Only the
`adiabat_grad` block is used at present; `specific_heat` is carried so the file
does not have to be regenerated if it is ever needed.

Coverage is `T ∈ [10, 3981] K` and `P ∈ [10⁻², 10³] bar`. Requests outside that
range are clamped to the table edge, with a warning.

If you use this table, cite Saumon, Chabrier & van Horn (1995), ApJS, 99, 713.
