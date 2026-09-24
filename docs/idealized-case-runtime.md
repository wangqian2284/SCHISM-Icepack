# Idealized MICE runtime checks

For the prescribed melt/freeze experiments, set `idealized_case=1/2` in
`namelist.icepack`. The configuration is read before initialization uses that
flag. `idealized_case=0` retains the realistic-case forcing requirements.

Cases 1 and 2 initialize both Icepack `rhow` and MICE `rhowat` to 1000 kg/m3.
Icepack derived constants (including `cprho`) and MICE `inv_rhowat`/`cc` are
updated with these densities before thermodynamic initialization. Realistic
cases retain the original densities (1026 and 1025 kg/m3 respectively).
This does not change initial ice concentration/thickness, ocean salinity,
SCHISM `rho0`/`shw`, or Icepack `cp_ocn`; configure the ocean inputs for a
freshwater experiment separately. Rebuild to apply this source change.

The `idealized_nml` group also configures sub-grid SST mixing for both realistic
and idealized cases. `sstn_kappa_e` is the dynamic mixing diffusivity in m2/s
(default 184; finite and nonnegative). Dynamic beta uses
`ice_dt*(relative_speed/mean_radius + sstn_kappa_e/mean_radius**2)` and is
clipped to [0,1] before mixing. Fixed mixing still uses `sstn_beta_fixed`
directly, with no timescale parameter. For a 100 s ice step and a 6 hour
relaxation timescale, set it to `0.00462962963`; recalculate it as
`ice_dt/21600` when the ice step changes. The sample includes this value;
the legacy default of 1.0 when fixed beta is omitted is unchanged.
Existing files that omit `sstn_kappa_e` now use 184 instead of the old
hard-coded 54. Set `sstn_kappa_e=54` explicitly to reproduce that diffusivity.

Idealized initial ice uses zero salinity for both the enthalpy calculation and
the ice salinity tracers, matching the original freshwater melt experiment.
With `ktherm=2`, ice at 0 C then has enthalpy `-rhoi*Lfresh`, retaining its
latent heat requirement. Using the default nonzero ice salinity at 0 C instead
initializes fully liquid mush even when the assigned thickness is 1 m.
Realistic initialization keeps its original salinity profile. Ocean salinity
is not changed by this initialization fix. No new input parameter is required;
rebuild and start from `ihot_mice=0` to regenerate the initial state. Existing
restart states are not repaired by this change.

The standalone prescribed-forcing setup can use `nws=0`, `ihconsv=0`, and
`isconsv=0` in `param.nml`. In idealized mode, the ocean still applies the heat,
shortwave and freshwater-derived salinity fluxes returned by Icepack. These
fluxes continue to apply to open-water nodes after ice disappears, allowing
subsequent cooling/freezing. This corrects a previously inactive ocean exchange
path; results can therefore differ from runs that simply removed the abort.

When `ihconsv=0`, idealized mode initializes ocean albedo to 0.06 and water type
to 1 (Jerlov I). These are explicit experiment defaults, not inferred from the
grid. When `ihconsv=1`, the normal `albedo.gr3` and `watertype.gr3` inputs are
still read. Prescribed cases use zero incident wave spectrum and atmospheric
pressure 101325 Pa, plus the existing ice-loading contribution on the hydro side.

Keep `ice_tests=0` and `ice_therm_on=1` in `mice.nml`: the legacy `ice_tests=1`
branch deliberately skips thermodynamics. For a fresh experiment use
`ihot_mice=0`; restart modes still require their restart files. For a stationary
column experiment, `ievp=0` and `ice_advection=0` disable ice dynamics/transport.

Other legitimate checks remain: `nstep_ice=1`, `tr_fsd=.true.`, `nfsdcat>=16`,
valid initial concentration/thickness, and a supported triangular ice grid.
If no geographic coordinates were loaded, idealized MICE uses `slam0/sfea0`;
other enabled forcing/modules may independently require `hgrid.ll`. In
particular, selecting `nws=2` still selects the external atmospheric input path.

Validation performed: compiled and executed the prescribed-atmosphere source
block with signaling-NaN inputs and floating-point traps for melt and freeze;
compiled and executed heat/salt gates with and without `USE_MICE`; checked that
the realistic-mode guard remains active, configuration precedes use, and optical
defaults follow allocation. Full MPI model execution has not been verified.

## Thin-ice sensitivity experiment

Set `idealized_case=1`, `idealized_hice=0.2` and `sstntest=2` for the
thin-ice dynamic-scheme sensitivity experiment. Keep all other baseline melt
inputs unchanged (including concentration, floe size, diffusivity, ocean state
and time step). Use `idealized_hice=1.0` for the original baseline. No separate
case number is needed; case 3 is not supported.
Start with `ihot_mice=0` so a restart does not overwrite the initial ice.
The initial thickness category is selected from 0.2 m automatically.
Enable `'meltt'` in `nml_list_icepack` to output surface melt by category.

## Freezing versus sub-grid melting

Set `subgrid_freezing` in `idealized_nml` (applies to realistic and idealized
cases, independently of `sstntest`). The default `.false.` selects original
`add_new_ice` using grid-mean `frzmlt`; FSD welding uses the same grid-mean
potential. Positive freezing heat is returned through the existing `hocnn`
coupling. `.true.` selects experimental `add_new_ice2` with local `frzmltn`
and its own heat feedback; welding then uses open-water `frzmltn(1)`.

With `.false.`, if grid SST or any occupied sub-region SST is below freezing,
all sub-region SSTs are reset to grid SST before computing melt/freeze
potentials. Empty ice categories and zero-area open water do not trigger this
fallback. Thus warm-grid/cold-open-water and cold-grid/warm-open-water states
both return to the grid-mean treatment; so do cold occupied ice categories.
The ocean grid SST itself is not changed. This is a temporary complete-mixing
fallback, not a completed sub-grid freezing or latent-heat redistribution model.
It preserves sub-region mean heat only when the pre-reset area-weighted SST
matches grid SST; pre-existing discrepancies are reset to the ocean state.
`beta=-1` marks this fallback rather than an actual mixing coefficient.

With `.true.`, the experimental behavior is retained: grid SST below freezing
still resets all sub-region SSTs, but warm-grid/local-supercooling states can
form local new ice while other categories melt. Keep this option off for the
current experiments. Surface melting and conductive basal growth remain
governed by their own energy balances; these are not disabled by the fallback.

Fallback logging: `mirror.out` contains one rank-0 line for each ice step with
at least one reset, e.g. `SST_FALLBACK step=100 wet_nodes=136 subgrid_freezing=F`.
`wet_nodes` is the global number of owned wet nodes reset during that call;
ghost nodes and initialization (`it_main<=0`) are excluded. These are per-step
counts, not unique nodes or accumulated counts over an output interval. Both
the grid-supercooling reset and the additional local-supercooling fallback are
included. No line is emitted on steps with zero counted resets.
