# Idealized MICE runtime checks

For the prescribed melt/freeze experiments, set `idealized_case=1/2` in
`namelist.icepack`. The configuration is read before initialization uses that
flag. `idealized_case=0` retains the realistic-case forcing requirements.

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
