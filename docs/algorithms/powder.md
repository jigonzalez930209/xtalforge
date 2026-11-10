# Powder pipeline: instrument, profiles and the automatic run

The individual fits have their own pages: [indexing](./indexing),
[Le Bail](./lebail), [Pawley](./pawley), [Rietveld](./rietveld) and
[extinction screening](./extinction). This page covers the parts that connect
them.

## Instrument and profile

`Instrument` holds the geometry (`BraggBrentano` or `DebyeScherrer`), the zero
shift, the specimen displacement `s` (mm), the goniometer radius `R` (mm, default
240), the `K alpha2 / K alpha1` intensity ratio (default 0.5) and a polarisation
fraction. Two aberrations move a line from its Bragg position:

```
2 theta_reported = 2 theta_nominal + zero_shift
delta(2 theta)   = -2 s cos(theta) / R          (Bragg-Brentano only)
```

A positive displacement moves the line to a lower angle; the capillary geometry
has no displacement term. A doublet (`doublet`, `doublet_checked`) places the
`K alpha2` line by Bragg's law from the `K alpha1` position at the intensity
ratio. The polarisation field is stored and validated but the PXRD profile
applies its own Lorentz-polarisation factor and does not read it yet.

`ProfileKind` selects Gaussian, Lorentzian or pseudo-Voigt
(`eta L + (1 - eta) G`, `eta` clamped to `[0, 1]`). The calculated PXRD pattern
uses one Gaussian of the requested FWHM for every reflection; a Caglioti
width (`U, V, W`) is **not** implemented. Le Bail and Pawley default to
pseudo-Voigt. The model is explicit: no axial divergence, specimen
transparency or receiving-slit correction, and no size/strain broadening beyond
the width parameters of the fits.

## Cells and extinction

`free_cell_parameters` returns the parameters a crystal system lets vary;
`free_cell_parameters_with_setting` adds the setting (unique axis, origin
choice). `extinction_symbol` and `space_group_candidates` enumerate the groups
compatible with the observed absences. The ranking described on the extinction
page uses a Le Bail BIC and a Pawley fit per leading class; the 1 % Rwp tie rule
orders the table.

## Automatic run

![Automatic powder run](/diagrams/powder-pipeline.svg)


`run_automatic` is the orchestrator:

1. find peaks and estimate the background;
2. index the peak list and keep the leading cells;
3. for each leading cell, enumerate the extinction classes;
4. Le Bail every (cell, class), score by BIC;
5. Pawley-fit every leading pair and rank by Rwp with the tie rule;
6. optionally refine the best phase with Rietveld.

The result carries the leading-cell table, the per-class evidence, the accepted
phase and the diagnostics that explain a rejection. It is deterministic: no
random restarts, and every ranking has a documented tie-break.

## Rietx bridge

`rietx_inventory` reports whether the optional Python `rietx` package is
available and what it provides. The bridge is an inventory and acceptance
check, not a hidden subprocess: a caller decides to run it. The license note is
kept with the bridge.

## Acceptance

`check_acceptance` compares a refinement with the reference tolerances:
maximum Rwp and maximum position error. A run that exceeds them is reported as
rejected with the metric that failed. The thresholds are constants in the
crate, not a command-line default.

## Validation

- Golden indexing, Le Bail, Pawley and Rietveld data are in
  `tests/golden/phase-7/`; the rutile pattern comes from Crystal-Builder and
  rietx.
- The automatic table golden checks the class ranking and the leading cell.
- The extinction golden checks the partition and the leading rutile class.
- The CoRE MOF subset checks that a real framework yields a bounded pattern;
  it is not a refinement reference.
