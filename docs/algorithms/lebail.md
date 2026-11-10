# Le Bail intensity extraction

## Iteration

![Le Bail iteration](/diagrams/lebail-iteration.svg)


Le Bail extraction partitions the observed counts among the allowed reflections
without a structural model. Starting from the PXRD reflection intensities
`I_k^0`, each pass rebuilds the calculated pattern

```text
y_calc(x) = sum_k I_k P(x; 2theta_k, w)
```

and updates every intensity with the observed-to-calculated ratio inside its
own window:

```text
I_k <- I_k * sum_window y_obs / sum_window y_calc
```

The window is `window * fwhm` degrees around the reflection. A window with no
calculated intensity keeps the previous value. After the requested passes the
final profile sum becomes `LebailResult::pattern` and the last update the
reported intensities.

## R-factors

- `Rp = sum |y_obs - y_calc| / sum |y_obs|`
- `Rwp = sqrt(sum w (y_obs - y_calc)^2 / sum w y_obs^2)`
- `chi^2 = sum w (y_obs - y_calc)^2 / (N_obs - N_reflections)`

with `w = 1 / max(y_obs, 1)`.

## Intensity esds

At the final intensities the profile design matrix `A` is assembled over the
observed points, `A[i][k] = P(x_i; 2theta_k, w)`, and the normal matrix is
`A^T W A`. Its inverse scaled by the reduced chi-squared gives the covariance;
the diagonal is `intensity_esds`. Reflections that overlap share covariance
terms, so a weak line next to a strong one receives a larger esd than the
diagonal approximation would suggest. Above 128 reflections the exact inverse
is skipped and only `variance / (A^T W A)_kk` is used, which is documented as
the diagonal approximation.

The esds describe the precision of the partition at the **fixed** cell,
positions and profile. They do not include the uncertainty of those
parameters; a Pawley fit with the cell and width free is the tool for that
(see [Pawley fit](./pawley)).

## Validation

`tests/lebail.rs` generates a noisy silicon pattern from the PXRD reflections
with a deterministic Poisson-like noise, runs Le Bail with the matching
Gaussian profile, and requires at least 90% of the recovered intensities to be
within four esds of the generating values. The other tests check the R-factors,
the esd length and positivity, and the input validation.

## Limits

- The profile is evaluated per reflection with one FWHM; size/strain broadening
  belongs to the Pawley and Rietveld stages.
- A reflection whose window contains no calculated intensity cannot be updated
  in that pass; the iteration keeps its previous value.
- The result is a partition, not a least-squares solution: the esds are the
  covariance of that partition, and overlapping reflections are correlated.
