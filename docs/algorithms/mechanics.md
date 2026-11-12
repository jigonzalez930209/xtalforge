# Mechanical and thermodynamic basics

**Tasks 9.6.1-9.6.3.** References: pymatgen `analysis.elasticity` and
`analysis.eos`; Birch (1947); Vinet, Ferrante, Smith and Rose (1987);
Murnaghan (1944); Voigt (1889), Reuss (1929) and Hill (1952).

## Elastic tensor

The stiffness is the derivative of the stress with respect to an engineering
strain. For each Voigt component `j`, two structures are built with
`+epsilon` and `-epsilon`:

```
C_ij = (sigma_i(+epsilon) - sigma_i(-epsilon)) / (2 epsilon)
```

The stress comes from the engine's analytic stress (GPa); UFF is the native
default, and any `Calculator` that reports stress can be used. The matrix is
symmetrized as `(C_ij + C_ji)/2` to remove finite-difference noise. The
compliance `S = C^-1` gives the Reuss moduli,

```
B_R = 1 / (S11 + S22 + S33 + 2 (S12 + S13 + S23))
G_R = 15 / (4 (S11 + S22 + S33) - 4 (S12 + S13 + S23)
            + 3 (S44 + S55 + S66))
```

and the Voigt moduli use the stiffness directly. The Hill values are the
arithmetic means; Young's modulus and Poisson's ratio follow from
`E = 9 B G / (3 B + G)` and `nu = (3 B - 2 G) / (2 (3 B + G))`. The stability
flag is a Cholesky factorization of the stiffness: positive definiteness is
the Born criterion in this approximation.

## Equation of state

Three models are fitted with Levenberg-Marquardt and numerical Jacobians:

- Birch-Murnaghan (third order):
  `E = E0 + 9 V0 B0 / 16 [(eta^2 - 1)^3 B0' + (eta^2 - 1)^2 (6 - 4 eta^2)]`
  with `eta = (V0/V)^(1/3)`.
- Vinet: the universal equation with `x = (V/V0)^(1/3)`.
- Murnaghan: the linear pressure derivative form.

`B0` is reported in GPa (1 eV/Angstrom^3 = 160.21766208 GPa), `B0'` is
dimensionless and the parameter errors come from `sigma^2 (J^T J)^-1` with
`sigma^2 = RSS/(n-4)`. The initial guess uses the minimum-energy point, a
three-point curvature for `B0` and `B0' = 4`.

## Formation energies and the convex hull

The formation energy per atom of a phase is

```
E_f = (E_total - sum_i n_i mu_i) / N
```

where the chemical potentials `mu_i` come from the pure-element entries or
from an explicit map. Every composition is converted to element fractions, so
all phases live in one simplex.

For each phase the hull energy is the solution of the linear program

```
minimize    sum_j w_j E_f,j
subject to  sum_j w_j x_j = x_target
            sum_j w_j = 1,  w_j >= 0
```

solved with a two-phase simplex. The decomposition is the set of weights above
`1e-9`; the distance above the hull is `E_f - hull` and a phase is stable when
it is at most `1e-6 eV/atom` above it. The formulation is independent of the
number of elements.
