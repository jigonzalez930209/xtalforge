# Reference: preparation

The full preparation pipeline lives in `xtalforge-ff` (`prepare`) and the CLI
`ff prepare`. The agent verbs cover the interactive subset.

## Pipeline

1. **Parse.** Read the CIF and keep the deposited file untouched.
2. **Symmetry.** Detect the space group at the requested tolerance; report the
   Hall number and the standardized cell.
3. **Occupancy.** Choose the wanted disorder component; never delete atoms
   silently.
4. **Hydrogens.** Add the missing H at element-specific X-H distances
   (C-H 1.09, N-H 1.01, O-H 0.96 A by default).
5. **Guests.** Remove solvent or guests only when the user asks; record it.
6. **Charges.** QEq, EQeq or EEQ when the model needs them.
7. **Report.** Every step, its tolerance and its edits.

## Agent verbs

| Step | Verb |
|---|---|
| symmetry | `symmetry` |
| occupancy | `info` then a deliberate edit |
| hydrogens | `hydrogens` |
| wrap | `wrap` |
| charges | `uff_energy` does not assign charges; use the CLI `ff charges` |
| review | `changes`, `history`, `undo` |

## Caveats

- A larger `symprec` can merge atoms that are genuinely distinct. Report the
  value and do not raise it until the structure becomes "right".
- Hydrogen placement assumes the perceived coordination; wrong bond rules give
  wrong hydrogens.
- The preparation pipeline does not decide which disorder component is
  chemically correct. That is a user decision.
