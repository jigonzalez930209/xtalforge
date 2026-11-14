//! Sample structures shipped with the docs playground.
//!
//! The files live in `docs/public/samples/` and are byte-identical copies of
//! the curated corpus in `tests/corpus/`; `tests/e2e/reference/generate.py`
//! checks that. Coordinates are CC0 (COD) or project models derived from CC0
//! data, see `docs/public/samples/README.md`.

export interface Sample {
  id: string
  label: string
  file: string
  note: string
}

export const SAMPLES: Sample[] = [
  { id: 'Si', label: 'Silicon', file: 'Si.cif', note: 'Diamond structure, Fd-3m' },
  { id: 'NaCl', label: 'Sodium chloride', file: 'NaCl.cif', note: 'Rock salt, Fm-3m' },
  { id: 'diamond', label: 'Diamond', file: 'diamond.cif', note: 'Carbon, Fd-3m' },
  { id: 'graphene', label: 'Graphene', file: 'graphene.cif', note: 'Hexagonal sheet, P6/mmm' },
  { id: 'TiO2', label: 'Rutile TiO2', file: 'TiO2.cif', note: 'Tetragonal, P4_2/mnm' },
  { id: 'MFI', label: 'MFI zeolite', file: 'MFI.cif', note: 'Silicalite-1 framework' },
  { id: 'MOF-5', label: 'MOF-5', file: 'MOF-5.cif', note: 'Zn4O nodes and BDC linkers' },
  { id: 'Fe', label: 'Iron', file: 'Fe.cif', note: 'Body-centred cubic, Im-3m' },
]

export async function loadSampleText(sample: Sample): Promise<string> {
  const base = import.meta.env.BASE_URL
  const response = await fetch(`${base}samples/${sample.file}`)
  if (!response.ok) {
    throw new Error(`could not load ${sample.file}: HTTP ${response.status}`)
  }
  return response.text()
}

export const TABS = [
  ['tour', 'Tour'],
  ['structure', 'Structure'],
  ['symmetry', 'Symmetry'],
  ['analysis', 'Analysis'],
  ['powder', 'Powder'],
  ['forcefield', 'Force field'],
  ['mof', 'MOF'],
  ['export', 'Export'],
] as const

