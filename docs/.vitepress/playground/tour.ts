//! Guided tour of the playground.
//!
//! Every step drives the real WASM functions: it loads a sample or runs an
//! operation, then opens the tab where the result is read. The numbers in the
//! notes are the ones the step produces, checked by `tests/e2e/`.

export interface TourStep {
  id: string
  title: string
  area: string
  note: string
  sample?: string
  operation?: string
  params?: Record<string, any>
  tab?: string
  style?: string
  scale?: number
  brillouin?: boolean
  select?: number[]
  /// Camera angles in radians.
  view?: { azimuth?: number; elevation?: number; zoom?: number }
}

export const TOUR: TourStep[] = [
  {
    id: 'read',
    title: 'Read a CIF',
    area: 'Structure',
    note: 'Silicon from COD: one asymmetric site expands to 8 atoms, 160.19 Å³ and 2.329 g/cm³.',
    sample: 'Si',
    tab: 'structure',
  },
  {
    id: 'symmetry',
    title: 'Detect the space group',
    area: 'Symmetry',
    note: 'moyo finds Fd-3m (No. 227) with 192 operations in the conventional cell; the Wyckoff table lists the 8-fold orbit.',
    sample: 'Si',
    tab: 'symmetry',
  },
  {
    id: 'zone',
    title: 'First Brillouin zone',
    area: 'Symmetry',
    note: 'The zone is built from the primitive fcc lattice: a truncated octahedron with 24 vertices and 14 faces.',
    sample: 'Si',
    tab: 'symmetry',
    brillouin: true,
  },
  {
    id: 'measure',
    title: 'Measure bonds and angles',
    area: 'Analysis',
    note: 'Three picked atoms give the tetrahedral Si–Si–Si angle of 109.47°; two give the 2.352 Å bond.',
    sample: 'Si',
    tab: 'analysis',
    select: [4, 0, 5],
  },
  {
    id: 'supercell',
    title: 'Supercell',
    area: 'Structure',
    note: 'A 2 × 2 × 2 rock-salt supercell: 64 atoms and eight times the volume.',
    sample: 'NaCl',
    operation: 'supercell',
    params: { na: 2, nb: 2, nc: 2 },
    tab: 'structure',
  },
  {
    id: 'slab',
    title: 'Cut a surface slab',
    area: 'Structure',
    note: 'Two (1 1 1) layers of silicon with 10 Å of vacuum, ready for adsorption models.',
    sample: 'Si',
    operation: 'slab',
    params: { h: 1, k: 1, l: 1, layers: 2, vacuum: 10 },
    tab: 'structure',
    view: { azimuth: 1.57, elevation: 0.08 },
  },
  {
    id: 'pxrd',
    title: 'Powder diffraction',
    area: 'Powder',
    note: 'Rutile with Cu Kα and International Tables form factors: (1 1 0) at 27.4° is the strongest line.',
    sample: 'TiO2',
    tab: 'powder',
  },
  {
    id: 'refine',
    title: 'Refine a measured pattern',
    area: 'Powder',
    note: 'Below the calculated pattern: load the reference rutile pattern and run Le Bail, Pawley, Rietveld (it starts from a Pawley cell) or the automatic run.',
    sample: 'TiO2',
    tab: 'powder',
  },
  {
    id: 'nanotube',
    title: 'Carbon nanotube',
    area: 'Builders',
    note: 'An armchair (5, 5) tube: 6.78 Å diameter, every carbon three-coordinated.',
    operation: 'nanotube',
    params: { n: 5, m: 5, length: 3 },
    tab: 'structure',
    // The tube axis is y: look from +x so the tube lies across the view.
    view: { azimuth: 0, elevation: 0.2 },
  },
  {
    id: 'twist',
    title: 'Twisted bilayer graphene',
    area: 'Builders',
    note: 'The commensurate (2, 1) moiré cell: 21.79° between the layers and 14 atoms per layer.',
    operation: 'twist',
    params: { n: 2, m: 2, theta: 21.79, layers: 2 },
    tab: 'structure',
    view: { azimuth: 0, elevation: 1.45 },
  },
  {
    id: 'c60',
    title: 'Fullerene C60',
    area: 'Builders',
    note: 'Built from the face spiral: 12 pentagons, 20 hexagons, shown in space-fill.',
    operation: 'c60',
    tab: 'structure',
    style: 'space-fill',
  },
  {
    id: 'mof',
    title: 'Assemble a MOF',
    area: 'MOF',
    note: 'A topology, a node and a linker from the embedded catalogue give a periodic framework.',
    operation: 'mof',
    tab: 'structure',
  },
  {
    id: 'carbon',
    title: 'Disordered porous carbon',
    area: 'Builders',
    note: 'A carbon sheet wrapped around the dia net, relaxed to 1.42 Å bonds with H, F and O edges.',
    operation: 'carbon',
    params: { net: 'dia', repeat: 1, density: 0.42, seed: 0 },
    tab: 'structure',
  },
  {
    id: 'polymer',
    title: 'Polymer chain',
    area: 'Builders',
    note: 'Polystyrene backbone from the SMILES repeat unit, periodic along a.',
    operation: 'polymer',
    params: { smiles: 'C(c1ccccc1)C', monomers: 4, density: 1.0 },
    tab: 'structure',
  },
  {
    id: 'porosity',
    title: 'Porosity of a zeolite',
    area: 'Analysis',
    note: 'MFI: run the distance grid in the Analysis tab for the void fraction and pore size.',
    sample: 'MFI',
    tab: 'analysis',
  },
  {
    id: 'uff',
    title: 'Force field',
    area: 'Force field',
    note: 'UFF energy, atom types, QEq charges, Γ-point vibrations and a minimization of rutile.',
    sample: 'TiO2',
    tab: 'forcefield',
  },
  {
    id: 'export',
    title: 'Export',
    area: 'Export',
    note: 'Download CIF, POSCAR, XYZ, a printable STL, the scene JSON or a PNG of the view.',
    tab: 'export',
  },
]
