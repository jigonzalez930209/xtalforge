// Typed object API for the xtalforge WASM package (Task 12.3.2).
// The legacy text/dict exports are documented as deprecated and keep their
// generated loose signatures until 2.0; this interface file is fully typed.

export interface StructureSummary {
  name: string;
  formula: string;
  hall: string;
  n_sites: number;
  n_atoms: number;
  volume: number;
  density: number;
  total_charge: number;
  is_p1: boolean;
  warnings: StructureWarning[];
}

export interface StructureWarning {
  code: string;
  message: string;
  atoms: number[];
  fields: string[];
  hint?: string;
}

export interface SceneBuffers {
  positions: Float64Array;
  radii: Float32Array;
  colors: Uint8Array;
  bond_positions: Float64Array;
  bond_colors: Uint8Array;
  cell_edges: Float64Array;
}

export interface PxrdBuffers {
  two_theta: Float64Array;
  intensity: Float64Array;
}

export interface MeshBuffers {
  vertices: Float64Array;
  normals: Float64Array;
  indices: Uint32Array;
}

export interface StructureHandle {
  readonly formula: string;
  readonly nSites: number;
  readonly nAtoms: number;
  readonly positions: Float64Array;
  readonly fracPositions: Float64Array;
  readonly cell: Float64Array;
  readonly numbers: Uint8Array;
  toJson(): string;
  toCif(): string;
  toPoscar(): string;
  toSummary(): StructureSummary;
  supercell(na: number, nb: number, nc: number): StructureHandle;
  slab(h: number, k: number, l: number, layers: number, vacuum: number): StructureHandle;
  p1(): StructureHandle;
  detectSymmetry(symprec: number): string;
  pxrd(): string;
  ffEnergy(): number;
  addHydrogens(distance: number): StructureHandle;
  setPositions(positions: Float64Array): StructureHandle;
  free(): void;
}

export interface StructureConstructor {
  fromText(text: string, formatHint: string): StructureHandle;
  fromCif(text: string): StructureHandle;
}

export function sceneBuffers(text: string, formatHint: string): SceneBuffers;
export function pxrdBuffers(text: string, formatHint: string): PxrdBuffers;
export function meshBuffers(text: string, formatHint: string): MeshBuffers;

export function fetch_optimade(url: string): Promise<string>;
export function fetch_cod(id: string): Promise<string>;

export function mesh_obj(text: string, style: string, scale: number): string;
export function mesh_gltf(text: string, style: string, scale: number): string;
