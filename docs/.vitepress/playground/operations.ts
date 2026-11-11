//! Structure-producing playground operations.
//!
//! Each operation calls one text-returning WASM export and returns the new
//! structure as CIF; the playground then loads it like a pasted file.

export type Params = Record<string, any>

export function runOperation(wasm: any, text: string, name: string, params: Params = {}): string {
  switch (name) {
    case 'supercell':
      return wasm.supercell_text(text, params.na ?? 2, params.nb ?? 2, params.nc ?? 2)
    case 'slab':
      return wasm.slab_text(text, params.h ?? 1, params.k ?? 0, params.l ?? 0, params.layers ?? 2, params.vacuum ?? 10)
    case 'p1':
      return wasm.to_p1_text(text)
    case 'asymmetrize':
      return wasm.asymmetrize_text(text, params.tolerance ?? 0.001)
    case 'primitive':
      return wasm.primitive_text(text, params.tolerance ?? 0.001)
    case 'standardize':
      return wasm.standardize_text(text, params.tolerance ?? 0.001)
    case 'reduce':
      return wasm.reduce_text(text, params.method ?? 'niggli')
    case 'graphene':
      return wasm.build_graphene_text(params.nx ?? 4, params.ny ?? 4, params.layers ?? 1, 3.35, 15)
    case 'nanotube':
      return wasm.build_nanotube_text(params.n ?? 5, params.m ?? 5, params.length ?? 2, 10)
    case 'c60':
      return wasm.build_c60_text()
    case 'c20':
      return wasm.build_c20_text()
    case 'ribbon':
      return wasm.build_ribbon_text(params.width ?? 6, params.length ?? 6, params.vacuum ?? 12)
    case 'twist':
      return wasm.build_twist_text(params.n ?? 1, params.m ?? 1, params.theta ?? 21.79, params.layers ?? 2)
    case 'ztc':
      return wasm.build_ztc_text(params.samples ?? 16, params.cell ?? 12, params.isovalue ?? 0)
    case 'schwarzite':
      return wasm.build_schwarzite_text(params.samples ?? 16, params.cell ?? 12, params.isovalue ?? 0)
    case 'carbon':
      return wasm.build_carbon_text(params.net ?? 'dia', params.repeat ?? 1, params.density ?? 0.42, params.seed ?? 0)
    case 'polymer':
      return wasm.build_polymer_text(params.smiles ?? 'CC', params.monomers ?? 8, params.density ?? 1.0)
    case 'mof':
      return wasm.build_mof_text(params.topology ?? 'pcu', params.node ?? 'Zn4O', params.linker ?? 'benzene')
    case 'hydrogens':
      return wasm.hydrogens_text(text, params.distance ?? 1.0)
    default:
      throw new Error(`unknown operation ${name}`)
  }
}
