//! Thin wrapper around the WASM package for the playground.
//!
//! Every number shown in the playground comes from the exported functions:
//! `read_cif`, `bonds`, `rings`, `detect_symmetry`, `wyckoff_info`,
//! `pxrd_info` and the `*_text` transforms. Nothing is recomputed in
//! JavaScript.

import type { SceneData } from './scene'

export interface RingStat {
  size: number
  count: number
}

export interface Inspection {
  summary: any
  bonds: any[]
  bondsCount: number
  rings: any[]
  ringStats: RingStat[]
  symmetry: any | null
  wyckoff: any[]
  error: string
}

export interface SceneStats {
  atoms: number
  bonds: number
  cellEdges: number
  data: SceneData
}

let pending: Promise<any> | null = null

export async function loadWasm(): Promise<any> {
  if (!pending) {
    pending = (async () => {
      const module = await import('xtalforge-wasm')
      await module.default()
      return module
    })()
  }
  return pending
}

export function inspect(wasm: any, text: string): Inspection {
  const empty: Inspection = {
    summary: null,
    bonds: [],
    bondsCount: 0,
    rings: [],
    ringStats: [],
    symmetry: null,
    wyckoff: [],
    error: '',
  }
  try {
    const summary = wasm.read_cif(text)
    let bonds: any[] = []
    let rings: any[] = []
    try {
      bonds = wasm.bonds(text)
    } catch {
      bonds = []
    }
    // Ring perception is quadratic on large frameworks; keep the initial
    // panel responsive and let the user ask for the full list.
    if ((summary?.n_atoms ?? 0) <= 600) {
      try {
        rings = wasm.rings(text, 8)
      } catch {
        rings = []
      }
    }
    const counts = new Map<number, number>()
    for (const ring of rings) counts.set(ring.size, (counts.get(ring.size) ?? 0) + 1)
    const ringStats = [...counts.entries()]
      .sort((a, b) => a[0] - b[0])
      .map(([size, count]) => ({ size, count }))
    let symmetry = null
    try {
      symmetry = wasm.detect_symmetry(text, 0.001)
    } catch {
      symmetry = null
    }
    let wyckoff: any[] = []
    try {
      wyckoff = wasm.wyckoff_info(text, 0.001)
    } catch {
      wyckoff = []
    }
    return { ...empty, summary, bonds, bondsCount: bonds.length, rings, ringStats, symmetry, wyckoff }
  } catch (caught) {
    return { ...empty, error: errorText(caught) }
  }
}

export function sceneStats(wasm: any, text: string, style: string, scale: number): SceneStats {
  const data: SceneData = wasm.scene_data(text, style, scale)
  return {
    atoms: data.atoms.length,
    bonds: data.bonds.length,
    cellEdges: data.cell_edges.length,
    data,
  }
}

export function errorText(caught: unknown): string {
  if (typeof caught === 'object' && caught !== null) {
    const value = caught as { code?: string; message?: string }
    if (value.code && value.message) return `${value.code}: ${value.message}`
    if (value.message) return value.message
  }
  return String(caught)
}
