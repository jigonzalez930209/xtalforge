//! Three.js renderer for the scene data produced by the Rust view crate.
//!
//! The playground never rebuilds geometry: it draws the atom positions,
//! radii, colours, bond endpoints and cell edges that `xtalforge-wasm`
//! returns from `scene_data`, so the 3D frame is the same one the Rust API,
//! the CLI (`mesh --scene`) and Python (`scene_data`) describe.

import * as THREE from 'three'
import { OrbitControls } from 'three/addons/controls/OrbitControls.js'

export interface SceneAtom {
  element: string
  position: [number, number, number]
  radius: number
  color: [number, number, number]
  index: number
  hidden: boolean
}

export interface SceneBond {
  i: number
  j: number
  start: [number, number, number]
  end: [number, number, number]
  radius: number
  color: [number, number, number]
  order: number
}

export interface SceneData {
  atoms: SceneAtom[]
  bonds: SceneBond[]
  cell_edges: Array<[[number, number, number], [number, number, number]]>
  settings: {
    style: string
    atom_scale: number
    bond_radius: number
    show_cell: boolean
    show_axes: boolean
    background: [number, number, number]
  }
}

export interface ProjectedAtom {
  index: number
  element: string
  x: number
  y: number
  depth: number
}

const SPHERE = new THREE.SphereGeometry(1, 24, 16)
const CYLINDER = new THREE.CylinderGeometry(1, 1, 1, 12, 1, true)
const Y_AXIS = new THREE.Vector3(0, 1, 0)

export class SceneRenderer {
  readonly canvas: HTMLCanvasElement
  private readonly renderer: THREE.WebGLRenderer
  private readonly scene = new THREE.Scene()
  private readonly camera: THREE.PerspectiveCamera
  private readonly controls: OrbitControls
  private readonly root = new THREE.Group()
  private readonly cellGroup = new THREE.Group()
  private readonly axesGroup = new THREE.Group()
  private atoms: SceneAtom[] = []
  private center = new THREE.Vector3()
  private baseDistance = 10
  private observer: ResizeObserver
  private running = true

  constructor(container: HTMLElement) {
    this.renderer = new THREE.WebGLRenderer({
      antialias: true,
      preserveDrawingBuffer: true,
      alpha: false,
    })
    this.renderer.setPixelRatio(1)
    this.renderer.outputColorSpace = THREE.SRGBColorSpace
    this.canvas = this.renderer.domElement
    this.canvas.setAttribute('data-testid', 'playground-canvas')
    container.appendChild(this.canvas)

    this.camera = new THREE.PerspectiveCamera(45, 1, 0.05, 20000)
    this.camera.up.set(0, 0, 1)
    this.controls = new OrbitControls(this.camera, this.canvas)
    this.controls.enableDamping = false
    this.controls.autoRotateSpeed = 1.4
    this.controls.addEventListener('change', () => this.renderOnce())

    this.scene.background = new THREE.Color(0xffffff)
    this.scene.add(this.root)
    this.root.add(this.cellGroup)
    this.root.add(this.axesGroup)
    this.scene.add(new THREE.HemisphereLight(0xffffff, 0x555555, 1.7))
    const key = new THREE.DirectionalLight(0xffffff, 1.6)
    key.position.set(4, -6, 8)
    this.scene.add(key)
    const fill = new THREE.DirectionalLight(0xffffff, 0.7)
    fill.position.set(-6, 4, -3)
    this.scene.add(fill)

    this.observer = new ResizeObserver(() => this.resize())
    this.observer.observe(container)
    this.resize()
    this.frame()
  }

  setScene(data: SceneData) {
    this.clear()
    this.atoms = data.atoms.filter((atom) => !atom.hidden)
    this.scene.background = new THREE.Color(
      data.settings.background[0] / 255,
      data.settings.background[1] / 255,
      data.settings.background[2] / 255,
    )
    this.addAtoms(this.atoms)
    this.addBonds(data.bonds)
    this.cellEdges = data.cell_edges
    if (data.settings.show_cell) {
      this.addCell(data.cell_edges)
    }
    if (data.settings.show_axes) {
      this.addAxes(data.cell_edges)
    }
    this.resetView()
  }

  private addAtoms(atoms: SceneAtom[]) {
    if (atoms.length === 0) return
    const material = new THREE.MeshStandardMaterial({
      roughness: 0.45,
      metalness: 0.05,
    })
    const mesh = new THREE.InstancedMesh(SPHERE, material, atoms.length)
    const matrix = new THREE.Matrix4()
    const color = new THREE.Color()
    const scale = new THREE.Vector3()
    atoms.forEach((atom, index) => {
      const radius = atom.radius > 0 ? atom.radius : 0.001
      scale.set(radius, radius, radius)
      matrix.compose(new THREE.Vector3(...atom.position), new THREE.Quaternion(), scale)
      mesh.setMatrixAt(index, matrix)
      color.setRGB(
        atom.color[0] / 255,
        atom.color[1] / 255,
        atom.color[2] / 255,
        THREE.SRGBColorSpace,
      )
      mesh.setColorAt(index, color)
    })
    mesh.instanceMatrix.needsUpdate = true
    if (mesh.instanceColor) mesh.instanceColor.needsUpdate = true
    this.root.add(mesh)
  }

  private addBonds(bonds: SceneBond[]) {
    const visible = bonds.filter((bond) => bond.radius > 0)
    if (visible.length === 0) return
    const material = new THREE.MeshStandardMaterial({
      roughness: 0.55,
      metalness: 0.05,
    })
    const mesh = new THREE.InstancedMesh(CYLINDER, material, visible.length)
    const matrix = new THREE.Matrix4()
    const color = new THREE.Color()
    const quaternion = new THREE.Quaternion()
    const midpoint = new THREE.Vector3()
    const direction = new THREE.Vector3()
    const scale = new THREE.Vector3()
    visible.forEach((bond, index) => {
      const start = new THREE.Vector3(...bond.start)
      const end = new THREE.Vector3(...bond.end)
      direction.subVectors(end, start)
      const length = direction.length()
      midpoint.addVectors(start, end).multiplyScalar(0.5)
      quaternion.setFromUnitVectors(Y_AXIS, direction.normalize())
      scale.set(bond.radius, Math.max(length, 1e-4), bond.radius)
      matrix.compose(midpoint, quaternion, scale)
      mesh.setMatrixAt(index, matrix)
      color.setRGB(
        bond.color[0] / 255,
        bond.color[1] / 255,
        bond.color[2] / 255,
        THREE.SRGBColorSpace,
      )
      mesh.setColorAt(index, color)
    })
    mesh.instanceMatrix.needsUpdate = true
    if (mesh.instanceColor) mesh.instanceColor.needsUpdate = true
    this.root.add(mesh)
  }

  private addCell(edges: SceneData['cell_edges']) {
    const points: number[] = []
    for (const [start, end] of edges) {
      points.push(...start, ...end)
    }
    const geometry = new THREE.BufferGeometry()
    geometry.setAttribute('position', new THREE.Float32BufferAttribute(points, 3))
    this.cellGroup.add(
      new THREE.LineSegments(
        geometry,
        new THREE.LineBasicMaterial({ color: 0x2f2f2f, transparent: true, opacity: 0.9 }),
      ),
    )
  }

  private addAxes(edges: SceneData['cell_edges']) {
    const colors = [0xc0392b, 0x27ae60, 0x2c6fbb]
    const origin = new THREE.Vector3(0, 0, 0)
    const vectors: THREE.Vector3[] = []
    for (const [start, end] of edges) {
      if (new THREE.Vector3(...start).distanceTo(origin) > 1e-6) continue
      const vector = new THREE.Vector3(...end)
      if (vectors.some((known) => known.distanceTo(vector) < 1e-6)) continue
      vectors.push(vector)
    }
    vectors.slice(0, 3).forEach((vector, index) => {
      const geometry = new THREE.BufferGeometry().setFromPoints([origin, vector])
      this.axesGroup.add(
        new THREE.Line(geometry, new THREE.LineBasicMaterial({ color: colors[index] })),
      )
    })
  }

  private clear() {
    for (const child of [...this.root.children]) {
      if (child === this.cellGroup || child === this.axesGroup) {
        for (const item of [...child.children]) {
          child.remove(item)
          item.traverse((object) => {
            if (object instanceof THREE.Mesh || object instanceof THREE.Line) {
              object.geometry.dispose()
              const material = object.material
              if (material instanceof THREE.Material) material.dispose()
            }
          })
        }
        continue
      }
      this.root.remove(child)
      child.traverse((object) => {
        if (object instanceof THREE.Mesh || object instanceof THREE.Line) {
          const geometry = object.geometry
          if (geometry !== SPHERE && geometry !== CYLINDER) geometry.dispose()
          const material = object.material
          if (material instanceof THREE.Material) material.dispose()
        }
      })
    }
  }

  setCellVisible(visible: boolean) {
    this.cellGroup.visible = visible
    this.renderOnce()
  }

  setAxesVisible(visible: boolean) {
    this.axesGroup.visible = visible
    this.renderOnce()
  }

  resetView() {
    const box = new THREE.Box3()
    let empty = true
    for (const atom of this.atoms) {
      box.expandByPoint(new THREE.Vector3(...atom.position))
      empty = false
    }
    for (const [start, end] of this.cellEdges) {
      box.expandByPoint(new THREE.Vector3(...start))
      box.expandByPoint(new THREE.Vector3(...end))
      empty = false
    }
    if (empty) box.expandByPoint(new THREE.Vector3(0, 0, 0))
    box.getCenter(this.center)
    const radius = Math.max(box.getSize(new THREE.Vector3()).length() / 2, 0.5)
    this.baseDistance = (radius / Math.sin((this.camera.fov * Math.PI) / 360)) * 1.15
    this.setView({ azimuth: 0.85, elevation: 0.42, zoom: 1 })
  }

  private cellEdges: SceneData['cell_edges'] = []

  setView({ azimuth = 0.85, elevation = 0.42, zoom = 1 }) {
    const distance = this.baseDistance / Math.max(zoom, 0.05)
    const cos = Math.cos(elevation)
    const direction = new THREE.Vector3(
      cos * Math.cos(azimuth),
      cos * Math.sin(azimuth),
      Math.sin(elevation),
    )
    this.camera.position.copy(this.center).addScaledVector(direction, distance)
    this.camera.near = Math.max(distance / 500, 0.01)
    this.camera.far = distance * 20
    this.camera.updateProjectionMatrix()
    this.controls.target.copy(this.center)
    this.controls.update()
    this.renderOnce()
  }

  setAutoRotate(enabled: boolean) {
    this.controls.autoRotate = enabled
  }

  resize() {
    const parent = this.canvas.parentElement
    if (!parent) return
    const width = Math.max(parent.clientWidth, 1)
    const height = Math.max(parent.clientHeight, 1)
    this.renderer.setSize(width, height, false)
    this.canvas.style.width = `${width}px`
    this.canvas.style.height = `${height}px`
    this.camera.aspect = width / height
    this.camera.updateProjectionMatrix()
    this.renderOnce()
  }

  renderOnce() {
    this.renderer.render(this.scene, this.camera)
  }

  snapshot(): string {
    this.renderOnce()
    return this.canvas.toDataURL('image/png')
  }

  projectAtoms(): ProjectedAtom[] {
    const width = this.renderer.domElement.width
    const height = this.renderer.domElement.height
    return this.atoms.map((atom) => {
      const vector = new THREE.Vector3(...atom.position).project(this.camera)
      return {
        index: atom.index,
        element: atom.element,
        x: (vector.x + 1) / 2,
        y: (1 - vector.y) / 2,
        depth: vector.z,
      }
    })
  }

  private frame() {
    if (!this.running) return
    this.controls.update()
    this.renderer.render(this.scene, this.camera)
    requestAnimationFrame(() => this.frame())
  }

  dispose() {
    this.running = false
    this.observer.disconnect()
    this.controls.dispose()
    this.clear()
    this.renderer.dispose()
  }
}
