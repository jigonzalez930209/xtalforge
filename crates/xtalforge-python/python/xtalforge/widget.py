"""anywidget viewer for structures and trajectories (Task 21.5.1).

The widget reuses the scene arrays from :func:`xtalforge.scene_data`, so a
notebook shows the same atom colours, radii and bonds as the desktop viewer.
It draws them on a canvas with a drag-to-rotate interaction; no Three.js or
network access is needed.

Install ``anywidget`` to use it::

    from xtalforge.widget import StructureWidget
    StructureWidget(open("framework.cif").read())

The widget is deliberately small: it is a preview for notebooks, not the full
WebGL playground.
"""

from __future__ import annotations

from typing import Any

try:
    import anywidget
    import traitlets
except ImportError:  # pragma: no cover - optional dependency
    anywidget = None
    traitlets = None

_ESM = """
function render({ model, el }) {
  const canvas = document.createElement("canvas");
  canvas.width = 480;
  canvas.height = 360;
  canvas.style.width = "100%";
  canvas.style.maxWidth = "480px";
  canvas.style.cursor = "grab";
  el.appendChild(canvas);
  const label = document.createElement("div");
  label.style.font = "12px sans-serif";
  el.appendChild(label);
  const ctx = canvas.getContext("2d");
  let angle = 0.6;
  let tilt = 0.3;
  function project(point) {
    const [x, y, z] = point;
    const ca = Math.cos(angle), sa = Math.sin(angle);
    const ct = Math.cos(tilt), st = Math.sin(tilt);
    const x1 = ca * x - sa * z;
    const z1 = sa * x + ca * z;
    const y1 = ct * y - st * z1;
    return [x1, y1];
  }
  function draw() {
    const scene = model.get("scene");
    const atoms = scene.atoms || [];
    const bonds = scene.bonds || [];
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    ctx.fillStyle = "#ffffff";
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    const scale = 70;
    const cx = canvas.width / 2;
    const cy = canvas.height / 2;
    for (const bond of bonds) {
      const a = atoms[bond.i], b = atoms[bond.j];
      if (!a || !b) continue;
      const p = project(a.position), q = project(b.position);
      ctx.strokeStyle = "#777";
      ctx.lineWidth = 2;
      ctx.beginPath();
      ctx.moveTo(cx + scale * p[0], cy - scale * p[1]);
      ctx.lineTo(cx + scale * q[0], cy - scale * q[1]);
      ctx.stroke();
    }
    for (const atom of atoms) {
      const p = project(atom.position);
      ctx.fillStyle = `rgb(${atom.color[0]},${atom.color[1]},${atom.color[2]})`;
      ctx.beginPath();
      ctx.arc(cx + scale * p[0], cy - scale * p[1], 5 * atom.radius, 0, 2 * Math.PI);
      ctx.fill();
    }
    label.textContent = `${atoms.length} atoms, ${bonds.length} bonds`;
  }
  canvas.addEventListener("mousedown", (event) => {
    const start = [event.clientX, event.clientY, angle, tilt];
    function move(moveEvent) {
      angle = start[2] + (moveEvent.clientX - start[0]) * 0.01;
      tilt = start[3] + (moveEvent.clientY - start[1]) * 0.01;
      draw();
    }
    function up() {
      window.removeEventListener("mousemove", move);
      window.removeEventListener("mouseup", up);
    }
    window.addEventListener("mousemove", move);
    window.addEventListener("mouseup", up);
  });
  model.on("change:scene", draw);
  draw();
}
export default { render };
"""


if anywidget is not None:

    class StructureWidget(anywidget.AnyWidget):
        """Interactive preview of one structure or a trajectory frame."""

        _esm = _ESM
        scene = traitlets.Dict({}).tag(sync=True)
        frame = traitlets.Int(0).tag(sync=True)
        frames = traitlets.Int(1).tag(sync=True)

        def __init__(self, structure: str, trajectory: list[str] | None = None) -> None:
            super().__init__()
            self._texts = [structure] + list(trajectory or [])
            self.frames = len(self._texts)
            self._update()

        def _update(self) -> None:
            from . import scene_data

            text = self._texts[self.frame]
            scene = scene_data(text, "ball-stick", 1.0)
            self.scene = {
                "atoms": [
                    {
                        "position": atom["position"],
                        "radius": atom["radius"],
                        "color": atom["color"],
                    }
                    for atom in scene["atoms"]
                ],
                "bonds": [{"i": bond["i"], "j": bond["j"]} for bond in scene["bonds"]],
            }

        def show_frame(self, index: int) -> None:
            """Shows one trajectory frame."""
            if not 0 <= index < len(self._texts):
                raise IndexError("frame index out of range")
            self.frame = index
            self._update()

else:

    class StructureWidget:  # type: ignore[no-redef]
        """Placeholder that explains the missing optional dependency."""

        def __init__(self, *_args: Any, **_kwargs: Any) -> None:
            raise ImportError(
                "StructureWidget needs the optional `anywidget` package; "
                "install it with `pip install anywidget`"
            )
