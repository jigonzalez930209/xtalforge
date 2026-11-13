// Optional Web Worker wrapper for the heavy operations (Task 12.3.4).
// Usage: const worker = new Worker(new URL("xtalforge-wasm/js/worker.js", import.meta.url), { type: "module" });
//        worker.postMessage({ id: 1, op: "supercell", text, hint: "cif", args: [2, 2, 2] });
import init, { Structure } from "../pkg/xtalforge_wasm.js";

let ready = init();

self.onmessage = async (event) => {
  const { id, op, text, hint, args } = event.data;
  await ready;
  try {
    const structure = Structure.fromText(text, hint || "");
    let result;
    switch (op) {
      case "supercell":
        result = structure.supercell(...args).toJson();
        break;
      case "slab":
        result = structure.slab(...args).toJson();
        break;
      case "p1":
        result = structure.p1().toJson();
        break;
      case "pxrd":
        result = structure.pxrd();
        break;
      case "symmetry":
        result = structure.detectSymmetry(args[0]);
        break;
      default:
        throw new Error(`unknown worker operation ${op}`);
    }
    self.postMessage({ id, result });
  } catch (error) {
    self.postMessage({ id, error: String(error && error.message ? error.message : error) });
  }
};
