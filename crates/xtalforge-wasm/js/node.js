// Node entry with filesystem helpers (Task 12.3.4).
import { readFile, writeFile } from "node:fs/promises";
import init, { Structure } from "../pkg/xtalforge_wasm.js";

let ready = init();

/** Reads a structure from a path. */
export async function read(path) {
  await ready;
  const text = await readFile(path, "utf8");
  return Structure.fromText(text, path.split(".").pop() || "");
}

/** Writes a structure in CIF, POSCAR or XYZ. */
export async function write(path, structure, format = "cif") {
  await ready;
  const text =
    format === "poscar"
      ? structure.toPoscar()
      : format === "cif"
        ? structure.toCif()
        : structure.toJson();
  await writeFile(path, text, "utf8");
  return text.length;
}
