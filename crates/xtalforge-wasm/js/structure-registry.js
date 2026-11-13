// Registers WASM Structure handles with FinalizationRegistry so the Rust
// allocation is released when JavaScript forgets the wrapper.
const registry =
  typeof FinalizationRegistry === "undefined"
    ? null
    : new FinalizationRegistry((handle) => {
        try {
          handle.free();
        } catch {
          // The handle may already have been freed explicitly.
        }
      });

/**
 * Tracks a Structure handle and returns it unchanged.
 * @template T
 * @param {T} handle
 * @returns {T}
 */
export function trackStructure(handle) {
  if (registry && handle && typeof handle.free === "function") {
    registry.register(handle, handle);
  }
  return handle;
}
