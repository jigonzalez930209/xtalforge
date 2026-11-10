# Reference: carbon builders

The carbon verbs are not agent verbs; use the CLI, Rust API, Python or WASM
builders, then `load` the result as text or a file.

| Builder | Name | Main options |
|---|---|---|
| Graphene sheet | `build graphene` | `nx`, `ny`, `layers`, `interlayer`, `vacuum` |
| Nanotube | `build nanotube` | `n`, `m`, `length` |
| Fullerene | `build c60`, `build c20` | fixed structures |
| Disordered carbon | `build carbon` | `net`, `repeat`, `density`, `radius_ratio`, `coverage`, `seed`, terminations |
| Gyroid / Schwarzite | `build ztc`, `build schwarzite` | `samples`, `cell`, `isovalue` |
| Nanoribbon | `build ribbon` | `width`, `length`, `vacuum` |
| Polymer | `build polymer` | `smiles`, `monomers`, `density` |

The disordered-carbon recipe uses a seeded pipeline: convolution surface over
a periodic net, split/collapse/flip remesh, tangential relaxation, projection
and H/F/O termination. The seed is explicit; the same seed and options give
the same structure. The result reports density, ring census, mean C-C and
periodicity.

Before analysis, run `info` and `diagnostics` on the built structure. A
termination that did not consume every edge shows up as dangling atoms or
unusual coordination.
