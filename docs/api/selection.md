# Selection and tracking

`Selection { atoms, bonds }` is a set of atom indices and normalised bond
pairs. `add_atom`, `remove_atom`, `toggle_atom`, `add_atoms`, `add_bond`,
`contains`, `union`, `intersect` and `clear` are the editing API.

`TrackingTable` assigns a stable `AtomId` to a site and keeps id-to-site and
site-to-id lookups. `remap` rewrites the current indices after an edit and
drops ids whose site disappeared, so a selection or measurement can follow an
atom across commands.
