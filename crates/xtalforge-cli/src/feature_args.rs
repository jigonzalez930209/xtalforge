//! Clap definitions for the Phase 9 commands.

use clap::Subcommand;

use crate::format::StructureFormat;

/// Actions under `xtalforge features`.
#[derive(Debug, Subcommand)]
pub enum FeatureAction {
    /// Random displacements with a fixed seed.
    Rattle {
        /// Input structure.
        path: String,
        /// Gaussian standard deviation in Angstrom.
        #[arg(long, default_value_t = 0.05)]
        amplitude: f64,
        /// Deterministic seed.
        #[arg(long, default_value_t = 0)]
        seed: u64,
        /// Draw bounded instead of Gaussian displacements.
        #[arg(long)]
        bounded: bool,
        /// Output path.
        #[arg(long)]
        out: Option<String>,
        /// Output format.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Add, remove, replace, move, sort or wrap sites.
    Edit {
        /// Input structure.
        path: String,
        /// Add `element:x,y,z` (fractional).
        #[arg(long = "add")]
        add: Vec<String>,
        /// Remove a site index.
        #[arg(long = "remove")]
        remove: Vec<usize>,
        /// Replace `index:element`.
        #[arg(long = "replace")]
        replace: Vec<String>,
        /// Move `index:x,y,z`.
        #[arg(long = "move")]
        move_site: Vec<String>,
        /// Replace species `old:new`.
        #[arg(long = "species")]
        species: Vec<String>,
        /// Fraction for partial species replacement.
        #[arg(long, default_value_t = 1.0)]
        fraction: f64,
        /// Seed for partial species replacement.
        #[arg(long, default_value_t = 0)]
        seed: u64,
        /// Sort sites by `element`, `electronegativity`, `mass` or `label`.
        #[arg(long)]
        sort: Option<String>,
        /// Wrap coordinates into [0, 1).
        #[arg(long)]
        wrap: bool,
        /// Output path.
        #[arg(long)]
        out: Option<String>,
        /// Output format.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Linear interpolation between two matched structures.
    Interpolate {
        /// First structure.
        first: String,
        /// Second structure.
        second: String,
        /// Number of intermediate frames.
        #[arg(long, default_value_t = 5)]
        images: u32,
        /// Output path; a `.xyz` path writes all frames.
        #[arg(long)]
        out: Option<String>,
        /// Output format.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Supercell closest to a cube for a target cell count.
    OptimalSupercell {
        /// Input structure.
        path: String,
        /// Number of original cells.
        cells: u32,
        /// Output path.
        #[arg(long)]
        out: Option<String>,
        /// Output format.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Periodic Voronoi cells.
    Voronoi {
        /// Input structure.
        path: String,
        /// Neighbor cutoff in Angstrom.
        #[arg(long, default_value_t = 6.0)]
        cutoff: f64,
    },
    /// Coordination shell and shape measures of one atom.
    Coordination {
        /// Input structure.
        path: String,
        /// P1 atom index.
        index: usize,
        /// Strategy: cutoff, voronoi or crystalnn.
        #[arg(long, default_value = "voronoi")]
        strategy: String,
        /// Cutoff in Angstrom.
        #[arg(long, default_value_t = 5.0)]
        cutoff: f64,
        /// Minimum Voronoi solid angle.
        #[arg(long, default_value_t = 0.1)]
        solid_angle: f64,
    },
    /// Bond-valence sums.
    Bvs {
        /// Input structure.
        path: String,
        /// Neighbor cutoff in Angstrom.
        #[arg(long, default_value_t = 5.0)]
        cutoff: f64,
    },
    /// Oxidation states from bond valence.
    Oxidation {
        /// Input structure.
        path: String,
        /// Neighbor cutoff in Angstrom.
        #[arg(long, default_value_t = 5.0)]
        cutoff: f64,
    },
    /// Total and partial radial distribution functions.
    Rdf {
        /// Input structure.
        path: String,
        /// Maximum distance in Angstrom.
        #[arg(long, default_value_t = 10.0)]
        rmax: f64,
        /// Bin width in Angstrom.
        #[arg(long, default_value_t = 0.05)]
        bin: f64,
        /// Smoothing window.
        #[arg(long, default_value_t = 0)]
        smoothing: usize,
    },
    /// Structure fingerprint and similarity.
    Fingerprint {
        /// First structure.
        first: String,
        /// Second structure.
        second: String,
        /// Maximum RDF distance.
        #[arg(long, default_value_t = 10.0)]
        rmax: f64,
        /// RDF bin width.
        #[arg(long, default_value_t = 0.1)]
        bin: f64,
    },
    /// Symmetry-unique Miller indices.
    Miller {
        /// Input structure.
        path: String,
        /// Maximum absolute index.
        #[arg(long, default_value_t = 2)]
        max_index: i32,
        /// Symmetry tolerance in Angstrom.
        #[arg(long, default_value_t = 1e-5)]
        symprec: f64,
    },
    /// All terminations of one surface.
    Terminations {
        /// Input structure.
        path: String,
        /// Miller index `h,k,l`.
        #[arg(long)]
        hkl: String,
        /// Number of layers.
        #[arg(long, default_value_t = 2)]
        layers: u32,
        /// Vacuum in Angstrom.
        #[arg(long, default_value_t = 10.0)]
        vacuum: f64,
        /// Mirror the slab into a symmetric one.
        #[arg(long)]
        symmetric: bool,
        /// Termination to write; defaults to the first.
        #[arg(long)]
        index: Option<usize>,
        /// Output path.
        #[arg(long)]
        out: Option<String>,
        /// Output format.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Wulff construction from surface energies.
    Wulff {
        /// Input structure.
        path: String,
        /// Surface energy `h,k,l,gamma`; repeat for more directions.
        #[arg(long = "surface")]
        surfaces: Vec<String>,
    },
    /// Adsorption sites on a slab and optional molecule placement.
    Adsorption {
        /// Input slab.
        path: String,
        /// Tolerance for the top layer in fractional units.
        #[arg(long, default_value_t = 0.05)]
        tolerance: f64,
        /// Optional molecule to place on the first site.
        #[arg(long)]
        molecule: Option<String>,
        /// Offset along the surface normal in Angstrom.
        #[arg(long, default_value_t = 2.0)]
        offset: f64,
        /// Output path.
        #[arg(long)]
        out: Option<String>,
        /// Output format.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Symmetry-unique vacancies.
    Vacancies {
        /// Input structure.
        path: String,
        /// Minimum distance between defect images in Angstrom.
        #[arg(long, default_value_t = 10.0)]
        min_image: f64,
        /// Defect to write; defaults to the first.
        #[arg(long)]
        index: Option<usize>,
        /// Output path.
        #[arg(long)]
        out: Option<String>,
        /// Output format.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Symmetry-unique substitutions.
    Substitutions {
        /// Input structure.
        path: String,
        /// Replacing species.
        species: String,
        /// Minimum distance between defect images in Angstrom.
        #[arg(long, default_value_t = 10.0)]
        min_image: f64,
        /// Defect to write; defaults to the first.
        #[arg(long)]
        index: Option<usize>,
        /// Output path.
        #[arg(long)]
        out: Option<String>,
        /// Output format.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Interstitial candidates from Voronoi vertices.
    Interstitials {
        /// Input structure.
        path: String,
        /// Inserted species.
        species: String,
        /// Minimum distance from existing atoms in Angstrom.
        #[arg(long, default_value_t = 1.0)]
        clearance: f64,
        /// Voronoi cutoff in Angstrom.
        #[arg(long, default_value_t = 5.0)]
        cutoff: f64,
        /// Minimum distance between defect images in Angstrom.
        #[arg(long, default_value_t = 10.0)]
        min_image: f64,
    },
    /// Enumerate and rank orderings of partial occupancies.
    Orderings {
        /// Input structure.
        path: String,
        /// Supercell `a,b,c`.
        #[arg(long, default_value = "1,1,1")]
        supercell: String,
        /// Maximum number of decorations.
        #[arg(long, default_value_t = 64)]
        max_orderings: usize,
        /// Ordering to write; defaults to the lowest energy.
        #[arg(long)]
        index: Option<usize>,
        /// Output path.
        #[arg(long)]
        out: Option<String>,
        /// Output format.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Special quasirandom structure for a binary disorder.
    Sqs {
        /// Input structure with partial occupancy.
        path: String,
        /// Supercell `a,b,c`.
        #[arg(long, default_value = "2,2,2")]
        supercell: String,
        /// Deterministic seed.
        #[arg(long, default_value_t = 0)]
        seed: u64,
        /// Pair cutoff in Angstrom.
        #[arg(long, default_value_t = 6.0)]
        cutoff: f64,
        /// Number of swap attempts.
        #[arg(long, default_value_t = 5000)]
        iterations: usize,
        /// Output path.
        #[arg(long)]
        out: Option<String>,
        /// Output format.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Commensurate 2D lattice matches between a film and a substrate.
    LatticeMatch {
        /// Film structure.
        film: String,
        /// Substrate structure.
        substrate: String,
        /// Maximum matrix entry.
        #[arg(long, default_value_t = 3)]
        max_index: i32,
        /// Maximum matched area in square Angstrom.
        #[arg(long, default_value_t = 200.0)]
        max_area: f64,
        /// Maximum relative length mismatch.
        #[arg(long, default_value_t = 0.05)]
        tolerance: f64,
    },
    /// Stack a film slab on a substrate slab.
    Stack {
        /// Film slab.
        film: String,
        /// Substrate slab.
        substrate: String,
        /// Gap in Angstrom.
        #[arg(long, default_value_t = 2.0)]
        gap: f64,
        /// Vacuum above the film in Angstrom.
        #[arg(long, default_value_t = 15.0)]
        vacuum: f64,
        /// In-plane shift `u,v` in substrate fractions.
        #[arg(long, default_value = "0,0")]
        shift: String,
        /// Output path.
        #[arg(long)]
        out: Option<String>,
        /// Output format.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Elastic tensor from stress finite differences.
    Elastic {
        /// Input structure.
        path: String,
        /// Applied engineering strain.
        #[arg(long, default_value_t = 0.005)]
        strain: f64,
    },
    /// Equation-of-state fits from energy-volume points.
    Eos {
        /// `volume,energy` point; repeat for more points.
        #[arg(long = "point")]
        points: Vec<String>,
        /// Model: birch, vinet or murnaghan.
        #[arg(long, default_value = "birch")]
        model: String,
    },
    /// Formation energies and convex-hull distances from a JSON phase list.
    Hull {
        /// JSON file with `[{ "name", "composition", "energy" }]`.
        path: String,
    },
}
