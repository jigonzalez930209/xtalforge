//! Element parameters for QEq and EQeq.
//!
//! **QEq.** The `chi` and `hardness` columns are the Rappé-Goddard
//! electronegativity and idempotential in eV.  They are the UFF table's
//! `chi` and twice its `Hard` column, because UFF writes its energy as
//! `chi Q + Hard Q^2` while the QEq paper writes `chi Q + 1/2 J Q^2`, and the
//! two must not be mixed.  This is the same conversion Crystal-Builder
//! `xtal/ff/uff/qeq.py` documents as `HARDNESS_TO_J`.
//!
//! **EQeq.** The rows are `[EA, IE1, IE2, IE3, IE4, IE5]` in eV, from NIST
//! Atomic Spectra Database ionisation energies and PubChem electron
//! affinities, as tabulated by Crystal-Builder `xtal/ff/charges/data/
//! ionization.csv`.  Missing values are `NaN`.  The table here covers the
//! common elements up to bismuth with enough ionisation stages for their
//! usual charge centres; an element outside it is an explicit
//! [`crate::Error::UnknownElement`], never a guessed number.

use crate::Error;

/// One QEq parameter pair.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct QeqParams {
    /// Electronegativity in eV.
    pub chi: f64,
    /// Idempotential (hardness) in eV.
    pub hardness: f64,
}

/// QEq parameters for common elements, from the UFF table.
const QEQ_TABLE: &[(&str, f64, f64)] = &[
    ("H", 4.528, 13.8904),
    ("He", 9.660, 29.8400),
    ("Li", 3.006, 4.7720),
    ("Be", 4.877, 8.8860),
    ("B", 5.110, 9.5000),
    ("C", 5.343, 10.1260),
    ("N", 6.899, 11.7600),
    ("O", 8.741, 13.3640),
    ("F", 10.874, 14.9480),
    ("Ne", 11.040, 21.1000),
    ("Na", 2.843, 4.5920),
    ("Mg", 3.951, 7.3860),
    ("Al", 4.060, 7.1800),
    ("Si", 4.168, 6.9740),
    ("P", 5.463, 8.0000),
    ("S", 6.928, 8.9720),
    ("Cl", 8.564, 9.8920),
    ("Ar", 9.465, 12.7100),
    ("K", 2.421, 3.8400),
    ("Ca", 3.231, 5.7600),
    ("Sc", 3.395, 6.1600),
    ("Ti", 3.470, 6.7600),
    ("V", 3.650, 6.8200),
    ("Cr", 3.415, 7.7300),
    ("Mn", 3.325, 8.2100),
    ("Fe", 3.760, 8.2800),
    ("Co", 4.105, 8.3500),
    ("Ni", 4.465, 8.4100),
    ("Cu", 4.200, 8.4400),
    ("Zn", 5.106, 8.5700),
    ("Ga", 3.641, 6.3200),
    ("Ge", 4.051, 6.8760),
    ("As", 5.188, 7.6180),
    ("Se", 6.428, 8.2620),
    ("Br", 7.790, 8.8500),
    ("Kr", 8.505, 11.4300),
    ("Rb", 2.331, 3.6920),
    ("Sr", 3.024, 4.8800),
    ("Y", 3.830, 5.6200),
    ("Zr", 3.400, 7.1000),
    ("Nb", 3.550, 6.7600),
    ("Mo", 3.465, 7.5100),
    ("Tc", 3.290, 7.9800),
    ("Ru", 3.575, 8.0300),
    ("Rh", 3.975, 8.0100),
    ("Pd", 4.320, 8.0000),
    ("Ag", 4.436, 6.2680),
    ("Cd", 5.034, 7.9140),
    ("In", 3.506, 5.7920),
    ("Sn", 3.987, 6.2480),
    ("Sb", 4.899, 6.6840),
    ("Te", 5.816, 7.0520),
    ("I", 6.822, 7.5240),
    ("Xe", 7.595, 9.9500),
    ("Cs", 2.183, 3.4220),
    ("Ba", 2.814, 4.7920),
    ("La", 2.836, 5.4830),
    ("Ce", 2.774, 5.3840),
    ("Pr", 2.858, 5.1280),
    ("Nd", 2.869, 5.2410),
    ("Pm", 2.881, 5.3460),
    ("Sm", 2.912, 5.4390),
    ("Eu", 2.878, 5.5750),
    ("Gd", 3.167, 5.9490),
    ("Tb", 3.018, 5.6680),
    ("Dy", 3.055, 5.7430),
    ("Ho", 3.127, 5.7820),
    ("Er", 3.187, 5.8290),
    ("Tm", 3.251, 5.8658),
    ("Yb", 3.289, 5.9300),
    ("Lu", 2.963, 4.9258),
    ("Hf", 3.700, 6.8000),
    ("Ta", 5.100, 5.7000),
    ("W", 4.630, 6.6200),
    ("Re", 3.960, 7.8400),
    ("Os", 5.140, 7.2600),
    ("Ir", 5.000, 8.0000),
    ("Pt", 4.790, 8.8600),
    ("Au", 4.894, 5.1720),
    ("Hg", 6.270, 8.3200),
    ("Tl", 3.200, 5.8000),
    ("Pb", 3.900, 7.0600),
    ("Bi", 4.690, 7.4800),
    ("Po", 4.210, 6.2300),
    ("At", 4.750, 6.0000),
    ("Rn", 5.370, 7.0000),
    ("Fr", 2.000, 4.0000),
    ("Ra", 2.843, 4.8680),
    ("Ac", 2.835, 5.6700),
    ("Th", 3.175, 5.8100),
    ("Pa", 2.985, 5.8100),
    ("U", 3.341, 5.7060),
    ("Np", 3.549, 5.4340),
    ("Pu", 3.243, 5.6380),
    ("Am", 2.990, 6.0070),
    ("Cm", 2.832, 6.3790),
    ("Bk", 3.193, 6.0710),
    ("Cf", 3.197, 6.2020),
    ("Es", 3.333, 6.1780),
    ("Fm", 3.400, 6.2000),
    ("Md", 3.470, 6.2200),
    ("No", 3.475, 6.3500),
    ("Lr", 3.500, 6.4000),
];

/// QEq parameters for one element.
pub(crate) fn qeq_params(element: &str) -> Result<QeqParams, Error> {
    let symbol = if element == "D" { "H" } else { element };
    QEQ_TABLE
        .iter()
        .find(|(name, _, _)| *name == symbol)
        .map(|(_, chi, hardness)| QeqParams {
            chi: *chi,
            hardness: *hardness,
        })
        .ok_or_else(|| Error::UnknownElement(symbol.to_owned()))
}

/// `NaN` in a compile-time table.
const NAN: f64 = f64::NAN;

/// EQeq rows: `EA, IE1, IE2, IE3, IE4, IE5`, in eV.
const EQEQ_TABLE: &[(&str, [f64; 6])] = &[
    ("H", [0.754, 13.598_434_599_702, NAN, NAN, NAN, NAN]),
    ("He", [NAN, 24.587_389_011, 54.417_765_528_2, NAN, NAN, NAN]),
    (
        "Li",
        [0.618, 5.391_714_996, 75.640_097_0, 122.454_359_13, NAN, NAN],
    ),
    (
        "Be",
        [NAN, 9.322_699, 18.211_15, 153.896_205, 217.718_584_59, NAN],
    ),
    (
        "B",
        [
            0.277,
            8.298_019,
            25.154_83,
            37.930_59,
            259.374_379,
            340.226_022_5,
        ],
    ),
    (
        "C",
        [
            1.263,
            11.260_288_0,
            24.383_143,
            47.887_78,
            64.493_52,
            392.090_56,
        ],
    ),
    (
        "N",
        [NAN, 14.534_13, 29.601_25, 47.445_3, 77.473_5, 97.890_1],
    ),
    (
        "O",
        [
            1.461, 13.618_055, 35.121_12, 54.935_54, 77.413_50, 113.899_0,
        ],
    ),
    (
        "F",
        [3.339, 17.422_82, 34.970_81, 62.707_98, 87.175, 114.249],
    ),
    (
        "Ne",
        [NAN, 21.564_541, 40.962_97, 63.423_3, 97.190_0, 126.247],
    ),
    (
        "Na",
        [0.548, 5.139_076_96, 47.286_36, 71.620_0, 98.936, 138.404],
    ),
    (
        "Mg",
        [NAN, 7.646_236, 15.035_271, 80.143_6, 109.265_4, 141.33],
    ),
    (
        "Al",
        [
            0.441, 5.985_769, 18.828_55, 28.447_642, 119.992_4, 153.825_2,
        ],
    ),
    (
        "Si",
        [1.385, 8.151_68, 16.345_85, 33.493_00, 45.141_79, 166.767],
    ),
    (
        "P",
        [
            0.746, 10.486_686, 19.769_49, 30.202_64, 51.443_87, 65.025_11,
        ],
    ),
    (
        "S",
        [2.077, 10.360_016_7, 23.337_88, 34.86, 47.222, 72.594_5],
    ),
    ("Cl", [3.617, 12.967_633, 23.813_64, 39.80, 53.24, 67.68]),
    ("Ar", [NAN, 15.759_611_9, 27.629_67, 40.735, 59.58, 74.84]),
    (
        "K",
        [0.501, 4.340_663_73, 31.625_00, 45.803_1, 60.917, 82.66],
    ),
    (
        "Ca",
        [NAN, 6.113_154_921_0, 11.871_719, 50.913_16, 67.273_2, 84.34],
    ),
    (
        "Sc",
        [0.188, 6.561_49, 12.799_77, 24.756_839, 73.489_4, 91.95],
    ),
    (
        "Ti",
        [0.079, 6.828_120, 13.575_5, 27.491_71, 43.267_17, 99.299],
    ),
    ("V", [0.525, 6.746_187, 14.634, 29.311_1, 46.709, 65.281_65]),
    ("Cr", [0.666, 6.766_51, 16.486_305, 30.959, 49.16, 69.46]),
    ("Mn", [NAN, 7.434_038_0, 15.639_99, 33.668, 51.21, 72.41]),
    ("Fe", [0.163, 7.902_468_1, 16.199_21, 30.651, 54.91, 75.00]),
    ("Co", [0.661, 7.881_01, 17.084_4, 33.50, 51.27, 79.50]),
    ("Ni", [1.156, 7.639_878, 18.168_838, 35.187, 54.92, 76.06]),
    ("Cu", [1.228, 7.726_380, 20.292_39, 36.841, 57.38, 79.8]),
    ("Zn", [NAN, 9.394_197, 17.964_39, 39.723_30, 59.573, 82.6]),
    (
        "Ga",
        [0.3, 5.999_302_0, 20.515_14, 30.725_76, 63.241, 86.01],
    ),
    (
        "Ge",
        [1.35, 7.899_435, 15.934_610, 34.057_6, 45.715_5, 90.500],
    ),
    ("As", [0.81, 9.788_55, 18.589_2, 28.349, 50.15, 62.77]),
    ("Se", [2.021, 9.752_368, 21.196, 31.697, 42.947, 68.30]),
    ("Br", [3.365, 11.813_81, 21.591, 34.871, 47.782, 59.595]),
    ("Kr", [NAN, 13.999_605_5, 24.359_84, 35.838, 50.85, 64.69]),
    ("Rb", [0.468, 4.177_128_1, 27.289_54, 39.247, 52.20, 68.44]),
    (
        "Sr",
        [NAN, 5.694_867_45, 11.030_276_5, 42.883_53, 56.280, 70.7],
    ),
    ("Y", [0.307, 6.217_26, 12.223_6, 20.524_41, 60.607_2, 75.35]),
    ("Zr", [0.426, 6.634_126, 13.13, 23.170, 34.418_36, 80.348]),
    ("Nb", [0.893, 6.758_85, 14.32, 25.04, 37.611, 50.572_8]),
    ("Mo", [0.746, 7.092_43, 16.16, 27.13, 40.33, 54.417]),
    ("Ag", [1.302, 7.576_234, 21.484_4, 34.8, 49.0, 65.0]),
    ("Cd", [NAN, 8.993_820, 16.908_313, 37.468, 51.0, 67.9]),
    ("In", [0.3, 5.786_355_8, 18.870_41, 28.044_15, 55.45, 69.3]),
    ("Sn", [1.2, 7.343_918, 14.633_07, 30.506, 40.74, 77.03]),
    ("Sb", [1.07, 8.608_389, 16.626, 25.323_5, 43.804, 55.00]),
    ("Te", [1.971, 9.009_808, 18.6, 27.84, 37.415_5, 59.3]),
    ("I", [3.059, 10.451_236, 19.131_26, 29.570, 40.357, 51.52]),
    ("Xe", [NAN, 12.129_843_7, 20.975, 31.05, 42.20, 54.1]),
    (
        "Cs",
        [0.472, 3.893_905_727_43, 23.157_45, 33.195, 43.0, 56.0],
    ),
    ("Ba", [NAN, 5.211_664_6, 10.003_826, 35.843_8, 47.0, 58.0]),
    ("La", [0.5, 5.576_9, 11.184_96, 19.177_3, 49.95, 61.6]),
    ("Ce", [0.5, 5.538_6, 10.956, 20.197_4, 36.906, 65.55]),
    ("Hf", [NAN, 6.825_070, 14.61, 22.55, 33.370, 68.37]),
    ("Ta", [0.322, 7.549_571, 16.2, 23.1, 35.0, 48.272]),
    ("W", [0.815, 7.864_03, 16.37, 26.0, 38.2, 51.6]),
    ("Re", [0.15, 7.833_52, 16.6, 27.0, 39.1, 51.9]),
    ("Os", [1.1, 8.438_23, 17.0, 25.0, 41.0, 55.0]),
    ("Ir", [1.565, 8.967_02, 17.0, 28.0, 40.0, 57.0]),
    ("Pt", [2.128, 8.958_83, 18.56, 29.0, 43.0, 56.0]),
    ("Au", [2.309, 9.225_554, 20.203, 30.0, 45.0, 60.0]),
    ("Hg", [NAN, 10.437_504, 18.756_88, 34.49, 48.55, 61.20]),
    ("Tl", [0.2, 6.108_287_3, 20.428_3, 29.852_0, 51.14, 62.6]),
    (
        "Pb",
        [0.36, 7.416_679_9, 15.032_499, 31.937_3, 42.332_56, 68.8],
    ),
    ("Bi", [0.946, 7.285_516, 16.703, 25.570_75, 45.37, 54.856]),
];

/// Hydrogen's electron affinity as EQeq uses it, in eV.
///
/// The measured 0.754 eV makes hydrogen far too ready to take an electron;
/// the published method substitutes -2.0 eV, and the reference implementation
/// documents it as a parameter of the method rather than an atomic property.
const HYDROGEN_AFFINITY: f64 = -2.0;

/// The common oxidation state each metal is expanded about in EQeq.
fn charge_centre(element: &str) -> i32 {
    match element {
        "Li" | "Na" | "K" | "Rb" | "Cs" | "Tl" => 1,
        "Be" | "Mg" | "Ca" | "Sr" | "Ba" | "Cr" | "Mn" | "Co" | "Ni" | "Cu" | "Zn" | "Cd"
        | "Hg" | "Pb" => 2,
        "Al" | "Sc" | "Ti" | "V" | "Fe" | "Ga" | "Y" | "Ru" | "Rh" | "La" | "Ce" | "Au" | "Bi" => 3,
        "Ge" | "Zr" | "Hf" | "Os" | "Ir" | "Pt" => 4,
        _ => 0,
    }
}

/// EQeq `(chi, hardness)` in eV for one element, moved to `Q = 0`.
///
/// `centre` overrides the table when the framework's oxidation state is known
/// to differ from the common one.
pub(crate) fn eqeq_params(element: &str, centre: Option<i32>) -> Result<(f64, f64), Error> {
    let symbol = if element == "D" { "H" } else { element };
    let row = EQEQ_TABLE
        .iter()
        .find(|(name, _)| *name == symbol)
        .map(|(_, values)| *values)
        .ok_or_else(|| Error::UnknownElement(symbol.to_owned()))?;
    let centre = centre.unwrap_or_else(|| charge_centre(symbol));
    if centre < 0 {
        return Err(Error::InvalidInput(format!(
            "EQeq cannot expand {symbol} about a negative charge"
        )));
    }
    let index = usize::try_from(centre).map_err(|_| {
        Error::InvalidInput(format!("EQeq charge centre for {symbol} is too large"))
    })?;
    let below_index = index;
    let above_index = index + 1;
    let mut below = *row.get(below_index).ok_or_else(|| {
        Error::InvalidInput(format!(
            "EQeq has no ionisation energy for {symbol} about charge {centre:+}"
        ))
    })?;
    let above = *row.get(above_index).ok_or_else(|| {
        Error::InvalidInput(format!(
            "EQeq has no ionisation energy for {symbol} about charge {centre:+}"
        ))
    })?;
    if symbol == "H" && centre == 0 {
        below = HYDROGEN_AFFINITY;
    }
    if !below.is_finite() || !above.is_finite() {
        return Err(Error::InvalidInput(format!(
            "EQeq has no ionisation energy for {symbol} about charge {centre:+}"
        )));
    }
    let hardness = above - below;
    Ok((
        0.5 * (above + below) - f64::from(centre) * hardness,
        hardness,
    ))
}
