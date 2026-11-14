//! Known engine failure patterns and their corrective actions (Task 18.4.3).

use super::FailurePattern;

/// One detection rule.
struct Rule {
    /// Stable code returned to every surface.
    code: &'static str,
    /// Lower-case fragments searched in the combined engine text.
    patterns: &'static [&'static str],
    /// Corrective action shown to the user.
    action: &'static str,
}

/// Rules in priority order.
///
/// A scheduler kill is reported before an SCF warning because a walltime or
/// memory kill often leaves a half-written convergence line behind. A missing
/// pseudopotential is reported before a numerical failure because the run
/// cannot succeed without it.
const RULES: &[Rule] = &[
    Rule {
        code: "walltime_exceeded",
        patterns: &[
            "due to time limit",
            "time limit exceeded",
            "walltime exceeded",
            "wall time exceeded",
            "job killed: walltime",
        ],
        action: "Resubmit with a longer --time/walltime and restart from the saved state \
                 (VASP: CONTCAR, CHGCAR and WAVECAR; QE: restart_mode='restart'; CP2K: RESTART). \
                 Split the calculation into shorter stages when the queue limit cannot grow.",
    },
    Rule {
        code: "out_of_memory",
        patterns: &[
            "out of memory",
            "out-of-memory",
            "oom-kill",
            "oom_kill",
            "oom killed",
            "cannot allocate memory",
            "std::bad_alloc",
            "bad_alloc",
            "memory limit exceeded",
        ],
        action: "Reduce the number of MPI ranks or threads per node, or request more memory per \
                 node. Lower the plane-wave cutoff or k-point grid, and avoid holding full \
                 WAVECAR/CHGCAR arrays in memory at once.",
    },
    Rule {
        code: "missing_pseudopotential",
        patterns: &[
            "missing pseudopotential",
            "pseudopotential not found",
            "pseudopotentials not found",
            "missing in the potcar",
            "missing potcar",
            "potcar not found",
            "potcar: no such file",
            "cannot open potcar",
            "could not find pseudopotential",
            "unable to find pseudopotential",
            "upf file not found",
            "cannot read pseudopotential",
        ],
        action: "Provide the missing pseudopotential files (VASP POTCAR, QE UPF/SSSP, CP2K GTH) \
                 and point the input at their directory. Check the species names and the POTCAR \
                 ordering, and never commit licensed POTCAR files to the repository.",
    },
    Rule {
        code: "zbrent",
        patterns: &["zbrent"],
        action: "VASP could not bracket a line-search minimum. Restart from CONTCAR with a \
                 smaller POTIM, use a gentler optimizer (IBRION=1 or 2), and check for atoms that \
                 are too close together or an inconsistent MAGMOM/LDAU setting.",
    },
    Rule {
        code: "scf_not_converged",
        patterns: &[
            "scf not converged",
            "scf convergence not achieved",
            "scf did not converge",
            "scf failed to converge",
            "electronic convergence not achieved",
            "electronic self-consistency was not achieved",
            "self-consistency not achieved",
            "convergence not achieved after",
        ],
        action: "Improve the SCF mixing (VASP: AMIX/BMIX; QE: mixing_beta; CP2K: MIXING), loosen \
                 the threshold slightly, and restart from the saved charge density (CHGCAR or \
                 WAVECAR, or QE restart_mode='restart'). If it still fails, switch to a more \
                 robust diagonalisation or smearing scheme.",
    },
];

/// The documented failure catalogue, in detection-priority order.
#[must_use]
pub fn failure_catalogue() -> Vec<FailurePattern> {
    RULES.iter().map(pattern_from_rule).collect()
}

/// The documented corrective action for a stable failure code.
#[must_use]
pub fn corrective_action(code: &str) -> Option<&'static str> {
    RULES
        .iter()
        .find(|rule| rule.code == code)
        .map(|rule| rule.action)
}

/// Detects the highest-priority known failure pattern in engine output.
///
/// The search is case-insensitive. It returns the first rule in the
/// documented priority order, so a scheduler kill wins over an SCF warning
/// left behind by the same run.
#[must_use]
pub fn detect_failure(text: &str) -> Option<FailurePattern> {
    let lower = text.to_ascii_lowercase();
    for rule in RULES {
        if let Some(pattern) = rule
            .patterns
            .iter()
            .find(|pattern| lower.contains(**pattern))
        {
            return Some(FailurePattern {
                code: rule.code.to_owned(),
                pattern: (*pattern).to_owned(),
                action: rule.action.to_owned(),
            });
        }
    }
    None
}

/// The canonical pattern for a code, used when the cause is known
/// independently of the text (a process timeout, for example).
pub(crate) fn pattern_for(code: &str) -> Option<FailurePattern> {
    RULES
        .iter()
        .find(|rule| rule.code == code)
        .map(pattern_from_rule)
}

fn pattern_from_rule(rule: &Rule) -> FailurePattern {
    FailurePattern {
        code: rule.code.to_owned(),
        pattern: rule
            .patterns
            .first()
            .copied()
            .unwrap_or(rule.code)
            .to_owned(),
        action: rule.action.to_owned(),
    }
}
