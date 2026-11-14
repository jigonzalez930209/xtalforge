//! Deterministic SLURM and PBS submission scripts (Task 18.4.2).

use std::path::{Path, PathBuf};

use crate::Error;

use super::util::{validate_relative_path, validate_token, write_atomic};
use super::{JobOptions, Scheduler};

/// Builds a deterministic submission script.
pub(crate) fn scheduler_script(options: &JobOptions) -> Result<String, Error> {
    if options.scheduler == Scheduler::Unknown {
        return Err(Error::Job(
            "unknown scheduler; expected slurm or pbs".to_owned(),
        ));
    }
    if options.command.trim().is_empty() {
        return Err(Error::Job("a job needs a command".to_owned()));
    }
    if options.nodes == 0 {
        return Err(Error::Job("a job needs at least one node".to_owned()));
    }
    validate_token(&options.job_name, "job name")?;
    if options.job_name.contains('/') || options.job_name.contains('\\') {
        return Err(Error::Job(format!(
            "job name {:?} must not contain path separators",
            options.job_name
        )));
    }
    validate_token(&options.queue, "queue")?;
    validate_token(&options.walltime, "walltime")?;
    validate_token(&options.output, "output path")?;
    validate_token(&options.error, "error path")?;
    validate_relative_path(&options.output, "output path")?;
    validate_relative_path(&options.error, "error path")?;
    if !valid_walltime(&options.walltime) {
        return Err(Error::Job(format!(
            "walltime {:?} must look like HH:MM:SS or D-HH:MM:SS",
            options.walltime
        )));
    }
    Ok(match options.scheduler {
        Scheduler::Slurm => slurm_script(options),
        Scheduler::Pbs => pbs_script(options),
        Scheduler::Unknown => {
            return Err(Error::Job(
                "unknown scheduler; expected slurm or pbs".to_owned(),
            ));
        }
    })
}

/// Writes the script into `work_dir`.
pub(crate) fn write_scheduler_script(
    options: &JobOptions,
    work_dir: impl AsRef<Path>,
) -> Result<PathBuf, Error> {
    let script = scheduler_script(options)?;
    let extension = options
        .scheduler
        .extension()
        .ok_or_else(|| Error::Job("unknown scheduler; expected slurm or pbs".to_owned()))?;
    let directory = work_dir.as_ref();
    let path = directory.join(format!("{}.{extension}", options.job_name));
    write_atomic(&path, &script)?;
    Ok(path)
}

fn valid_walltime(value: &str) -> bool {
    !value.is_empty()
        && value
            .chars()
            .all(|character| character.is_ascii_digit() || character == ':' || character == '-')
        && value.chars().any(|character| character.is_ascii_digit())
}

fn slurm_script(options: &JobOptions) -> String {
    let mut script = String::from("#!/bin/bash\n");
    script.push_str(&format!("#SBATCH --job-name={}\n", options.job_name));
    script.push_str(&format!("#SBATCH --partition={}\n", options.queue));
    script.push_str(&format!("#SBATCH --nodes={}\n", options.nodes));
    script.push_str(&format!("#SBATCH --time={}\n", options.walltime));
    script.push_str(&format!("#SBATCH --output={}\n", options.output));
    script.push_str(&format!("#SBATCH --error={}\n", options.error));
    script.push_str("set -euo pipefail\n");
    script.push_str("cd \"${SLURM_SUBMIT_DIR:-$PWD}\"\n");
    script.push_str(&options.command);
    script.push('\n');
    script
}

fn pbs_script(options: &JobOptions) -> String {
    let mut script = String::from("#!/bin/bash\n");
    script.push_str(&format!("#PBS -N {}\n", options.job_name));
    script.push_str(&format!("#PBS -q {}\n", options.queue));
    script.push_str(&format!("#PBS -l nodes={}\n", options.nodes));
    script.push_str(&format!("#PBS -l walltime={}\n", options.walltime));
    script.push_str(&format!("#PBS -o {}\n", options.output));
    script.push_str(&format!("#PBS -e {}\n", options.error));
    script.push_str("set -euo pipefail\n");
    script.push_str("cd \"${PBS_O_WORKDIR:-$PWD}\"\n");
    script.push_str(&options.command);
    script.push('\n');
    script
}
