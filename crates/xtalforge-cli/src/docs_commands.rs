//! Shell completions and the manual page (Task 8.5.4).

use clap::CommandFactory;
use clap_complete::Shell;

use crate::{Cli, CliError};

pub(crate) fn completions(shell: Shell) -> Result<(), CliError> {
    let mut command = Cli::command();
    let name = command.get_name().to_owned();
    clap_complete::generate(shell, &mut command, name, &mut std::io::stdout());
    Ok(())
}

pub(crate) fn man() -> Result<(), CliError> {
    let command = Cli::command();
    let man = clap_mangen::Man::new(command);
    man.render(&mut std::io::stdout())
        .map_err(|error| CliError::new("io", error.to_string()))
}
