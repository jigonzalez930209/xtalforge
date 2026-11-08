//! Command dispatch for the CLI.

use crate::{
    Cli, CliError, analysis_commands, analyze, args, args::Command, basic, build, build_structures,
    commands_extra, convert, ff_commands, measure, mesh, mof_commands, powder_commands, pxrd,
    spacegroup, symmetry, symmetry::SymmetryAction, symmetry_analysis,
};

pub(crate) fn run(cli: &Cli) -> Result<(), CliError> {
    if let Some(threads) = cli.threads
        && threads == 0
    {
        return Err(CliError::input("--threads must be a positive number"));
    }
    match &cli.command {
        Command::Hydrogens {
            path,
            distance,
            out,
            format,
        } => build_structures::hydrogens(cli, path, *distance, out.as_deref(), *format),
        Command::Powder { action } => powder_commands::run(cli, action),
        Command::Mof { action } => mof_commands::run(cli, action),
        Command::Analysis { action } => analysis_commands::run(cli, action),
        Command::Ff { action } => ff_commands::run(cli, action),
        Command::Formats => commands_extra::formats(cli),
        Command::Params { name } => commands_extra::params(cli, name.as_deref()),
        Command::Project { action } => match action {
            args::ProjectAction::Pack { path, out } => commands_extra::project_pack(cli, path, out),
            args::ProjectAction::Info { path } => commands_extra::project_info(cli, path),
            args::ProjectAction::Extract { path, dir } => {
                commands_extra::project_extract(cli, path, dir)
            }
        },
        Command::Trajectory {
            path,
            frame,
            out,
            format,
        } => commands_extra::trajectory(cli, path, *frame, out.as_deref(), *format),
        Command::Script { path } => commands_extra::script(cli, path),
        Command::Version => {
            let version = xtalforge::core::version();
            if cli.json {
                println!("{}", serde_json::json!({ "version": version }));
            } else if cli.quiet {
                println!("{version}");
            } else {
                println!("xtalforge {version}");
            }
            Ok(())
        }
        Command::Element { input } => {
            let element = xtalforge::core::element_info(input).map_err(CliError::core)?;
            if cli.json {
                println!(
                    "{}",
                    serde_json::json!({
                        "z": element.z,
                        "symbol": element.symbol,
                        "name": element.name,
                        "mass": element.mass,
                        "covalent_radius": element.covalent_radius,
                        "vdw_radius": element.vdw_radius,
                        "color": element.color,
                        "is_metal": element.is_metal,
                        "valence": element.valence,
                        "electronegativity": element.electronegativity,
                        "pyykko_single": element.pyykko_single,
                        "pyykko_double": element.pyykko_double,
                        "pyykko_aromatic": element.pyykko_aromatic,
                    })
                );
            } else if cli.quiet {
                println!("{}", element.symbol);
            } else {
                let [red, green, blue] = element.color;
                println!(
                    "{} {} (Z={}) mass={:.4} covalent={:.3} vdw={:.3} color=#{red:02X}{green:02X}{blue:02X}",
                    element.symbol,
                    element.name,
                    element.z,
                    element.mass,
                    element.covalent_radius,
                    element.vdw_radius
                );
            }
            Ok(())
        }
        Command::Info { path } => basic::info(cli, path),
        Command::Cell { params } => basic::cell(cli, params),
        Command::CellVolume {
            a,
            b,
            c,
            alpha,
            beta,
            gamma,
        } => basic::cell_volume(cli, *a, *b, *c, *alpha, *beta, *gamma),
        Command::Convert { path, format, out } => convert::run(cli, path, *format, out.as_deref()),
        Command::Pxrd {
            path,
            radiation,
            max_two_theta,
            step,
            out,
        } => pxrd::run(cli, path, radiation, *max_two_theta, *step, out.as_deref()),
        Command::Build { action } => match action {
            build_structures::BuildAction::Graphene {
                nx,
                ny,
                layers,
                interlayer,
                vacuum,
                out,
                format,
            } => build_structures::graphene(
                cli,
                build_structures::GrapheneOptions {
                    nx: *nx,
                    ny: *ny,
                    layers: *layers,
                    interlayer: *interlayer,
                    vacuum: *vacuum,
                    out: out.as_deref(),
                    format: *format,
                },
            ),
            build_structures::BuildAction::Nanotube {
                n,
                m,
                length,
                vacuum,
                out,
                format,
            } => build_structures::nanotube(cli, *n, *m, *length, *vacuum, out.as_deref(), *format),
            build_structures::BuildAction::C20 { out, format } => {
                build_structures::c20(cli, out.as_deref(), *format)
            }
            build_structures::BuildAction::C60 { out, format } => {
                build_structures::c60(cli, out.as_deref(), *format)
            }
            build_structures::BuildAction::Twist {
                n,
                m,
                theta,
                layers,
                out,
                format,
            } => build_structures::twist(cli, *n, *m, *theta, *layers, out.as_deref(), *format),
            build_structures::BuildAction::Ribbon {
                width,
                length,
                vacuum,
                out,
                format,
            } => build_structures::ribbon(cli, *width, *length, *vacuum, out.as_deref(), *format),
            build_structures::BuildAction::Ztc {
                samples,
                cell,
                isovalue,
                out,
                format,
            } => build_structures::ztc(cli, *samples, *cell, *isovalue, out.as_deref(), *format),
            build_structures::BuildAction::Schwarzite {
                samples,
                cell,
                isovalue,
                out,
                format,
            } => build_structures::schwarzite(
                cli,
                *samples,
                *cell,
                *isovalue,
                out.as_deref(),
                *format,
            ),
            build_structures::BuildAction::Polymer {
                smiles,
                monomers,
                density,
                out,
                format,
            } => {
                build_structures::polymer(cli, smiles, *monomers, *density, out.as_deref(), *format)
            }
        },
        Command::Mesh {
            path,
            style,
            scale,
            out,
        } => mesh::run(cli, path, style, *scale, out.as_deref()),
        Command::Supercell {
            path,
            na,
            nb,
            nc,
            out,
            format,
        } => build::supercell(cli, path, [*na, *nb, *nc], out.as_deref(), *format),
        Command::Slab {
            path,
            h,
            k,
            l,
            layers,
            vacuum,
            out,
            format,
        } => build::slab(
            cli,
            path,
            [*h, *k, *l],
            *layers,
            *vacuum,
            out.as_deref(),
            *format,
        ),
        Command::Bonds { path } => analyze::bonds(cli, path),
        Command::Rings { path, max_size } => analyze::rings(cli, path, *max_size),
        Command::Groups { path } => analyze::groups(cli, path),
        Command::Neighbors { path, cutoff } => analyze::neighbors(cli, path, *cutoff),
        Command::Distance { path, i, j } => measure::distance(cli, path, *i, *j),
        Command::Angle { path, i, j, k } => measure::angle(cli, path, *i, *j, *k),
        Command::Spacegroup { action } => match action {
            spacegroup::SpacegroupAction::Show { symbol } => spacegroup::show(cli, symbol),
            spacegroup::SpacegroupAction::List {
                system,
                centering,
                point_group,
            } => spacegroup::list(cli, system.as_deref(), *centering, point_group.as_deref()),
            spacegroup::SpacegroupAction::Search { text } => spacegroup::search(cli, text),
        },
        Command::Symmetry { action } => match action {
            SymmetryAction::Detect {
                path,
                symprec,
                angle_tolerance,
                hall_number,
            } => symmetry::detect(
                cli,
                path,
                *symprec,
                angle_tolerance.map(f64::to_radians),
                *hall_number,
            ),
            SymmetryAction::P1 { path, out, format } => {
                symmetry::p1(cli, path, out.as_deref(), *format)
            }
            SymmetryAction::Asymmetrize {
                path,
                symprec,
                out,
                format,
            } => symmetry::asymmetrize(cli, path, *symprec, out.as_deref(), *format),
            SymmetryAction::Set {
                path,
                symbol,
                mode,
                out,
                format,
            } => symmetry::set(cli, path, symbol, mode, out.as_deref(), *format),
            SymmetryAction::Standardize {
                path,
                symprec,
                out,
                format,
            } => symmetry_analysis::standardize(cli, path, *symprec, out.as_deref(), *format),
            SymmetryAction::Primitive {
                path,
                symprec,
                out,
                format,
            } => symmetry_analysis::primitive(cli, path, *symprec, out.as_deref(), *format),
            SymmetryAction::Reduce {
                path,
                algorithm,
                out,
                format,
            } => symmetry_analysis::reduce(cli, path, algorithm, out.as_deref(), *format),
            SymmetryAction::Wyckoff { path, symprec } => {
                symmetry_analysis::wyckoff(cli, path, *symprec)
            }
            SymmetryAction::Layer { path, symprec } => {
                symmetry_analysis::layer(cli, path, *symprec)
            }
            SymmetryAction::Magnetic { path, symprec } => {
                symmetry_analysis::magnetic(cli, path, *symprec)
            }
            SymmetryAction::Pointgroup { path, tolerance } => {
                symmetry_analysis::pointgroup(cli, path, *tolerance)
            }
            SymmetryAction::Kpath { path } => symmetry_analysis::kpath(cli, path),
            SymmetryAction::Setting {
                path,
                hall_number,
                out,
                format,
            } => symmetry_analysis::setting(cli, path, *hall_number, out.as_deref(), *format),
            SymmetryAction::Brillouin { path } => symmetry_analysis::brillouin(cli, path),
            SymmetryAction::Scan { path } => symmetry_analysis::scan(cli, path),
            SymmetryAction::Diagnose { path, symbol } => {
                symmetry_analysis::diagnose(cli, path, symbol)
            }
            SymmetryAction::Subgroups {
                hall_number,
                epsilon,
            } => symmetry::subgroups(cli, *hall_number, *epsilon),
            SymmetryAction::Lower {
                path,
                subgroup,
                tolerance,
                out,
                format,
            } => symmetry::lower(cli, path, *subgroup, *tolerance, out.as_deref(), *format),
        },
        Command::Torsion { path, i, j, k, l } => measure::torsion(cli, path, *i, *j, *k, *l),
    }
}

pub(crate) fn report_error(cli: &Cli, error: &CliError) {
    if cli.json {
        let payload = serde_json::json!({
            "error": {
                "code": error.code,
                "message": error.message,
            }
        });
        eprintln!("{payload}");
    } else {
        eprintln!("error: {}", error.message);
    }
}
