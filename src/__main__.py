"""
Entrypoint for LOGAN CLI
"""

import click
import os
import pathlib

import ccbr_tools.pkg_util
import ccbr_tools.pipeline.util
import ccbr_tools.pipeline.nextflow


def repo_base(*paths):
    basedir = pathlib.Path(__file__).absolute().parent.parent
    return basedir.joinpath(*paths)


def print_citation_flag(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    ccbr_tools.pkg_util.print_citation(
        citation_file=repo_base("CITATION.cff"), output_format="bibtex"
    )
    ctx.exit()


@click.group(
    cls=ccbr_tools.pkg_util.CustomClickGroup,
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.version_option(
    ccbr_tools.pkg_util.get_version(repo_base=repo_base),
    "-v",
    "--version",
    is_flag=True,
)
@click.option(
    "--citation",
    is_flag=True,
    callback=print_citation_flag,
    expose_value=False,
    is_eager=True,
    help="Print the citation in bibtex format and exit.",
)
def cli():
    """whoLe genOme-sequencinG Analysis pipeliNe

    For more options, run:
    logan [command] --help"""
    pass


help_msg_extra = """
\b
Nextflow options:
-profile <profile>    Nextflow profile to use (e.g. test)
-params-file <file>   Nextflow params file to use (e.g. assets/params.yml)
-preview              Preview the processes that will run without executing them

\b
EXAMPLES:
Execute with slurm:
  logan run --output path/to/outdir --mode slurm
Preview the processes that will run:
  logan run --output path/to/outdir --mode local -preview
Add nextflow args (anything supported by `nextflow run`):
  logan run --output path/to/outdir --mode slurm -profile test
  logan run --output path/to/outdir --mode slurm -profile test -params-file assets/params.yml
"""


# DEVELOPER NOTE: cannot use single-hyphen options e.g. -m, -o or else it may clash with nextflow's cli options
# e.g. -profile clashed with -o (--output) and caused the command to be parsed as "-pr -o file"
@click.command(
    epilog=help_msg_extra,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option(
    "--main",
    "main_path",
    help="Path to the logan main.nf file or the GitHub repo (CCBR/LOGAN). Defaults to the version installed in the $PATH.",
    type=str,
    default=repo_base("main.nf"),
    show_default=True,
    hidden=True,
)
@click.option(
    "--output",
    help="Output directory path for logan init & run. Equivalient to nextflow launchDir. Defaults to your current working directory.",
    type=click.Path(file_okay=False, dir_okay=True, writable=True),
    default=pathlib.Path.cwd(),
    show_default=False,
)
@click.option(
    "--mode",
    "_mode",
    help="Run mode (slurm, local)",
    type=str,
    default="slurm",
    show_default=True,
)
@click.option(
    "--forceall",
    "-F",
    "force_all",
    help="Force all processes to run (i.e. do not use nextflow -resume)",
    is_flag=True,
    default=False,
    show_default=True,
)
@click.argument("nextflow_args", nargs=-1)
def run(main_path, output, _mode, force_all, **kwargs):
    """
    Run the workflow
    
    Note: you must first run `logan init --output <output_dir>` to initialize
    the output directory.
    """
    if (  # this is the only acceptable github repo option for logan
        main_path != "CCBR/LOGAN"
    ):
        # make sure the path exists
        if not os.path.exists(main_path):
            raise FileNotFoundError(
                f"Path to the logan main.nf file not found: {main_path}"
            )
    output_dir = output if isinstance(output, pathlib.Path) else pathlib.Path(output)
    ccbr_tools.pkg_util.msg_box("Output Directory", errmsg=str(output_dir))
    if not output_dir.is_dir() or not (output_dir / "nextflow.config").exists():
        raise FileNotFoundError(
            f"output directory not initialized: {output_dir}. Hint: you must initialize the output directory with `logan init --output {output_dir}`"
        )
    current_wd = os.getcwd()
    try:
        os.chdir(output_dir)
        ccbr_tools.pipeline.nextflow.run(
            nextfile_path=main_path,
            mode=_mode,
            force_all=force_all,
            pipeline_name="LOGAN",
            **kwargs,
        )
    finally:
        os.chdir(current_wd)


@click.command()
@click.option(
    "--output",
    help="Output directory path for logan init & run. Equivalient to nextflow launchDir. Defaults to your current working directory.",
    type=click.Path(file_okay=False, dir_okay=True, writable=True),
    default=pathlib.Path.cwd(),
    show_default=False,
)
def init(output):
    """Initialize the launch directory by copying the system default config files"""
    output_dir = output if isinstance(output, pathlib.Path) else pathlib.Path(output)
    ccbr_tools.pkg_util.msg_box(f"Initializing LOGAN in {output_dir}")
    (output_dir / "log/").mkdir(parents=True, exist_ok=True)
    paths = ("nextflow.config", "conf/", "assets/")
    ccbr_tools.pipeline.util.copy_config(paths, repo_base=repo_base, outdir=output_dir)


cli.add_command(run)
cli.add_command(init)


def main():
    cli()


cli(prog_name="logan")

if __name__ == "__main__":
    main()
