import argparse
import sys
from pathlib import Path

from stdatamodels.jwst.datamodels import MultiSlitModel

from jwst.exp_to_source import exp_to_source

__all__ = []  # type: ignore[var-annotated]


class Main:
    """
    Convert exposure-based slits data to source-based data.

    Command-line interface to the exp_to_source step. For help, use::

        exp_to_source -h

    See also the exp_to_source readthedocs page: :ref:`exp_to_source`

    Parameters
    ----------
    args : str or [str, ...] or None
        The command-line arguments. This is passed to
        :py:func:`argparse.parse_args`. If None, ``sys.argv``
        is used.
    """  # fmt: skip

    def __init__(self, args=None):
        if args is None:
            args = sys.argv[1:]
        if isinstance(args, str):
            args = args.split(" ")

        parser = argparse.ArgumentParser(
            description="Convert exposure-based data to source-based data",
            usage="python -m jwst.exp_to_source.main files",
        )
        parser.add_argument("files", type=str, nargs="+", help="Files to convert")

        parser.add_argument(
            "-o",
            "--output-path",
            type=str,
            default=".",
            help='Folder to save results in. Default: "%(default)s"',
        )

        parser.add_argument(
            "--dry-run",
            action="store_true",
            dest="dry_run",
            help="Execute but do not save results.",
        )

        try:
            parsed = parser.parse_args(args=args)
        except SystemExit:
            return

        # self.sources is a dict keyed on source name whose value is
        # the corresponding MultiExposureModel,
        # e.g. {source_name: MultiExposureModel, ...}
        exposures = [MultiSlitModel(f) for f in parsed.files]
        self.sources = exp_to_source(exposures)
        if not parsed.dry_run:
            for source in self.sources:
                out_path = ".".join([source, "fits"])
                out_path = Path(parsed.output_path) / out_path
                self.sources[source].save(out_path)


if __name__ == "__main__":
    Main()
