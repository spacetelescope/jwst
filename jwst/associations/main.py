"""Main entry for the association generator."""

import sys
import argparse
import logging

import numpy as np

from pathlib import Path
from jwst import __version__
from jwst.associations import AssociationPool
from jwst.associations import config
from jwst.associations.exceptions import AssociationError
from jwst.associations.lib.log_config import log_config, DMS_config

__all__ = ["Main", "main"]

# Configure logging
logger = log_config(name=__package__)


class Main:
    """
    Generate Associations from an Association Pool.

    Parameters
    ----------
    args : [str, ...], or None
        The command line arguments. Can be one of
        - `None`: `sys.argv` is then used.
        - `[str, ...]`: A list of strings which create the command line
        with the similar structure as `sys.argv`

    pool : None or AssociationPool
        If `None`, a pool file must be specified in the `args`.
        Otherwise, an `AssociationPool`

    Attributes
    ----------
    pool : `AssociationPool`
        The pool read in, or passed in through the parameter `pool`
    rules : `AssociationRegistry`
        The rules used for association creation.
    associations : [`Association`, ...]
        The list of generated associations.

    Notes
    -----
    Refer to the :ref:`Association Generator <associations>`
    documentation for a full description.
    """

    def __init__(self, args=None, pool=None):
        self.configure(args=args, pool=pool)

    @classmethod
    def cli(cls, args=None, pool=None):
        """
        Run the full association generation process.

        Parameters
        ----------
        args : [str, ...], or None
            The command line arguments. Can be one of
            - `None`: `sys.argv` is then used.
            - `[str, ...]`: A list of strings which create the command line
            with the similar structure as `sys.argv`

        pool : None or AssociationPool
            If `None`, a pool file must be specified in the `args`.
            Otherwise, an `AssociationPool`

        Returns
        -------
        generator : Main
            A fully executed association generator.
        """
        generator_cli = cls(args=args, pool=pool)
        generator_cli.generate()
        generator_cli.save()
        return generator_cli

    @property
    def orphaned(self):
        """
        The pool of exposures that do not belong to any association.

        Returns
        -------
        pool : `jwst.associations.pool.AssociationPool`
            The pool of orphaned exposures.
        """
        not_in_asn = np.ones((len(self.pool),), dtype=bool)
        for asn in self.associations:
            try:
                indexes = [item.index for item in asn.from_items]
            except AttributeError:
                continue
            not_in_asn[indexes] = False

        orphaned = self.pool[not_in_asn]
        return orphaned

    def configure(self, args=None, pool=None):
        """
        Configure to prepare for generation.

        Parameters
        ----------
        args : [str, ...], or None
            The command line arguments. Can be one of

            - `None`: `sys.argv` is then used.
            - `[str, ...]`: A list of strings which create the command line
              with the similar structure as `sys.argv`

        pool : None or AssociationPool
            If `None`, a pool file must be specified in the `args`.
            Otherwise, an `AssociationPool`
        """
        self.parse_args(args, has_pool=pool)
        parsed = self.parsed

        # Configure logging
        logging_config = None
        if parsed.DMS_enabled:
            logging_config = DMS_config
        logger = log_config(name=__package__, config=logging_config)
        logger.setLevel(parsed.loglevel)
        config.DEBUG = (parsed.loglevel != 0) and (parsed.loglevel <= logging.DEBUG)

        # Preamble
        logger.info("Command-line arguments: %s", parsed)
        logger.context.set("asn_candidate_ids", parsed.asn_candidate_ids)

        if pool is None:
            logger.info(f"Reading pool {parsed.pool}")
            pool = AssociationPool.read(
                parsed.pool,
                delimiter=parsed.delimiter,
                fmt=parsed.pool_format,
            )
        self.pool = pool

        # DMS: Add further info to logging.
        try:
            logger.context.set("program", self.pool[0]["PROGRAM"])
        except KeyError:
            pass

        # Determine mode of operation. Options are
        #  1) Only specified candidates
        #  2) Only discovered associations that do not match
        #     candidate associations
        #  3) Both discovered and all candidate associations.
        if not parsed.discover and not parsed.all_candidates and parsed.asn_candidate_ids is None:
            parsed.discover = True
            parsed.all_candidates = True

    def generate(self):
        """Generate the associations."""
        logger.info("Generating associations.")
        parsed = self.parsed
        if parsed.per_pool_algorithm:
            from jwst.associations.generator.generate_per_pool import generate_per_pool

            self.associations = generate_per_pool(
                self.pool,
                rule_defs=parsed.rules,
                candidate_ids=parsed.asn_candidate_ids,
                all_candidates=parsed.all_candidates,
                discover=parsed.discover,
                version_id=parsed.version_id,
                finalize=not parsed.no_finalize,
                merge=parsed.merge,
                ignore_default=parsed.ignore_default,
            )
        else:
            from jwst.associations.generator.generate_per_candidate import generate_per_candidate

            self.associations = generate_per_candidate(
                self.pool,
                rule_defs=parsed.rules,
                candidate_ids=parsed.asn_candidate_ids,
                all_candidates=parsed.all_candidates,
                discover=parsed.discover,
                version_id=parsed.version_id,
                finalize=not parsed.no_finalize,
                merge=parsed.merge,
                ignore_default=parsed.ignore_default,
                dms_enabled=parsed.DMS_enabled,
            )

        logger.debug(self.__str__())

    def parse_args(self, args=None, has_pool=False):
        """
        Set command line arguments.

        Parameters
        ----------
        args : list, str, or None
            List of command-line arguments.
            If a string, spaces separate the arguments.
            If None, `sys.argv` is used.

        has_pool : bool-like
            Do not require `pool` from the command line if a pool is already in hand.
        """
        if args is None:
            args = sys.argv[1:]
        if isinstance(args, str):
            args = args.split(" ")

        parser = argparse.ArgumentParser(
            description="Generate Association Data Products", usage="asn_generate pool"
        )
        if not has_pool:
            parser.add_argument("pool", type=str, help="Association Pool")
        op_group = parser.add_mutually_exclusive_group()
        op_group.add_argument(
            "-i",
            "--ids",
            nargs="+",
            dest="asn_candidate_ids",
            help="space-separated list of association candidate IDs to operate on.",
        )
        op_group.add_argument(
            "--discover", action="store_true", help="Produce discovered associations"
        )
        op_group.add_argument(
            "--all-candidates",
            action="store_true",
            dest="all_candidates",
            help="Produce all association candidate-specific associations",
        )
        parser.add_argument(
            "-p",
            "--path",
            type=Path,
            default=".",
            help='Folder to save the associations to. Default: "%(default)s"',
        )
        parser.add_argument(
            "--save-orphans",
            dest="save_orphans",
            nargs="?",
            const=Path("orphaned.csv"),
            default=False,
            help='Save orphaned items into the specified table. Default: "%(default)s"',
        )
        parser.add_argument(
            "--version-id",
            dest="version_id",
            nargs="?",
            const=True,
            default=None,
            help=(
                "Version tag to add into association name and products."
                " If not specified, no version will be used."
                " If specified without a value, the current time is used."
                " Otherwise, the specified string will be used."
            ),
        )
        parser.add_argument("-r", "--rules", action="append", help="Association Rules file.")
        parser.add_argument(
            "--ignore-default",
            action="store_true",
            help="Do not include default rules. -r should be used if set.",
        )
        parser.add_argument(
            "--dry-run",
            action="store_true",
            dest="dry_run",
            help="Execute but do not save results.",
        )
        parser.add_argument(
            "-d",
            "--delimiter",
            type=str,
            default="|",
            help="""Delimiter
            to use if pool files are comma-separated-value
            (csv) type files. Default: "%(default)s"
            """,
        )
        parser.add_argument(
            "--pool-format",
            type=str,
            default="ascii",
            help=(
                "Format of the pool file."
                " Any format allowed by the astropy"
                " Unified File I/O interface is allowed."
                ' Default: "%(default)s"'
            ),
        )
        parser.add_argument(
            "-v",
            "--verbose",
            action="store_const",
            dest="loglevel",
            const=logging.INFO,
            default=logging.NOTSET,
            help="Output progress and results.",
        )
        parser.add_argument(
            "-D",
            "--debug",
            action="store_const",
            dest="loglevel",
            const=logging.DEBUG,
            help="Output detailed debugging information.",
        )
        parser.add_argument(
            "--DMS",
            action="store_true",
            dest="DMS_enabled",
            help="Running under DMS workflow conditions.",
        )
        parser.add_argument(
            "--format",
            default="json",
            help='Format of the association files. Default: "%(default)s"',
        )
        parser.add_argument(
            "--version",
            action="version",
            version=f"%(prog)s {__version__}",
            help="Version of the generator.",
        )
        parser.add_argument(
            "--no-finalize",
            action="store_true",
            help="Do not run the finalization methods on the interim associations",
        )
        parser.add_argument(
            "--merge",
            action="store_true",
            help="Merge associations into single associations with multiple products",
        )
        parser.add_argument(
            "--no-merge",
            action=DeprecateNoMerge,
            help='Deprecated: Default is to not merge. See "--merge".',
        )
        parser.add_argument(
            "--per-pool-algorithm",
            action="store_true",
            help="Use the original, per-pool, algorithm that does not "
            "segment pools based on candidates",
        )

        self.parsed = parser.parse_args(args=args)

    def save(self):
        """Save the associations to disk."""
        if self.parsed.dry_run:
            return

        for asn in self.associations:
            try:
                (fname, serialized) = asn.dump(format=self.parsed.format)
            except AssociationError as exception:
                logger.warning("Cannot serialize association %s", asn)
                logger.warning("Reason:", exc_info=exception)
                continue
            with Path(self.parsed.path / Path(fname)).open("w") as f:
                f.write(serialized)

        if self.parsed.save_orphans:
            self.orphaned.write(
                self.parsed.path / self.parsed.save_orphans,
                format="ascii",
                delimiter="|",
            )

    def __str__(self):
        result = []
        result.append(
            f"There were {len(self.associations):d} associations "
            f"and {len(self.orphaned):d} orphaned items found.\n"
            "Associations found are:"
        )
        for assocs in self.associations:
            result.append(assocs.__str__())

        return "\n".join(result)


def main(args=None, pool=None):
    """
    Command-line entrypoint for the association generator.

    Wrapper around `Main.cli` so that the return is either True or an exception.

    Parameters
    ----------
    args : [str, ...], or None
        The command line arguments. Can be one of

        - `None`: `sys.argv` is then used.
        - `[str, ...]`: A list of strings which create the command line
        with the similar structure as `sys.argv`

    pool : None or AssociationPool
        If `None`, a pool file must be specified in the `args`.
        Otherwise, an `AssociationPool`
    """
    Main.cli(args, pool)


# #########
# Utilities
# #########
class DeprecateNoMerge(argparse.Action):
    """Deprecate the `--no-merge` option."""

    def __init__(self, option_strings, dest, **kwargs):
        super(DeprecateNoMerge, self).__init__(option_strings, dest, const=True, nargs=0, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):  # noqa: ARG002
        logger.warning(
            'The "--no-merge" option is now the default and deprecated.'
            ' Use "--merge" to force merging.'
        )
        setattr(namespace, self.dest, values)


if __name__ == "__main__":
    Main()
