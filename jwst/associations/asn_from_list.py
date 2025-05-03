"""Create an association from a list."""

import argparse
import sys
from pathlib import Path

from . import AssociationRegistry
from .lib.rules_level3_base import DMS_Level3_Base

__all__ = ["asn_from_list"]


def asn_from_list(items, rule=DMS_Level3_Base, **kwargs):
    """
    Create an association from a list.

    Parameters
    ----------
    items : [object [, ...]]
        List of items to add.
    rule : `Association` rule
        The association rule to use.
    **kwargs : dict
        Other named parameters required or pertinent to adding
        the items to the association.

    Returns
    -------
    association : `Association`-based instance
        The association with the items added.

    Notes
    -----
    This is a lower-level tool for artificially creating
    an association. As such, the association created may not be valid.
    It is presumed the user knows what they are doing.
    """
    asn = rule()
    asn._add_items(items, **kwargs)  # noqa: SLF001
    return asn


class Main:
    """Command-line interface for asn_from_list."""

    def __init__(self, args=None):
        """
        Initialize the class for asn_from_list.

        Parameters
        ----------
        args : list or None
            The command line arguments. Can be one of
                - `None`: `sys.argv` is then used.
                - `[str, ...]`: A list of strings which create the command line
                  with the similar structure as `sys.argv`
        """
        self.configure(args)

    def asn_from_list(self):
        """Create the associaton from the list."""
        parsed = self.parsed
        self.asn = asn_from_list(
            parsed.filelist, rule=self.rule, product_name=parsed.product_name, acid=parsed.acid
        )

    @classmethod
    def cli(cls, args=None):
        """
        Command-line interface to creating an association from a list.

        Parameters
        ----------
        args : list or None
            The command line arguments. Can be one of
                - `None`: `sys.argv` is then used.
                - `[str, ...]`: A list of strings which create the command line
                  with the similar structure as `sys.argv`

        Returns
        -------
        ~jwst.associations.association.Association
            The association created from the list.
        """
        association = cls(args=args)
        association.asn_from_list()
        association.save()
        return association

    def configure(self, args=None):
        """
        Configure the asn_from_list CLI object.

        Parameters
        ----------
        args : list or None
            The command line arguments. Can be one of
                - `None`: `sys.argv` is then used.
                - `[str, ...]`: A list of strings which create the command line
                  with the similar structure as `sys.argv`
        """
        if args is None:
            args = sys.argv[1:]
        if isinstance(args, str):
            args = args.split(" ")

        parser = argparse.ArgumentParser(
            description="Create an association from a list of files",
            usage="asn_from_list -o mosaic_asn.json --product-name my_mosaic *.fits",
        )

        parser.add_argument(
            "-o", "--output-file", type=str, required=True, help="File to write association to"
        )

        parser.add_argument(
            "-f",
            "--format",
            type=str,
            default="json",
            help='Format of the association files. Default: "%(default)s"',
        )

        parser.add_argument(
            "--product-name", type=str, help="The product name when creating a Level 3 association"
        )

        parser.add_argument(
            "-r",
            "--rule",
            type=str,
            default="DMS_Level3_Base",
            help=('The rule to base the association structure on. Default: "%(default)s"'),
        )
        parser.add_argument(
            "--ruledefs",
            action="append",
            help=(
                "Association rules definition file(s)"
                " If not specified, the default rules will be searched."
            ),
        )
        parser.add_argument(
            "-i",
            "--id",
            type=str,
            default="o999",
            help='The association candidate id to use. Default: "%(default)s"',
            dest="acid",
        )

        parser.add_argument(
            "filelist", type=str, nargs="+", help="File list to include in the association"
        )

        self.parsed = parsed = parser.parse_args(args=args)

        # Get the rule
        self.rule = AssociationRegistry(parsed.ruledefs, include_bases=True)[parsed.rule]

    def save(self):
        """Save association."""
        parsed = self.parsed
        with Path(parsed.output_file).open("w") as outfile:
            name, serialized = self.asn.dump(format=parsed.format)
            outfile.write(serialized)


def main(args=None):
    """
    Command-line entrypoint for asn_from_list.

    Wrapper around `Main.cli` so that the return is either True or an exception

    Parameters
    ----------
    args : list or None
        The command line arguments. Can be one of

        - `None`: `sys.argv` is then used.
        - `[str, ...]`: A list of strings which create the command line
          with the similar structure as `sys.argv`
    """
    Main.cli(args=args)
