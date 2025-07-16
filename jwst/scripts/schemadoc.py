#!/usr/bin/env python

"""Create documentation from the schema file of a datamodel class."""

# Licensed under a 3-clause BSD style license - see LICENSE

import argparse
import logging
from pathlib import Path

from stdatamodels.jwst.datamodels import _defined_models as defined_models
from stdatamodels.schema import build_docstring


def get_docstrings(template, model_names, all_models=False):
    """Get the docstring for every model class."""
    logger = logging.getLogger()
    if all_models:
        klasses = defined_models
    else:
        klasses = {}
        for model_name in model_names:
            klasses[model_name] = defined_models[model_name]

    for klassname, klass in klasses.items():
        try:
            docstring = build_docstring(klass, template)
        except Exception as err:
            logger.error(f"{klassname} : {err}")
        else:
            logger.info(f".. {klassname} ..")
            logger.info(docstring)


def main():
    """Get docstrings from a specific datamodel or all datamodel(s)."""
    long_description = """
    Create documentation from the schema file of a datamodel class.
    """
    parser = argparse.ArgumentParser(description=long_description)
    parser.add_argument(
        "-a", "--all_models", action="store_true", help="generate docstring for all models"
    )
    parser.add_argument("-t", "--template", type=str, help="input file containing templates")
    parser.add_argument("models", nargs="*", help="Models to generate docstrings for")
    args = parser.parse_args()

    # Set logger to only print to screen
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    if args.template is None:
        template = "{path} : {title} ({datatype})\n"
    else:
        with Path.open(args.template, "r") as fd:
            template = fd.readlines()
            template = "".join(template)

    get_docstrings(template, args.models, all_models=args.all_models)


if __name__ == "__main__":
    main()
