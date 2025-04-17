#!/usr/bin/env python

import sys
import logging
from pathlib import Path

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs


logger = logging.getLogger()
logger.setLevel(logging.INFO)


def main():
    """
    Get copies of all the cfg files currently in use by the jwst pipeline software.

    Use from terminal as follows:
    $ collect_pipeline_cfgs .
    """
    if len(sys.argv) < 2:
        logging.error("ERROR: missing argument (destination directory")
        sys.exit(1)

    if len(sys.argv) > 2:
        logging.error("ERROR: too many arguments")
        sys.exit(1)

    dst = sys.argv[1]

    if Path(dst).exists():
        # if the destination directory already exists, make sure the user wants to
        # overwrite its contents
        count = 0
        while True:
            msg = (
                f"WARNING: {dst} already exists, do you want to overwrite the \n"
                f"contents of this directory?  Enter yes or no [yes]: "
            )

            answer = input(msg).strip().lower()
            if answer not in ["yes", "no", ""]:
                logging.info(f"\n{answer} is not a valid response\n")
            else:
                if answer in ["yes", ""]:
                    break
                else:
                    sys.exit(1)

            # exit after 10 tries
            count += 1
            if count > 10:
                sys.exit(1)

    collect_pipeline_cfgs(dst)


if __name__ == "__main__":
    main()
