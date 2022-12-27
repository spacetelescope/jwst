#!/usr/bin/env python

import os
import sys

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs


def main():
    if len(sys.argv) < 2:
        print('ERROR: missing argument (destination directory')
        sys.exit(1)

    if len(sys.argv) > 2:
        print('ERROR: too many arguments')
        sys.exit(1)

    dst = sys.argv[1]

    if os.path.exists(dst):
        # if the destination directory already exists, make sure the user wants to
        # overwrite it's contents
        count = 0
        while True:
            msg = (
                f'WARNING: {dst} already exists, do you want to overwrite the \n'
                f'contents of this directory?  Enter yes or no [yes]: '
            )

            sys_vrs = sys.version[:3]
            if sys_vrs == '2.7':
                answer = raw_input(msg).strip().lower()  # noqa
            else:
                answer = input(msg).strip().lower()
            if answer not in ['yes', 'no', '']:
                print('%s is not a valid response\n' % answer)
            else:
                print('')
                if answer in ['yes', '']:
                    break
                else:
                    sys.exit(1)

            # exit after 10 tries
            count += 1
            if count > 10:
                print()
                sys.exit(1)

    collect_pipeline_cfgs(dst)


if __name__ == '__main__':
    main()
