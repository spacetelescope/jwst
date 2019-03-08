#!/usr/bin/env python

import sys

from jwst import jwpsf


def run():
    jwpsf.jwpsf()

if __name__ == '__main__':
    if '--version' in sys.argv:
        print(jwpsf.__version__)
        sys.exit(0)
    run()
