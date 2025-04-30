#!/usr/bin/env python

from jwst.associations import asn_gather


def main():
    """Copy or move data that is listed in an association."""
    kwargs = asn_gather.from_cmdline()
    asn_gather.asn_gather(**kwargs)


if __name__ == "__main__":
    main()
