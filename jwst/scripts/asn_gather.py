#!/usr/bin/env python

"""Copy or Move data that is listed in an association."""

from jwst.associations import asn_gather

if __name__ == "__main__":
    kwargs = asn_gather.from_cmdline()
    asn_gather.asn_gather(**kwargs)
