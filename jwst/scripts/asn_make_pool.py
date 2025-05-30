#!/usr/bin/env python

from jwst.associations import mkpool


def main():
    """
    Create an association pool from a list of FITS files.

    Implement the command line functionality of `mkpool`. For full information on
    the associtation tool use `asn_make_pool --help`.

    Returns
    -------
    pool : jwst.associations.AssociationPool
       The association pool.
    """
    kwargs = mkpool.from_cmdline()
    pool_file = kwargs["pool"]
    del kwargs["pool"]

    pool = mkpool.mkpool(**kwargs)

    pool.write(pool_file, overwrite=True)


if __name__ == "__main__":
    main()
