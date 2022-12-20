#!/usr/bin/env python

from jwst.associations import mkpool


def main():
    kwargs = mkpool.from_cmdline()
    pool_file = kwargs['pool']
    del kwargs['pool']

    pool = mkpool.mkpool(**kwargs)

    pool.write(pool_file, overwrite=True)


if __name__ == '__main__':
    main()
