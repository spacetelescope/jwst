#!/usr/bin/env python

import argparse

from jwst.associations import asn_edit


def main():
    """
    Parse command line, read, edit, write association file
    """

    # Parse command line arguments
    description_text = """
Edit Association File

This script adds or removes filenames from an association file. The
first argument is the name of the association file. Subsequent
arguments are the filenames. Options determine which operation is
performed: --add or --remove. If adding files the --type option sets
the exposure type for the new files. If removing file, the --ignore
option will not use the filename suffix when matching the filenames for
removal. Normally the output association file is the same as the input.
The --output option allows you to set a different filename. All options
can be abbreviated to their first letter.

When adding files to an association, the file must exist on the disk.
When removing files, the filename must be found in the association. If
not, no change will be made to the file.
    """
    parser = argparse.ArgumentParser(
        description=description_text,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    arg_group = parser.add_mutually_exclusive_group()
    arg_group.add_argument(
        '-a', '--add', action='store_true', help='Add files to association'
    )
    arg_group.add_argument(
        '-r', '--remove', action='store_true', help='Remove files to association'
    )

    parser.add_argument(
        '-t', '--type', default='science', help='Exptype, if adding filenames'
    )
    parser.add_argument(
        '-i',
        '--ignore',
        action='store_true',
        help='Ignore suffix on filename when matching',
    )

    parser.add_argument(
        '-o', '--output', help='Output association name if different than input'
    )

    parser.add_argument('association', help='The association file name')

    parser.add_argument('filenames', nargs='+', help='The filenames to process')

    args = parser.parse_args()

    # Read the association into memory
    asn = asn_edit.reader(args.association)

    # Either add or remove filenames
    if args.add:
        asn_edit.add(asn, args.filenames, args.type)
    elif args.remove:
        asn_edit.remove(asn, args.filenames, args.ignore)

    # Write the edited file out again
    if args.output:
        output_file = args.output
    else:
        output_file = args.association
    asn_edit.writer(asn, output_file)


if __name__ == '__main__':
    main()
