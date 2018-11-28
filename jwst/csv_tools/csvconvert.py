import sys
import argparse
import inspect

from .csv_to_table import csv_to_table
from .table_to_hdulist import table_to_hdulist
from .table_to_json import table_to_json


class CSVConvertScript():
    """Convert a CSV file to other formats
    version of the CSV file.
    """

    formats = ['json', 'fits']
    json_separators = (',', ': ')
    json_indent = 4

    def __init__(self, args=None):

        if args is None:
            args = sys.argv[1:]

        if isinstance(args, str):
            args = args.split(' ')

        # Get some defaults from other routines
        defaults = get_default_args(csv_to_table)

        # Define the command line.
        parser = argparse.ArgumentParser(description='Convert CSV with headers to other formats',
                                         usage='python -m jwst_tools.csv_tools.csvconvert csvfile outfile')
        parser.add_argument('csvfile', type=str,
                            help='Input CSV file to convert')
        parser.add_argument('outfile', type=str,
                            help='Output file')
        parser.add_argument('-f', '--format', type=str,
                            default='json',
                            help='Output format. Default is "%(default)s". Options are: {}'.format(self.formats))
        parser.add_argument('-c', '--comments', type=str,
                            default=defaults['comments'],
                            help='Regular expression of initial characters that indicate comment lines. Default: "%(default)s"')
        parser.add_argument('-s', '--header-search', type=str,
                            default=defaults['header_search'],
                            help='Regular expression to pull the keyword and value from a comment. Default: "%(default)s"')
        parser.add_argument('-d', '--delimiter', type=str,
                            default=defaults['delimiter'],
                            help='Single character delimiter to distinguish columns. Default: "%(default)s"')
        parsed = parser.parse_args(args=args)

        # Build the named parameter arguments
        kwargs = {}
        kwargs['comments'] = parsed.comments
        kwargs['header_search'] = parsed.header_search
        kwargs['delimiter'] = parsed.delimiter

        # Do the conversion
        self.table = csv_to_table(parsed.csvfile, **kwargs)
        self.format = None
        if parsed.format == 'json':
            self.output = table_to_json(self.table,
                                        separators=self.json_separators,
                                        indent=self.json_indent)
            with open(parsed.outfile, 'w') as outfile:
                outfile.write(self.output)

        elif parsed.format == 'fits':
            self.output = table_to_hdulist(self.table)
            self.output.writeto(parsed.outfile, clobber=True)
        else:
            raise NotImplementedError('Format {} is not implemented.'.format(parsed.format))
        self.format = parsed.format

#*********
# Utilities
#*********

def get_default_args(func):
    '''Retrieve the defaults for arguments of a given fucntion.'''

    args = inspect.getargspec(func)
    return dict(list(zip(args.args[-len(args.defaults):], args.defaults)))


#******
# Main
#******

if __name__ == '__main__':
    CSVConvertScript()
