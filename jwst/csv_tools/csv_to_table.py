import re
from astropy.io import ascii


def csv_to_table(handle,
                 delimiter='|',
                 comments=r'\s*(#|//|--)',
                 header_search=r'(?P<keyword>\w+)\s*=\s*(?P<value>\w+.*?)(/(?P<comment>.*))?$'):
    '''Produce an astropy table, and related keyword/value pairs from a CSV file.

    Parameters
    ----------
    handle: File handle
            The file object being read from.

    delimiter: str
               The single character delimiter that distinguish columns.

    comments: string
              A regular expression for determining whether a line
              is a comment or not.

    header_search: string
                   A regular expression the finds and parses
                   keyword/value pairs from comment lines.
                   The regex must contain two named groups:
                   'keyword' for the keyword match and 'value'
                   for the value. An optional 'comment'
                   group may also be present to put a comment with
                   the keyword/value pair, as in standard FITS fashion.
    Returns
    -------
    astropy.Table:
        An astropy Table is returned.
        Of note is one extra parameter on the table, 'meta',
        which is the dict of header keywords found.
    '''

    # Setup the astropy reader.
    reader = ascii.get_reader(Reader=CSVKeywords)
    reader.header.start_line = 0
    reader.header.splitter.delimiter = delimiter
    reader.data.splitter.delimiter = reader.header.splitter.delimiter
    reader.header.comment = comments
    reader.data.comment = reader.header.comment
    reader.header.keywords = header_search

    # All setup, return the table
    return reader.read(handle)


class CSVHeaderKeywords(ascii.BasicHeader):
    '''Class to search the csv file for keywords.'''

    def update_meta(self, lines, meta):
        pattern = re.compile(r'\s+'.join((self.comment, self.keywords)))
        for line in lines:
            matches = pattern.match(line)
            if matches is not None:
                groups = matches.groupdict()
                value = groups['value'].strip()
                if 'comment' in groups and groups['comment'] is not None:
                    value = (value, groups['comment'].strip())
                meta['table'][groups['keyword']] = value


class CSVKeywords(ascii.Basic):
    header_class = CSVHeaderKeywords
