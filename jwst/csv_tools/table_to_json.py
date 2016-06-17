from collections import OrderedDict
import json

def table_to_json(table, **json_kwargs):
    """Create a JSON string based on the table

    Parameters
    ----------
    table: astropy.table.Table
           The table. Note that the extra parameter, table.meta,
           if present, is expected to be a dictionary of keyword/value pairs.

    json_kwargs: dict
                 Arguments to the json.dumps call.

    Notes
    -----
    The JSON will consist of an object where keyword/value pairs are located
    in the header object and the table will be in the table object.

    Returns
    -------
    str: JSON
    """

    data = {}

    try:
        data['header'] = table.meta
    except AttributeError:
        pass

    data['columns'] = OrderedDict((name, list(table[name])) for name in table.colnames)

    return json.dumps(data, **json_kwargs)
