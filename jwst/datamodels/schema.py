# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, unicode_literals, print_function


from astropy.extern import six
from collections import OrderedDict


# return_result included for backward compatibility
def find_fits_keyword(schema, keyword, return_result=False):
    """
    Utility function to find a reference to a FITS keyword in a given
    schema.  This is intended for interactive use, and not for use
    within library code.

    Parameters
    ----------
    schema : JSON schema fragment
        The schema in which to search.

    keyword : str
        A FITS keyword name

    Returns
    -------
    locations : list of str
    """
    def find_fits_keyword(subschema, path, combiner, ctx, recurse):
        if len(path) and path[0] == 'extra_fits':
            return True
        if subschema.get('fits_keyword') == keyword:
            results.append('.'.join(path))

    results = []
    walk_schema(schema, find_fits_keyword, results)

    return results

def build_fits_dict(schema):
    """
    Utility function to create a dict that maps FITS keywords to their
    metadata attribute in a input schema.  
    
    Parameters
    ----------
    schema : JSON schema fragment
        The schema in which to search.
 
    Returns
    -------
    results : dict
        Dictionary with FITS keywords as keys and schema metadata 
        attributes as values
   
    """
    def build_fits_dict(subschema, path, combiner, ctx, recurse):
        if len(path) and path[0] == 'extra_fits':
            return True
        kw = subschema.get('fits_keyword')
        if kw is not None:
            results[kw] = '.'.join(path)

    results = {}
    walk_schema(schema, build_fits_dict, results)

    return results


class SearchSchemaResults(list):
    def __repr__(self):
        import textwrap

        result = []
        for path, description in self:
            result.append(path)
            result.append(
                textwrap.fill(
                    description, initial_indent='    ',
                    subsequent_indent='    '))
        return '\n'.join(result)


def search_schema(schema, substring):
    """
    Utility function to search the metadata schema for a particular
    phrase.

    This is intended for interactive use, and not for use within
    library code.

    The searching is case insensitive.

    Parameters
    ----------
    schema : JSON schema fragment
        The schema in which to search.

    substring : str
        The substring to search for.

    Returns
    -------
    locations : list of tuples
    """
    substring = substring.lower()

    def find_substring(subschema, path, combiner, ctx, recurse):
        matches = False
        for param in ('title', 'description'):
            if substring in schema.get(param, '').lower():
                matches = True
                break

        if substring in '.'.join(path).lower():
            matches = True

        if matches:
            description = '\n\n'.join([
                schema.get('title', ''),
                schema.get('description', '')]).strip()
            results.append(('.'.join(path), description))

    results = SearchSchemaResults()
    walk_schema(schema, find_substring, results)
    results.sort()
    return results


def walk_schema(schema, callback, ctx={}):
    """
    Walks a JSON schema tree in breadth-first order, calling a
    callback function at each entry.

    Parameters
    ----------
    schema : JSON schema

    callback : callable
        The callback receives the following arguments at each entry:

        - subschema: The subschema for the entry
        - path: A list of strings defining the path to the entry
        - combiner: The current combiner in effect, will be 'allOf',
          'anyOf', 'oneOf', 'not' or None
        - ctx: An arbitrary context object, usually a dictionary
        - recurse: A function to call to recurse deeper on a node.

        If the callback returns `True`, the subschema will not be
        further recursed.

    ctx : object, optional
        An arbitrary context object
    """
    def recurse(schema, path, combiner, ctx):
        if callback(schema, path, combiner, ctx, recurse):
            return

        for c in ['allOf', 'not']:
            for sub in schema.get(c, []):
                recurse(sub, path, c, ctx)

        for c in ['anyOf', 'oneOf']:
            for i, sub in enumerate(schema.get(c, [])):
                recurse(sub, path + [i], c, ctx)

        if schema.get('type') == 'object':
            for key, val in six.iteritems(schema.get('properties', {})):
                recurse(val, path + [key], combiner, ctx)

        if schema.get('type') == 'array':
            items = schema.get('items', {})
            if isinstance(items, list):
                for i, item in enumerate(items):
                    recurse(item, path + [i], combiner, ctx)
            elif len(items):
                recurse(items, path + ['items'], combiner, ctx)

    recurse(schema, [], None, ctx)


def flatten_combiners(schema):
    """
    Flattens the allOf and anyOf operations in a JSON schema.

    TODO: Write caveats -- there's a lot
    """
    newschema = OrderedDict()

    def add_entry(path, schema, combiner):
        # TODO: Simplify?
        cursor = newschema
        for i in range(len(path)):
            part = path[i]
            if isinstance(part, int):
                cursor = cursor.setdefault('items', [])
                while len(cursor) <= part:
                    cursor.append({})
                cursor = cursor[part]
            elif part == 'items':
                cursor = cursor.setdefault('items', OrderedDict())
            else:
                cursor = cursor.setdefault('properties', OrderedDict())
                if i < len(path) - 1 and isinstance(path[i + 1], int):
                    cursor = cursor.setdefault(part, [])
                else:
                    cursor = cursor.setdefault(part, OrderedDict())

        cursor.update(schema)

    def callback(schema, path, combiner, ctx, recurse):
        type = schema.get('type')
        schema = OrderedDict(schema)
        if type == 'object':
            del schema['properties']
        elif type == 'array':
            del schema['items']
        if 'allOf' in schema:
            del schema['allOf']
        if 'anyOf' in schema:
            del schema['anyOf']

        add_entry(path, schema, combiner)

    walk_schema(schema, callback)

    return newschema
