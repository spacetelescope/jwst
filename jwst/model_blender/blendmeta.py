""" blendmeta - Merge metadata from multiple models to create
                   a new metadata instance and table

"""
from collections import OrderedDict

import numpy as np
from astropy.io import fits

from .. import datamodels
from .. import associations
from ..stpipe import crds_client
from ..datamodels import fits_support

from .blendrules import KeywordRules

__version__ = '0.8.0'
__vdate__ = '14-Feb-2018'

empty_list = [None, '', ' ', 'INDEF', 'None']


# Primary functional interface for the code
def blendmodels(product, inputs=None, output=None,
                rules_file=None, verbose=False):
    """ Blend models that went into creating the original drzfile into a
    new metadata instance with a table that contains attribute values from
    all input datamodels.

    The drzfile will be used to determine the names of the input models,
    should no filenames be provided in the 'inputs' parameter.

    The drzfile will be updated 'in-place' with the new metadata attributes
    and table.

    Parameters
    ----------
    product : str
        Name of combined product with metadata that needs updating. This can
        be specified as a single filename.
        When no value for 'inputs' has been provided, this file
        will be used to determine the names of the input datamodels
        whose metadata need to be blended to create the new combined metadata.

    inputs : list, optional
        If provided, the filenames provided in this list will be used to get
        the metadata which will be blended into the final output metadata.

    rules_file : str, optional [Default: None]
        Filename of a rules file to be used.  If None, look for default file
        specified in the input models metadata as the 'BLENDRULES' reference
        file.

    verbose : bool, optional [Default: False]
        Print out additional messages during processing when specified.

    """

    if inputs in empty_list:
        inputs = extract_filenames_from_product(product)
        input_filenames = inputs
    else:
        input_filenames = [i.meta.filename for i in inputs]

    if verbose:
        print('Creating blended metadata from: ')
        for i in input_filenames:
            print('\t{}'.format(i))

    newmeta, newtab = get_blended_metadata(inputs, verbose=verbose,
                                           rules_file=rules_file)

    # open product to update metadata from input models used to create product
    output_model = datamodels.open(product)
    output_model.meta = newmeta
    """
        # This is code used by resample.blend to merge blended attributes
        # into output_model metadata.
        #
        # Now merge the keyword values from new_hdrs into the metatdata for the
        # output datamodel
        #
        # start by building dict which maps all FITS keywords in schema to
        # attribute in the schema
        fits_dict = schema.build_fits_dict(output_model.schema)
        # Need to insure that output_model does not already have an instance
        # of hdrtab from previous processing, an instance that would be
        # incompatible with the new table generated now...
        if hasattr(output_model, 'hdrtab'):
            # If found, remove it to be replaced by new instance
            del output_model.hdrtab
        # Now assign values from new_hdrs to output_model.meta using
        #  fits_dict map
        for hdr in new_hdrs:
            for kw in hdr:
                if kw in fits_dict:
                    output_model[fits_dict[kw]] = hdr[kw]
    """
    # Now, append HDRTAB as new element in datamodel
    newtab_schema = build_tab_schema(newtab)
    output_model.add_schema_entry('hdrtab', newtab_schema)
    output_model.hdrtab = fits_support.from_fits_hdu(newtab, newtab_schema)

    print('Updated metadata in ', output, ' with blended metadata.')

    # Clean up for the next run
    del newmeta, newtab


def get_blended_metadata(input_models, verbose=False, rules_file=None):
    """
    Return a blended metadata instance and table based on the input datamodels.
    This will serve as the primary interface for blending datamodels.

    Parameters
    ----------
    input_models : list
        Either a single list of filenames from which to extract the metadata to
        be blended, or a list of `datamodels.DataModel` objects to be blended.

    rules_file : str, optional [Default: None]
        Filename of a rules file to be used.  If None, look for default file
        installed with the code.

    Returns
    -------
    metadata : list
        A list of blended metadata instances, one for each i

    new_table : object
        Single fits.TableHDU object that contains the combined results from
        all input headers(extension). Each row will correspond to an image,
        and each column corresponds to a single keyword listed in the rules.

    """
    hdrlist = None
    if not isinstance(input_models, list):
        input_models = [input_models]

    # Look for rules reference file from CRDS
    if rules_file is None:
        rules_file = crds_client.get_reference_file(input_models[0],
                                                    'blendrules')

    # Turn input filenames into a set of header objects
    if isinstance(input_models[0], str):
        hdrlist = []
        for fname in input_models:
            model = datamodels.open(fname)
            hdrlist.append(model.meta.copy())
    else:
        hdrlist = input_models

    num_files = len(hdrlist)

    # Determine what blending rules need to be merged to create the final
    # blended headers. There will be a separate set of rules for each
    # instrument, and all rules get merged into a composite set of rules that
    # get applied to all input headers regardless of instrument.
    #
    # Instrument identification will be extracted from the INSTRUME keyword
    # from the PRIMARY header of each input
    #
    icache = {}
    for i in range(num_files):
        hlist = hdrlist[i]
        inst = hlist.instrument.name.lower()
        tel = hlist.telescope.lower()

        if inst not in icache:
            # initialize the appropriate class for this data's instrument
            inst_class = KeywordRules(inst, telescope=tel,
                                      rules_file=rules_file)
            if verbose:
                print("Found RULEFILE for {}/{} of: {}".format(tel, inst,
                      inst_class.rules_file))
            # Interpret rules for this class based on image that
            # initialized this instrument's rules
            inst_class.interpret_rules(hlist)

            # Now add this interpreted class to the cache
            icache[inst] = inst_class

    # Create final merged set of rules
    final_rules = None
    for inst in icache:
        if final_rules is None:
            final_rules = icache[inst]
        else:
            final_rules.merge(icache[inst])

    # Apply rules to each set of input headers
    new_headers = []

    # apply rules to PRIMARY headers separately, since there is only
    # 1 PRIMARY header per image, yet many extension headers
    newhdr, newtab = final_rules.apply(hdrlist)
    final_rules.add_rules_kws(newhdr)
    new_headers.append(newhdr)

    if len(newtab) > 0:
        # Now merge the results for all the tables into single table extension
        new_table = fits.BinTableHDU.from_columns(newtab)
        new_table.header['EXTNAME'] = 'HDRTAB'
    else:
        new_table = None
    return new_headers, new_table


def merge_tables_by_cols(tables):
    """
    Merge all input tables provided as a list of np.ndarray objects into a
    single np.ndarray object
    """
    # build new list of dtypes for all cols from all tables
    # However, skip duplicate columns from subsequent tables
    # Only the values from the last table will be kept for duplicate columns
    new_dtypes = []
    new_cnames = []
    num_rows = len(tables[0])
    print('#'*40)
    print('')
    print('New merge...')
    print('')
    print('#'*40)
    for t in tables:
        for cname in t.dtype.names:
            if cname not in new_cnames:
                new_cnames.append(cname)
                new_dtypes.append((cname, t[cname].dtype.str))
        print('#'*40)
        print(t.dtype)
        print('#'*40)
    # create an empty table with the combined dtypes from all input tables
    new_table = np.zeros((num_rows,), dtype=new_dtypes)
    print(new_dtypes)
    # copy data column-by-column from each input to new output table
    for t in tables:
        for cname in t.dtype.names:
            print('CNAME: ', cname, ' with dtype: ', t[cname].dtype,
                  new_table[cname].dtype)
            new_table[cname] = t[cname]

    return new_table


def cat_headers(hdr1, hdr2):
    """
    Create new `astropy.io.fits.Header` object from
    concatenating 2 input Headers
    """
    nhdr = hdr1.copy()
    for c in hdr2.cards:
        nhdr.append(c)

    return fits.Header(nhdr)


def extract_filenames_from_product(product):
    """
    Returns the list of filenames with extensions of input observations that
    were used to generate the product.
    """
    asn_table = product.meta.asn.table_name
    asn = associations.load_asn(asn_table)
    prod = asn['products'][0]
    fnames = [m['expname'] for m in prod['members']
              if m['exptype'].upper() == 'SCIENCE']
    return fnames


def build_tab_schema(new_table):
    """
    Return new schema definition that describes the input table.

    """
    hdrtab = OrderedDict()
    hdrtab['title'] = 'Combined header table'
    hdrtab['fits_hdu'] = 'HDRTAB'
    datatype = []
    for col in new_table.columns:
        cname = col.name
        ctype = convert_dtype(str(col.dtype))
        c = OrderedDict()
        c['name'] = cname
        c['datatype'] = ctype
        datatype.append(c)
    hdrtab['datatype'] = datatype

    return hdrtab


def convert_dtype(value):
    """
    Convert numarray column dtype into YAML-compatible format description
    """
    if 'S' in value:
        # working with a string description
        str_len = int(value[value.find('S') + 1:])
        new_dtype = ['ascii', str_len]
    else:
        new_dtype = str(value)

    return new_dtype
