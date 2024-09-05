"""blendmeta - Merge metadata from multiple models.

    This module will create a new metadata instance and table from a
    list of input datamodels or filenames.
"""
from collections import OrderedDict

from astropy.io import fits

from stdatamodels import fits_support
from stdatamodels import schema as dm_schema

from jwst.datamodels import ModelContainer

from .. import associations

from .blendrules import KeywordRules


__doctest_skip__ = ['blendmodels']

# Primary functional interface for the code


def blendmodels(product, inputs, ignore=None):
    """
    Run main interface for blending metadata from multiple models.

    Blend models that went into creating the original drzfile into a
    new metadata instance with a table that contains attribute values from
    all input datamodels.

    The product will be updated 'in-place' with the new metadata attributes
    and FITS BinTableHDU table.  The blended FITS table, with extname=HDRTAB,
    has 1 column for each metadata attribute recorded from the input models,
    one row for each input model, and column names are the FITS keywords for
    that metadata attribute.  For example, values from `meta.observation.time`
    would be stored in the `TIME-OBS` column.

    Rules for what function to use to determine the blended output attribute
    value and what metadata attributes should be used as columns in the
    blended FITS table are defined in the datamodel schema.


    NOTE
    ====
    Custom rules for a metadata value should be computed by the calling routine
    and used to update the metadata in the output model AFTER
    calling this function.

    Parameters
    ----------
    product : str
        Name of combined product with metadata that needs updating. This can
        be specified as a single filename.

    inputs : list
        This can be either a list of filenames or a list of DataModels objects.
        If provided, the filenames provided in this list will be used to get
        the metadata which will be blended into the final output metadata.

    ignore : list of str, None, optional
        A list of string the meta attribute names which, if provided,
        will show which attributes should not be blended.

    Example
    -------
    This example shows how to blend the metadata from a set of DataModels
    already read in memory for the product created by the `resample` step.
    This example relies on the Association file used as the input to the
    `resample` step to specify all the inputs for blending using the
    following syntax::

    >>> from stdatamodels.jwst import datamodels
    >>> asnfile = "jw99999-a3001_20170327t121212_coron3_001_asn.json"
    >>> data = datamodels.open(asnfile)
    >>> input_models = [data[3], data[4]]  # we know the last datasets are SCIENCE
    >>> blendmodels(data.meta.asn_table.products[0].name, input_models)
    """
    newmeta, newtab = get_blended_metadata(inputs)

    # open product to update metadata from input models used to create product
    output_model = product

    '''
    NOTE 17-Jan-2017:
     Effort needs to be made to insure that new blended values conform
     to the definitions of the attribute as provided by the schema.
     This will address #1650 for the jwst package.
       https://github.com/STScI-JWST/jwst/issues/1650
    '''

    # Start by identifying elements of the model which need to be ignored
    ignore_list = ['meta.wcs']  # Necessary since meta.wcs is not in schema
    if ignore:
        ignore_list.extend(ignore)

    # Now assign values from new_hdrs to output_model.meta
    flat_new_metadata = newmeta.to_flat_dict()

    for attr in flat_new_metadata:
        attr_use = not [attr.startswith(i) for i in ignore_list].count(True)
        if attr.startswith('meta') and attr_use:
            try:
                output_model[attr] = newmeta[attr]
            except KeyError:
                # Ignore keys that are in the asdf tree but not in the schema
                pass

    # Now, append HDRTAB as new element in datamodel
    newtab_schema = build_tab_schema(newtab)
    if hasattr(output_model, 'hdrtab'):
        del output_model.hdrtab
    output_model.add_schema_entry('hdrtab', newtab_schema)
    output_model.hdrtab = fits_support.from_fits_hdu(newtab, newtab_schema)

    # Clean up for the next run
    del newmeta, newtab


def get_blended_metadata(input_models):
    """
    Return a blended metadata instance and table based on the input datamodels.
    This will serve as the primary interface for blending datamodels.

    Parameters
    ----------
    input_models : list
        Either a single list of filenames from which to extract the metadata to
        be blended, or a list of `datamodels.JwstDataModel` objects to be blended.
        The input models are assumed to have the blending rules defined as
        an integral part of the schema definition for the model.

    Returns
    -------
    metadata : list
        A list of blended metadata instances, one for each i

    new_table : object
        Single fits.TableHDU object that contains the combined results from
        all input headers(extension). Each row will correspond to an image,
        and each column corresponds to a single keyword listed in the rules.

    """
    if not isinstance(input_models, list) and \
       not isinstance(input_models, ModelContainer):
        input_models = [input_models]

    num_files = len(input_models)

    # Determine what blending rules need to be merged to create the final
    # blended headers. There will be a separate set of rules for each
    # instrument, and all rules get merged into a composite set of rules that
    # get applied to all input headers regardless of instrument.
    #
    # Instrument identification will be extracted from the INSTRUME keyword
    # from the PRIMARY header of each input
    #
    icache = {}
    model_type = None
    for i in range(num_files):
        model = input_models[i]
        if model_type is None:
            model_type = type(model)
        elif type(model) != model_type:
            raise ValueError("blending only works for a single model type")
        inst = model.meta.instrument.name.lower()

        if inst not in icache:
            # initialize the appropriate class for this data's instrument
            inst_class = KeywordRules(model)
            # Interpret rules for this class based on image that
            # initialized this instrument's rules
            inst_class.interpret_rules(model)

            # Now add this interpreted class to the cache
            icache[inst] = inst_class

    # Create final merged set of rules
    final_rules = None
    for inst in icache:
        if final_rules is None:
            final_rules = icache[inst]
        else:
            final_rules.merge(icache[inst])

    # apply rules to datamodels metadata
    new_meta, newtab = final_rules.apply(input_models)

    if len(newtab) > 0:
        # Now merge the results for all the tables into single table extension
        new_table = fits.BinTableHDU.from_columns(newtab)
        new_table.header['EXTNAME'] = 'HDRTAB'
    else:
        new_table = None

    return new_meta, new_table


def build_tab_schema(new_table):
    """Return new schema definition that describes the input table."""
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
    """Convert numarray column dtype into YAML-compatible format description"""
    if 'S' in value:
        # working with a string description
        str_len = int(value[value.find('S') + 1:])
        new_dtype = ['ascii', str_len]
    else:
        new_dtype = str(value)

    return new_dtype
