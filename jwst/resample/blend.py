from collections import OrderedDict
import six

from fitsblender import blendheaders

from .. import datamodels
from ..datamodels import schema
from ..datamodels import fits_support

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def blendfitsdata(input_list, output_model):
    """
    Primary interface for JWST pipeline use of fitsblender.blendheaders

    This function will update the output_model datamodel with the blended metadata
    from the list of FITS objects generated from the input_list of filenames.
    """
    new_hdrs, new_table = blendheaders.get_blended_headers(input_list)

    # Now merge the keyword values from new_hdrs into the metatdata for the
    # output datamodel
    #
    # start by building dict which maps all FITS keywords in schema to their
    # attribute in the schema
    fits_dict = schema.build_fits_dict(output_model.schema)
    # Need to insure that output_model does not already have an instance
    # of hdrtab from previous processing, an instance that would be
    # incompatible with the new table generated now...
    if hasattr(output_model, 'hdrtab'):
        # If found, remove it to be replaced by new instance
        del output_model.hdrtab
    # Now assign values from new_hdrs to output_model.meta using fits_dict map
    for hdr in new_hdrs:
        for kw in hdr:
            if kw in fits_dict:
                output_model[fits_dict[kw]] = hdr[kw]

    # Now, append HDRTAB as new element in datamodel
    new_schema = build_tab_schema(new_table)

    output_model.add_schema_entry('hdrtab', new_schema)
    output_model.hdrtab = fits_support.from_fits_hdu(new_table, new_schema)


def blendmetadata(input_models, output_model):
    final_rules = build_meta_rules(input_models)

    # Apply rules to each set of input headers
    new_headers = []
    i = 0
    # apply rules to PRIMARY headers separately, since there is only
    # 1 PRIMARY header per image, yet many extension headers
    newphdr, newtab = final_rules.apply(phdrlist)
    final_rules.add_rules_kws(newphdr) # Adds HISTORY comments on rules used
    new_headers.append(newphdr)
    for hdrs in hdrlist[1:]:
        newhdr, newtab = final_rules.apply(hdrs)
        new_headers.append(newhdr)

    # create list of combined PRIMARY/SCI headers for use in creating
    # the new table extensions
    tabhdrs = []
    for phdr, scihdr in zip(hdrlist[0], hdrlist[1]):
        tabhdrs.append(cat_headers(phdr, scihdr))
    # Create extension table from list of all combined PRI/SCI headers
    tabhdr, newtab = final_rules.apply(tabhdrs)

    # Now merge the keyword values from new_hdrs into the metatdata for the
    # output datamodel
    #
    # start by building dict which maps all FITS keywords in schema to their
    # attribute in the schema
    fits_dict = schema.build_fits_dict(output_model.schema)

    # Need to insure that output_model does not already have an instance
    # of hdrtab from previous processing, an instance that would be
    # incompatible with the new table generated now...
    if hasattr(output_model, 'hdrtab'):
        # If found, remove it to be replaced by new instance
        del output_model.hdrtab

    # Now assign values from new_hdrs to output_model.meta using fits_dict map
    for hdr in new_hdrs:
        for kw in hdr:
            if kw in fits_dict:
                output_model[fits_dict[kw]] = hdr[kw]

    # Now, append HDRTAB as new element in datamodel
    new_schema = build_tab_schema(new_table)
    output_model.add_schema_entry('hdrtab', new_schema)
    output_model.hdrtab = fits_support.from_fits_hdu(new_table, new_schema)


def build_meta_rules(input_models, rules_file=None):
    # Determine what blending rules need to be merged to create the final
    # blended headers. There will be a separate set of rules for each
    # instrument, and all rules get merged into a composite set of rules that
    # get applied to all input headers regardless of instrument.
    #
    # Instrument identification will be extracted from the INSTRUME keyword
    # from the PRIMARY header of each input

    icache = {}
    for model in input_models:
        inst = model.meta.instrument.name.lower()
        tel = model.meta.telescope.lower()
        if inst not in icache:
            # initialize the appropriate class for this data's instrument
            inst_class = blendheaders.KeywordRules(inst, telescope=tel,
                rules_file=rules_file)
            log.debug("Found blendheaders RULEFILE for {}/{} of: {}".format(
                tel, inst, inst_class.rules_file))
            # Interpret rules for this class based on image that
            # initialized this instrument's rules
            inst_class.interpret_rules(model.meta)
            # Now add this interpreted class to the cache
            icache[inst] = inst_class

    # Create final merged set of rules
    final_rules = None
    for inst in icache:
        if final_rules is None:
            final_rules = icache[inst]
        else:
            final_rules.merge(icache[inst])

    return final_rules


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
        new_dtype = [u'ascii', str_len] ## CHANGED
    else:
        new_dtype = unicode(str(value))

    return new_dtype
