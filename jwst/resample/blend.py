from collections import OrderedDict
import six

from fitsblender import blendheaders

from .. import datamodels
from ..datamodels import schema
from ..datamodels import fits_support ## NEW

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def blendmetadata(input_list, output_model):
    """ Primary interface for JWST pipeline use of fitsblender.blendheaders

    This function will update the output_model datamodel with the blended metadata
    from the list of FITS objects generated from the input_list of filenames.
     
    """
    new_hdrs,new_table = blendheaders.get_blended_headers(input_list)
    
    # Now merge the keyword values from new_hdrs into the metatdata for the
    # output datamodel
    # 
    # start by building dict which maps all FITS keywords in schema to their 
    # attribute in the schema
    fits_dict = schema.build_fits_dict(output_model.schema)
    # Now assign values from new_hdrs to output_model.meta using fits_dict map
    for hdr in new_hdrs:
        for kw in hdr:
            if kw in fits_dict:
                output_model[fits_dict[kw]] = hdr[kw]

    # Now, append HDRTAB as new element in datamodel
    new_schema = build_tab_schema(new_table)
    output_model.add_schema_entry('hdrtab',new_schema)
    output_model.hdrtab = fits_support.from_fits_hdu(new_table, new_schema) ## CHANGED
        
def build_tab_schema(new_table):
    """
    Return new schema definition that describes the input table.
    
    """
    hdrtab = OrderedDict()
    hdrtab['title']='Combined header table'
    hdrtab['fits_hdu'] = 'HDRTAB'
    datatype = []
    for col in new_table.columns:
        cname = col.name
        ctype = convert_dtype(str(col.dtype))
        c = OrderedDict()
        c['name'] = cname
        c['datatype'] = ctype
        datatype.append(c)
    hdrtab['datatype']=datatype
    
    return hdrtab
    
def convert_dtype(value):
    """
    Convert numarray column dtype into YAML-compatible format description
    """
    if 'S' in value:
        # working with a string description
        str_len = int(value[value.find('S')+1:])
        new_dtype = [u'ascii', str_len] ## CHANGED
    else:
        new_dtype = unicode(str(value))

    return new_dtype
