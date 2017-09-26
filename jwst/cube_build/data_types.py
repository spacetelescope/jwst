# Routines used for building cubes
from __future__ import absolute_import, print_function

import sys
import time
import numpy as np
import math
import json
import os

#from astropy.io import fits
from ..associations import load_asn
from .. import datamodels


import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

#********************************************************************************
class DataTypes(object):
#********************************************************************************

    """
    Class to handle reading the input to the processing, which
    can be a single science exposure or an IFU cube association table.
    The input and output member info is loaded into an ASN table model.
    """

    template = {"asn_rule": "",
              "target": "",
              "asn_pool": "",
              "asn_type": "",
              "products": [
                  {"name": "",
                   "members": [
                      {"exptype": "",
                       "expname": ""}
                      ]
                  }
                ]
              }

    def __init__(self, input,single,output_file):

        self.input_models = []
        self.filenames = []
        self.output_name = None
        self.data_type = None # singleton, multi
        self.input_type = None # Model, File, ASN, Container

        # IF a single model or a single file  is passed in then
        # self.filename & self.input_model hold the values for this singe dataset
        self.InputType  = ''
        if isinstance(input, datamodels.ImageModel):
#            print('this is a single file passed as a Model')
            # It's a single image that's been passed in as a model
            # input is a model
            self.filenames.append(input.meta.filename)
            self.input_models.append(input)
            self.input_type = 'Model'
            self.data_type = 'singleton'
            self.output_name = self.build_product_name(self.filenames[0])
#            print('the name of the output file', self.output_name) 
        elif isinstance(input,datamodels.ModelContainer):
#            print('this is a model container type')
            self.input_type='Container'
            self.data_type = 'multi'
            self.output_name  = 'Temp'
            if not single:  # find the name of the output file from the association 
                with datamodels.ModelContainer(input) as input_model:
                    self.output_name =input_model.meta.asn_table.products[0].name

            for i in range(len(input)):
                model = input[i]
                self.input_models.append(model)
                self.filenames.append(model.meta.filename)
#            print('number of models',len(self.filenames))

        elif isinstance(input, str):
            try:
                # The name of an association table
                # for associations - use Association.load
                # in cube_build_io.SetFileTable - set up:
                # input_model & filename lists
                iproduct = 0 # only one product found in association table
                with open(input, 'r') as input_fh:
#                    print('read in association table')
                    asn_table = load_asn(input_fh)
                    self.input_type = 'ASN'
                    self.data_type = 'multi'
                    self.output_name =  asn_table['products'][0]['name']
                    for m in asn_table['products'][iproduct]['members']:
                        self.filenames.append(m['expname'])
                        self.input_models.append(datamodels.ImageModel(m['expname']))
            except:
                # The name of a single image file
#                print(' this is a single file  read in filename')
                self.input_type = 'File'
                self.data_type = 'singleton'
                self.filenames.append(input)
                self.input_models.append(datamodels.ImageModel(input))
                self.output_name = self.build_product_name(self.filenames[0])

        else:
            raise TypeError

# if the user has set the output name - strip out *.fits 
# later suffixes will be added to this name to designate the
# channel, subchannel or grating,filter the data is covers.

        if output_file !=None :
            basename,ext = os.path.splitext(os.path.basename(output_file))
#            print('basename',basename)
#            root, ext = os.path.splitext(output_file)
#            default = root.find('cube_build') # the user has not provided a name
            self.output_name = basename
        

    def build_product_name(self, filename):
        indx = filename.rfind('.fits')
        single_product = filename[:indx]
        return single_product


        
    

# TODO:  Routines not used below - saved just in case we need them later - if not
# remove. 

    def interpret_image_model(self, model):
        """ Interpret image model as single member association data product.
            Currently this routien is not used by cube_build - it was left
            if needed down the road
        """

        # An in-memory ImageModel for a single exposure was provided as input
        self.asn_table = self.template
        self.asn_table['target'] = model.meta.target.catalog_name
        self.asn_table['asn_rule'] = 'singleton'
        self.asn_table['asn_type'] = 'singleton'
        self.asn_table['products'][0]['name'] = self.build_product_name(self.filenames[0])
        self.rootname = self.filename[:self.filename.rfind('_')]
        self.asn_table['products'][0]['members'][0]['expname'] = self.filenames[0]

    def get_inputs(self, product=0):
        members = []
        for p in self.asn_table['products'][product]['members']:
            members.append(p['expname'])
        return members
    def get_outputs(self, product=0):
        return self.asn_table['products'][product]['name']


