#!/usr/bin/env python
from ..stpipe import Pipeline
from .. import datamodels

import os
import json

# calniriss ami step imports
from ..ami import ami_analyze_step
from ..ami import ami_average_step
from ..ami import ami_normalize_step


__version__ = "1.0"

# Define logging
import logging
log = logging.getLogger()
log.setLevel(logging.DEBUG)

class NirissAmiPipeline(Pipeline):
    """

    NirissAmiPipeline: Apply all level-3 calibration steps to one or
    more level-2 NIRISS AMI exposures. Included steps are:
    ami_analyze (fringe detection)
    ami_average (average results of fringe detection)
    ami_normalize (normalize results by reference target)

    """

    spec = """
    """

    # Define alias to steps
    step_defs = {'ami_analyze': ami_analyze_step.AmiAnalyzeStep,
                 'ami_average': ami_average_step.AmiAverageStep,
                 'ami_normalize': ami_normalize_step.AmiNormalizeStep
                 }

    def process(self, input):

        log.info('Starting calniriss_ami ...')

        # Retrieve the input(s)
        input_table = AmiInput(input)

        # Process all 'LGAVG' output products that are present in the ASN
        for prod in input_table.asn['products']:
            if prod['prodtype'][:5] == 'LGAVG':

                # Loop over all the members of this product
                lg_files = []
                for member in prod['members']:

                    # Do the LG analysis for this member
                    input_file = member['expname']
                    self.log.debug(' Do LG processing for member %s', input_file)

                    result = self.ami_analyze(input_file)

                    # Construct and save the output file name for this member
                    output_file = lg_product_name(input_file)
                    lg_files.append(output_file)

                    # Save the output file and close the models for this member
                    self.log.debug(' Saving LG results to %s', output_file)
                    result.save(output_file)
                    result.close()

                # Now produce the LGAVG product from these member results
                if len(lg_files) > 1:
                    self.log.debug(' Calling average_LG ...')
                    result = self.ami_average(lg_files)

                    self.log.debug(' Saving results to %s', prod['name'])
                    result.save(prod['name'])


        # Now that all LGAVG products have been produced, do normalization of
        # the target results by reference results, if reference results exist
        for prod in input_table.asn['products']:

            # Find the LGNORM product in the asn table
            if prod['prodtype'] == 'LGNORM':

                # Find the target and reference members for this product
                target_name = None
                ref_name = None
                for member in prod['members']:
                    if member['exptype'] == 'LGAVGT':
                        target_name = member['expname']
                    if member['exptype'] == 'LGAVGR':
                        ref_name = member['expname']

                # Apply the normalization
                if target_name is not None and ref_name is not None:

                    self.log.debug(' Call normalize_LG for %s and %s', target_name, ref_name)
                    result = self.ami_normalize(target_name, ref_name)

                    # Save the result
                    self.log.debug(' Saving result to %s', prod['name'])
                    result.save(prod['name'])
                    result.close()


        # We're done
        log.info('... ending calniriss_ami')

        return


class AmiInput(object):
    """Class to handle reading the input to the AMI processing, which
       can be a single science exposure or an AMI-mode association table.
       The input and output member info is loaded into an ASN table model.
    """

    template = {"asn_rule": "",
              "targname": "",
              "asn_pool": "",
              "program": "",
              "asn_type": "ami",
              "products": [
                  {"name": "",
                   "prodtype": "",
                   "members": [
                      {"exptype": "",
                       "expname": ""}
                      ]
                   }
                ]
              }

    def __init__(self, input):
        self.input = input # keep a record of original input name for later

        if isinstance(input, datamodels.ImageModel):
            # It's a single image that's been passed in as a model
            self.interpret_image_model(input)
        elif isinstance(input, str):
            try:
                # The name of an association table
                self.asn = json.load(open(input, 'r'))
            except:
                # The name of a single image file
                self.interpret_image_model(datamodels.ImageModel(input))
        else:
            raise TypeError

    def interpret_image_model(self, model):
        """ Interpret image model as single member association data product.
        """
        self.input_models = []
        self.filename = model.meta.filename

        # An in-memory ImageModel for a single exposure was provided as input
        self.asn = self.template
        self.asn['targname'] = model.meta.target.catalog_name
        self.asn['program'] = model.meta.observation.program_number
        self.asn['asn_rule'] = 'Singleton_{0}'.format(model.meta.instrument.name)

        self.asn['products'][0]['name'] = lg_product_name(self.filename)
        self.asn['products'][0]['prodtype'] = 'LGAVGT'

        self.rootname = self.filename[:self.filename.rfind('_')]
        self.asn['products'][0]['members'][0]['expname'] = self.filename
        self.input_models.append(model)

def lg_product_name(filename):
    return os.path.splitext(filename)[0] + "_{0}.fits".format("lg")
