#!/usr/bin/env python
from ..associations import Association
from ..stpipe import Pipeline
from .. import datamodels

# step imports
from ..assign_wcs import assign_wcs_step
from ..background import background_step
from ..imprint import imprint_step
#from jwst.msaflagging import msa_flag_step
from ..extract_2d import extract_2d_step
from ..flatfield import flat_field_step
from ..srctype import srctype_step
from ..straylight import straylight_step
from ..fringe import fringe_step
from ..photom import photom_step

__version__ = "2.0"

# Define logging
import logging
log = logging.getLogger()
log.setLevel(logging.DEBUG)

class Spec2Pipeline(Pipeline):
    """
    Spec2Pipeline: Processes JWST spectroscopic exposures from Level 2a to 2b.
    Accepts a single exposure or an association as input.

    Included steps are:
    assign_wcs, background subtraction, NIRSpec MSA imprint subtraction,
    NIRSpec MSA bad shutter flagging, 2-D subwindow extraction, flat field,
    source type decision, straylight, fringe, and photom.
    """

    spec = """
        save_bsub = boolean(default=False)
    """

    # Define aliases to steps
    step_defs = {'assign_wcs': assign_wcs_step.AssignWcsStep,
                 'bkg_subtract': background_step.BackgroundStep,
                 'imprint_subtract': imprint_step.ImprintStep,
                 #'msa_flagging' : msa_flag_step.MsaFlagStep,
                 'extract_2d': extract_2d_step.Extract2dStep,
                 'flat_field': flat_field_step.FlatFieldStep,
                 'srctype': srctype_step.SourceTypeStep,
                 'straylight': straylight_step.StraylightStep,
                 'fringe': fringe_step.FringeStep,
                 'photom': photom_step.PhotomStep
                }

    # The main process method
    def process(self, input):

        log.info('Starting calwebb_spec2 ...')

        # Retrieve the input(s)
        input_table = Lvl2Input(input)

        # Loop over all the members
        for member in input_table.asn['members']:

            input_file = member['expname']
            self.log.debug(' Working on %s ...', input_file)
            input = datamodels.open(input_file)

            # Apply WCS info
            input = self.assign_wcs(input)

            # Do background processing
            if len(member['bkgexps']) > 0:

                # Assemble the list of background exposures to use
                bkg_list = []
                for bkg in member['bkgexps']:
                    bkg_list.append(bkg['expname'])

                # Call the background subtraction step
                input = self.bkg_subtract(input, bkg_list)

                # Save the background-subtracted product, if requested
                if self.save_bsub:
                    if isinstance(input, datamodels.CubeModel):
                        self.save_model(input, "bsubints")
                    else:
                        self.save_model(input, "bsub")

            # Apply NIRSpec MSA imprint subtraction
            if input.meta.exposure.type in ['NRS_MSASPEC', 'NRS_IFU']:
                if len(member['imprint']) > 0:
                    imprint_filename = member['imprint'][0]['expname']
                    input = self.imprint_subtract(input, imprint_filename)

            # Apply NIRSpec MSA bad shutter flagging
            # Stubbed out as placeholder until step module is created
            #if input.meta.exposure.type in ['NRS_MSASPEC', 'NRS_IFU']:
            #    input = self.msa_flagging(input)

            # Extract 2D sub-windows for NIRSpec slit and MSA
            if input.meta.exposure.type in ['NRS_FIXEDSLIT', 'NRS_MSASPEC']:
                input = self.extract_2d(input)

            # Apply flat-field correction
            input = self.flat_field(input)

            # Apply the source type decision step
            input = self.srctype(input)

            # Apply the straylight correction for MIRI MRS
            if input.meta.exposure.type == 'MIR_MRS':
                input = self.straylight(input)

            # Apply the fringe correction for MIRI MRS
            if input.meta.exposure.type == 'MIR_MRS':
                input = self.fringe(input)

            # Apply flux calibration
            input = self.photom(input)

            # Save the calibrated exposure
            if input_table.poolname:
                input.meta.asn.pool_name = input_table.poolname
                input.meta.asn.table_name = input_table.filename
            else:
                input.meta.asn.pool_name = ' '
                input.meta.asn.table_name = ' '

            if isinstance(input, datamodels.CubeModel):
                self.save_model(input, 'calints')
            else:
                self.save_model(input, "cal")
            log.info('Save calibrated product to %s' % input.meta.filename)

            input.close()

        # We're done
        log.info('... ending calwebb_spec2')

        return


import json

class Lvl2Input(object):

    """
    Class to handle reading the input to the processing, which
    can be a single science exposure or an association table.
    The input and output member info is loaded into an ASN table model.
    """

    template = {"asn_rule": "",
              "target": "",
              "asn_pool": "",
              "asn_type": "",
              "members": [
                  {"exptype": "",
                   "expname": "",
                   "bkgexps": [],
                   "imprint": []}
                ]
              }

    def __init__(self, input):

        if isinstance(input, str):
            self.filename = input
            try:
                # The name of an association table
                with open(input, 'r') as input_fh:
                    self.asn = Association.load(input_fh)
            except:
                # The name of a single image file
                self.interpret_image_model(datamodels.open(input))
            self.poolname = self.asn['asn_pool']
        else:
            raise TypeError

    def interpret_image_model(self, model):
        """ Interpret image model as a single member association data product.
        """

        # A single exposure was provided as input
        self.asn = self.template
        self.asn['target'] = model.meta.target.catalog_name
        self.asn['asn_rule'] = 'singleton'
        self.asn['asn_type'] = 'singleton'
        self.asn['asn_pool'] = ''

        self.rootname = self.filename[:self.filename.rfind('_')]
        self.asn['members'][0]['expname'] = self.filename
