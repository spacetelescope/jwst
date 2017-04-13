#!/usr/bin/env python
from ..associations import load_asn
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
from ..pathloss import pathloss_step
from ..photom import photom_step
from ..cube_build import cube_build_step
from ..extract_1d import extract_1d_step
from ..resample import resample_spec_step

__version__ = "3.0"

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
    source type decision, straylight, fringe, pathloss, photom, resample_spec,
    cube_build, and extract_1d.
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
                 'pathloss': pathloss_step.PathLossStep,
                 'photom': photom_step.PhotomStep,
                 'resample_spec': resample_spec_step.ResampleSpecStep,
                 'cube_build': cube_build_step.CubeBuildStep,
                 'extract_1d': extract_1d_step.Extract1dStep
                }

    # The main process method
    def process(self, input):

        log.info('Starting calwebb_spec2 ...')

        # Retrieve the input(s)
        input_table = Lvl2Input(input)

        # Loop over all the members
        for member in input_table.asn['members']:

            input_file = member['expname']

            # Skip processing of non-science members
            if member['exptype'].upper() != 'SCIENCE':
                self.log.info('Skipping non-science input %s', input_file)
                continue

            self.log.info('Working on input %s ...', input_file)
            input = datamodels.open(input_file)
            exp_type = input.meta.exposure.type

            # Apply WCS info
            input = self.assign_wcs(input)

            # If assign_wcs was skipped, abort the rest of processing,
            # because so many downstream steps depend on the WCS
            if input.meta.cal_step.assign_wcs == 'SKIPPED':
                log.error('Assign_wcs processing was skipped')
                log.error('Aborting remaining processing for this exposure')
                log.error('No output product will be created')
                continue

            # Do background processing, if necessary
            if 'bkgexps' in member and len(member['bkgexps']) > 0:

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
            if exp_type in ['NRS_MSASPEC', 'NRS_IFU']:
                if 'imprint' in member and len(member['imprint']) > 0:
                    imprint_filename = member['imprint'][0]['expname']
                    input = self.imprint_subtract(input, imprint_filename)

            # Apply NIRSpec MSA bad shutter flagging
            # Stubbed out as placeholder until step module is created
            #if exp_type in ['NRS_MSASPEC', 'NRS_IFU']:
            #    input = self.msa_flagging(input)

            # Extract 2D sub-windows for NIRSpec slit and MSA
            if exp_type in ['NRS_FIXEDSLIT', 'NRS_BRIGHTOBJ', 'NRS_MSASPEC']:
                input = self.extract_2d(input)

            # Apply flat-field correction
            input = self.flat_field(input)

            # Apply the source type decision step
            input = self.srctype(input)

            # Apply the straylight correction for MIRI MRS
            if exp_type == 'MIR_MRS':
                input = self.straylight(input)

            # Apply the fringe correction for MIRI MRS
            if exp_type == 'MIR_MRS':
                input = self.fringe(input)

            # Apply pathloss correction to NIRSpec exposures
            if exp_type in ['NRS_FIXEDSLIT', 'NRS_BRIGHTOBJ', 'NRS_MSASPEC',
                            'NRS_IFU']:
                input = self.pathloss(input)

            # Apply flux calibration
            input = self.photom(input)

            # Record ASN pool and table names in output
            if input_table.poolname:
                input.meta.asn.pool_name = input_table.poolname
                input.meta.asn.table_name = input_table.filename
            else:
                input.meta.asn.pool_name = ' '
                input.meta.asn.table_name = ' '

            # Save the calibrated exposure
            if isinstance(input, datamodels.CubeModel):
                self.save_model(input, 'calints')
            else:
                self.save_model(input, "cal")
            log.info('Saved calibrated product to %s' % input.meta.filename)

            # Produce a resampled product, either via resample_spec for
            # "regular" spectra or cube_build for IFU data. No resampled
            # product is produced for time-series modes.
            if input.meta.exposure.type in ['NRS_FIXEDSLIT', 'NRS_BRIGHTOBJ',
                'NRS_MSASPEC', 'NIS_WFSS', 'NRC_GRISM']:

                # Call the resample_spec step
                resamp = self.resample_spec(input)

                # Save the resampled product
                if self.resample_spec.skip != True:
                    self.save_model(resamp, 's2d')
                    log.info('Saved resampled product to %s' % resamp.meta.filename)

                # Pass the resampled data to 1D extraction
                x1d_input = resamp.copy()
                resamp.close()

            elif exp_type in ['MIR_MRS', 'NRS_IFU']:

                # Call the cube_build step for IFU data
                cube = self.cube_build(input)

                # Save the cube product
                if self.cube_build.skip != True:
                    self.save_model(cube, 's3d')
                    log.info('Saved IFU cube product to %s' % cube.meta.filename)

                # Pass the cube along for input to 1D extraction
                x1d_input = cube.copy()
                cube.close()

            else:
                # Pass the unresampled cal product to 1D extraction
                x1d_input = input

            # Extract a 1D spectrum from the 2D/3D data
            x1d_output = self.extract_1d(x1d_input)
            x1d_input.close()

            # Save the extracted spectrum
            if self.extract_1d.skip != True:
                if isinstance(input, datamodels.CubeModel):
                    self.save_model(x1d_output, 'x1dints')
                else:
                    self.save_model(x1d_output, 'x1d')
                log.info('Saved extracted spectrum to %s' % x1d_output.meta.filename)

            input.close()
            x1d_output.close()

        # We're done
        log.info('Ending calwebb_spec2')

        return


class Lvl2Input(object):

    """
    Class to handle reading the input to the processing, which
    can be a single science exposure or an association table.
    The input member info is loaded into an ASN table model.
    """

    template = {"asn_pool": "",
                "members": [
                  {"expname": "",
                   "exptype": ""}
                 ]
                }

    def __init__(self, input):

        if isinstance(input, str):
            self.filename = input
            try:
                # The name of an association table
                with open(input, 'r') as input_fh:
                    self.asn = load_asn(input_fh)
            except:
                # The name of a single image file
                self.interpret_image_model(datamodels.open(input))
            self.poolname = self.asn['asn_pool']

        elif isinstance(input, datamodels.DataModel):
            # A single data model
            self.filename = input.meta.filename
            self.interpret_image_model(input)
            self.poolname = self.asn['asn_pool']

        else:
            raise TypeError

    def interpret_image_model(self, model):
        """ Interpret image model as a single member association
        """

        # A single exposure was provided as input
        self.asn = self.template
        self.asn['asn_pool'] = ''
        self.asn['members'][0]['expname'] = self.filename
        self.asn['members'][0]['exptype'] = 'science'
