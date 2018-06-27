#!/usr/bin/env python
from collections import defaultdict
import os.path as op

from .. import datamodels
from ..associations.load_as_asn import LoadAsLevel2Asn
from ..stpipe import Pipeline

# calwebb IMAGE2 step imports
from ..background import background_step
from ..assign_wcs import assign_wcs_step
from ..flatfield import flat_field_step
from ..photom import photom_step
from ..resample import resample_step


__version__ = '0.9.3'


class Image2Pipeline(Pipeline):
    """
    Image2Pipeline: Processes JWST imaging-mode slope data from Level-2a to
    Level-2b.

    Included steps are:
    background_subtraction, assign_wcs, flat_field, photom and resample.
    """

    spec = """
        save_bsub = boolean(default=False) # Save background-subracted science
    """

    # Define alias to steps
    step_defs = {
        'bkg_subtract': background_step.BackgroundStep,
        'assign_wcs': assign_wcs_step.AssignWcsStep,
        'flat_field': flat_field_step.FlatFieldStep,
        'photom': photom_step.PhotomStep,
        'resample': resample_step.ResampleStep
        }

    # List of normal imaging exp_types
    image_exptypes = ['MIR_IMAGE', 'NRC_IMAGE', 'NIS_IMAGE']

    def process(self, input):

        self.log.info('Starting calwebb_image2 ...')

        # Retrieve the input(s)
        asn = LoadAsLevel2Asn.load(input, basename=self.output_file)

        # Each exposure is a product in the association.
        # Process each exposure.
        for product in asn['products']:
            self.log.info('Processing product {}'.format(product['name']))
            self.output_file = product['name']
            result = self.process_exposure_product(
                product,
                asn['asn_pool'],
                asn.filename
            )

            # Save result
            suffix = 'cal'
            if isinstance(result, datamodels.CubeModel):
                suffix = 'calints'
            self.save_model(result, suffix)

            self.closeout(to_close=[result])

        self.log.info('... ending calwebb_image2')

    # Process each exposure
    def process_exposure_product(
            self,
            exp_product,
            pool_name=' ',
            asn_file=' '
    ):
        """Process an exposure found in the association product

        Parameters
        ---------
        exp_product: dict
            A Level2b association product.
        """
        # Find all the member types in the product
        members_by_type = defaultdict(list)
        for member in exp_product['members']:
            members_by_type[member['exptype'].lower()].append(member['expname'])

        # Get the science member. Technically there should only be
        # one. We'll just get the first one found.
        science = members_by_type['science']
        if len(science) != 1:
            self.log.warn(
                'Wrong number of science files found in {}'.format(
                    exp_product['name']
                )

            )
            self.log.warn('    Using only first one.')
        science = science[0]

        self.log.info('Working on input %s ...', science)
        if isinstance(science, datamodels.DataModel):
            input = science
        else:
            input = datamodels.open(science)

        # Record ASN pool and table names in output
        input.meta.asn.pool_name = op.basename(pool_name)
        input.meta.asn.table_name = op.basename(asn_file)

        # Do background processing, if necessary
        if len(members_by_type['background']) > 0:

            # Setup for saving
            self.bkg_subtract.suffix = 'bsub'
            if isinstance(input, datamodels.CubeModel):
                self.bkg_subtract.suffix = 'bsubints'

            # Backwards compatibility
            if self.save_bsub:
                self.bkg_subtract.save_results = True

            # Call the background subtraction step
            input = self.bkg_subtract(input, members_by_type['background'])

        # work on slope images
        input = self.assign_wcs(input)
        input = self.flat_field(input)
        input = self.photom(input)

        # Resample individual exposures, but only if it's one of the
        # regular science image types.
        if input.meta.exposure.type.upper() in self.image_exptypes:
            result = self.resample(input)
            if result:
                # write out resampled exposure
                self.save_model(result, suffix='i2d')
                result.close()

        # That's all folks
        self.log.info(
            'Finished processing product {}'.format(exp_product['name'])
        )
        return input
