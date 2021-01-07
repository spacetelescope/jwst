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


__all__ = ['Image2Pipeline']


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
    image_exptypes = ['MIR_IMAGE', 'NRC_IMAGE', 'NIS_IMAGE', 'FGS_IMAGE']

    def process(self, input):

        self.log.info('Starting calwebb_image2 ...')

        # Retrieve the input(s)
        asn = LoadAsLevel2Asn.load(input, basename=self.output_file)

        # Each exposure is a product in the association.
        # Process each exposure.
        results = []
        for product in asn['products']:
            self.log.info('Processing product {}'.format(product['name']))
            if self.save_results:
                self.output_file = product['name']
            try:
                getattr(asn, 'filename')
            except AttributeError:
                asn.filename = "singleton"

            result = self.process_exposure_product(
                product,
                asn['asn_pool'],
                op.basename(asn.filename)
            )

            # Save result
            suffix = 'cal'
            if isinstance(result, datamodels.CubeModel):
                suffix = 'calints'
            result.meta.filename = self.make_output_path(suffix=suffix)
            result.meta.filetype = 'calibrated'
            results.append(result)

        self.log.info('... ending calwebb_image2')

        self.output_use_model = True
        self.suffix = False
        return results

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

        pool_name: str
            The pool file name. Used for recording purposes only.

        asn_file: str
            The name of the association file.
            Used for recording purposes only.
        """
        # Find all the member types in the product
        members_by_type = defaultdict(list)
        for member in exp_product['members']:
            members_by_type[member['exptype'].lower()].append(member['expname'])

        # Get the science member. Technically there should only be
        # one. We'll just get the first one found.
        science = members_by_type['science']
        if len(science) != 1:
            self.log.warning(
                'Wrong number of science files found in {}'.format(
                    exp_product['name']
                )

            )
            self.log.warning('    Using only first one.')
        science = science[0]

        self.log.info('Working on input %s ...', science)
        if isinstance(science, datamodels.DataModel):
            input = science
        else:
            input = datamodels.open(science)

        # Record ASN pool and table names in output
        input.meta.asn.pool_name = pool_name
        input.meta.asn.table_name = asn_file

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
        # regular 2D science image types
        if input.meta.exposure.type.upper() in self.image_exptypes and \
        len(input.data.shape) == 2:
            self.resample.save_results = self.save_results
            self.resample.suffix = 'i2d'
            self.resample(input)

        # That's all folks
        self.log.info(
            'Finished processing product {}'.format(exp_product['name'])
        )
        return input
