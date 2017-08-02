#!/usr/bin/env python
from collections import defaultdict
import logging

from .. import datamodels
from ..associations.load_as_asn import LoadAsLevel2Asn
from ..stpipe import Pipeline

# calwebb IMAGE2 step imports
from ..assign_wcs import assign_wcs_step
from ..flatfield import flat_field_step
from ..photom import photom_step


__version__ = "3.0"


class Image2Pipeline(Pipeline):
    """
    Image2Pipeline: Processes JWST imaging-mode slope data from Level-2a to
    Level-2b.

    Included steps are:
    assign_wcs, flat_field, and photom.
    """

    # Define alias to steps
    step_defs = {
        'assign_wcs': assign_wcs_step.AssignWcsStep,
        'flat_field': flat_field_step.FlatFieldStep,
        'photom': photom_step.PhotomStep,
        }

    def process(self, input):

        self.log.info('Starting calwebb_image2 ...')

        # Retrieve the input(s)
        asn = LoadAsLevel2Asn.load(input)

        # If `output_file` is specified and there is only
        # one product, go ahead an apply
        if self.output_file:
            if len(asn['products']) > 1:
                self.log.warn(
                    '"output_file" specified, but more than one product'
                    'will be created. Ignoring.'
                    '\nConsider using "output_dir" instead'
                )
            else:
                asn['products'][0]['name'] = self.output_file

        # Each exposure is a product in the association.
        # Process each exposure.
        for product in asn['products']:
            self.log.info('Processing product {}'.format(product['name']))
            self.output_basename = product['name']
            result = self.process_exposure_product(
                product,
                asn['asn_pool'],
                asn.filename
            )

            # Save result
            suffix = 'cal'
            self.save_model(result, suffix)

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
        input.meta.asn.pool_name = pool_name
        input.meta.asn.table_name = asn_file

        # work on slope images
        input = self.assign_wcs(input)
        input = self.flat_field(input)
        input = self.photom(input)

        # That's all folks
        self.log.info(
            'Finished processing product {}'.format(exp_product['name'])
        )
        return input
