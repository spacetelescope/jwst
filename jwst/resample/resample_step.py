import logging

import numpy as np

from ..stpipe import Step
from .. import datamodels
from . import resample
from ..assign_wcs.util import update_s_region


__all__ = ["ResampleStep"]


class ResampleStep(Step):
    """
    Resample input data onto a regular grid using the drizzle algorithm.

    Parameters
    -----------
    input : DataModel or Association
        Single filename for either a single image or an association table.
    """

    spec = """
        pixfrac = float(default=None)
        kernel = string(default=None)
        fillval = string(default=None)
        weight_type = option('exptime', default=None)
        good_bits = integer(min=0, default=4)
        single = boolean(default=False)
        blendheaders = boolean(default=True)
    """

    reference_file_types = ['drizpars']

    def process(self, input):

        input = datamodels.open(input)

        # If single input, wrap in a ModelContainer
        if not isinstance(input, datamodels.ModelContainer):
            input_models = datamodels.ModelContainer([input])
            input_models.meta.resample.output = input.meta.filename
            self.blendheaders = False
        else:
            input_models = input

        for reftype in self.reference_file_types:
            ref_filename = self.get_reference_file(input_models[0], reftype)
        if ref_filename is not None:
            self.log.info('Drizpars reference file: {}'.format(ref_filename))
            kwargs = self.get_drizpars(ref_filename, input_models)

        # Call the resampling routine
        resamp = resample.ResampleData(input_models, **kwargs)
        resamp.do_drizzle()

        for model in resamp.output_models:
            model.meta.cal_step.resample = "COMPLETE"
            update_s_region(model)
            model.meta.asn.pool_name = input_models.meta.pool_name
            model.meta.asn.table_name = input_models.meta.table_name


        if len(resamp.output_models) == 1:
            result = resamp.output_models[0]
        else:
            result = resamp.output_models

        return result


    def get_drizpars(self, ref_filename, input_models):
        """
        Extract drizzle parameters from reference file.
        """
        drizpars_table = datamodels.DrizParsModel(ref_filename).data

        num_groups = len(input_models.group_names)
        filtname = input_models[0].meta.instrument.filter
        row = None
        filter_match = False
        # look for row that applies to this set of input data models
        for n, filt, num in zip(
            range(0, len(drizpars_table)),
            drizpars_table['filter'],
            drizpars_table['numimages']
            ):
            # only remember this row if no exact match has already been made for
            # the filter. This allows the wild-card row to be anywhere in the
            # table; since it may be placed at beginning or end of table.

            if str(filt) == "ANY" and not filter_match and num_groups >= num:
                row = n
            # always go for an exact match if present, though...
            if filtname == filt and num_groups >= num:
                row = n
                filter_match = True

        # With presence of wild-card rows, code should never trigger this logic
        if row is None:
            self.log.error("No row found in %s matching input data.", ref_filename)
            raise ValueError

        # Define the keys to then pull from drizpars reffile table.  Not the
        # step param 'weight_type' is 'wht_type' in the FITS binary table
        drizpars = dict(
            pixfrac=self.pixfrac,
            kernel=self.kernel,
            fillval=self.fillval,
            wht_type=self.weight_type
            )

        # For parameters that are set in drizpars table but not set by the
        # user, use these.  Otherwise, use values set by user.
        reffile_drizpars = {k:v for k,v in drizpars.items() if v is None}
        user_drizpars = {k:v for k,v in drizpars.items() if v is not None}

        # read in values from that row for each parameter
        for k in reffile_drizpars:
            if k in drizpars_table.names:
                reffile_drizpars[k] = drizpars_table[k][row]
        # Convert the 'wht_type' key to a 'weight_type' key
        reffile_drizpars['weight_type'] = reffile_drizpars.pop('wht_type')

        # Convert the strings in the FITS binary table from np.bytes_ to str
        for k,v in reffile_drizpars.items():
            if isinstance(v, np.bytes_):
                reffile_drizpars[k] = v.decode('UTF-8')

        kwargs = dict(
            good_bits=self.good_bits,
            single=self.single,
            blendheaders=self.blendheaders
            )
        kwargs.update({**reffile_drizpars, **user_drizpars})

        for k,v in kwargs.items():
            self.log.debug('   {}={}'.format(k, v))

        return kwargs

