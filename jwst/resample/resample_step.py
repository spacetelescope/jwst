import logging

import numpy as np

from stpipe.extern.configobj.validate import Validator
from stpipe.extern.configobj.configobj import ConfigObj

from . import resample
from ..stpipe import Step
from .. import datamodels
from ..assign_wcs import util

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["ResampleStep"]


# Force use of all DQ flagged data except for DO_NOT_USE and NON_SCIENCE
GOOD_BITS = '~DO_NOT_USE+NON_SCIENCE'


class ResampleStep(Step):
    """
    Resample input data onto a regular grid using the drizzle algorithm.

    Parameters
    -----------
    input :  ~jwst.datamodels.DataModel or ~jwst.associations.Association
        Single filename for either a single image or an association table.
    """

    spec = """
        pixfrac = float(default=1.0)
        kernel = string(default='square')
        fillval = string(default='INDEF')
        weight_type = option('ivm', 'exptime', default='ivm')
        pixel_scale_ratio = float(default=1.0) # Ratio of output to input pixel scale
        single = boolean(default=False)
        blendheaders = boolean(default=True)
        allowed_memory = float(default=None)  # Fraction of memory to use for the combined image.
    """

    reference_file_types = ['drizpars']

    def process(self, input):

        input = datamodels.open(input)

        if isinstance(input, datamodels.ModelContainer):
            input_models = input
            try:
                output = input_models.meta.asn_table.products[0].name
            except AttributeError:
                # coron data goes through this path by the time it gets to
                # resampling.
                # TODO: figure out why and make sure asn_table is carried along
                output = None
        else:
            input_models = datamodels.ModelContainer([input])
            output = input.meta.filename
            self.blendheaders = False

        # Check that input models are 2D images
        if len(input_models[0].data.shape) != 2:
            # resample can only handle 2D images, not 3D cubes, etc
            raise RuntimeError("Input {} is not a 2D image.".format(input_models[0]))

        # Get drizzle parameters reference file
        for reftype in self.reference_file_types:
            ref_filename = self.get_reference_file(input_models[0], reftype)

        if ref_filename != 'N/A':
            self.log.info('Drizpars reference file: {}'.format(ref_filename))
            kwargs = self.get_drizpars(ref_filename, input_models)
        else:
            # If there is no drizpars reffile
            self.log.info("No DIRZPARS reffile")
            kwargs = self._set_spec_defaults()

        kwargs['allowed_memory'] = self.allowed_memory

        # Call the resampling routine
        resamp = resample.ResampleData(input_models, output=output, **kwargs)
        result = resamp.do_drizzle()

        for model in result:
            model.meta.cal_step.resample = 'COMPLETE'
            util.update_s_region_imaging(model)
            model.meta.asn.pool_name = input_models.asn_pool_name
            model.meta.asn.table_name = input_models.asn_table_name
            self.update_phot_keywords(model)
            model.meta.filetype = 'resampled'

        if len(result) == 1:
            result = result[0]

        return result

    def update_phot_keywords(self, model):
        """Update pixel scale keywords"""
        if model.meta.photometry.pixelarea_steradians is not None:
            model.meta.photometry.pixelarea_steradians *= self.pixel_scale_ratio**2
        if model.meta.photometry.pixelarea_arcsecsq is not None:
            model.meta.photometry.pixelarea_arcsecsq *= self.pixel_scale_ratio**2
        model.meta.resample.pixel_scale_ratio = self.pixel_scale_ratio
        model.meta.resample.pixfrac = self.pixfrac

    def get_drizpars(self, ref_filename, input_models):
        """
        Extract drizzle parameters from reference file.

        This method extracts parameters from the drizpars reference file and
        uses those to set defaults on the following ResampleStep configuration
        parameters:

        pixfrac = float(default=None)
        kernel = string(default=None)
        fillval = string(default=None)

        Once the defaults are set from the reference file, if the user has
        used a resample.cfg file or run ResampleStep using command line args,
        then these will overwerite the defaults pulled from the reference file.
        """
        with datamodels.DrizParsModel(ref_filename) as drpt:
            drizpars_table = drpt.data

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

        # Define the keys to pull from drizpars reffile table.
        # All values should be None unless the user set them on the command
        # line or in the call to the step
        drizpars = dict(
            pixfrac=self.pixfrac,
            kernel=self.kernel,
            fillval=self.fillval,
            pscale_ratio=self.pixel_scale_ratio,
        )

        # For parameters that are set in drizpars table but not set by the
        # user, use these.  Otherwise, use values set by user.
        reffile_drizpars = {k: v for k, v in drizpars.items() if v is None}
        user_drizpars = {k: v for k, v in drizpars.items() if v is not None}

        # read in values from that row for each parameter
        for k in reffile_drizpars:
            if k in drizpars_table.names:
                reffile_drizpars[k] = drizpars_table[k][row]

        # Convert the strings in the FITS binary table from np.bytes_ to str
        for k, v in reffile_drizpars.items():
            if isinstance(v, np.bytes_):
                reffile_drizpars[k] = v.decode('UTF-8')

        all_drizpars = {**reffile_drizpars, **user_drizpars}

        kwargs = dict(
            good_bits=GOOD_BITS,
            single=self.single,
            blendheaders=self.blendheaders
        )

        kwargs.update(all_drizpars)

        for k, v in kwargs.items():
            self.log.debug('   {}={}'.format(k, v))

        return kwargs

    def _set_spec_defaults(self):
        """NIRSpec currently has no default drizpars reference file, so default
        drizzle parameters are not set properly.  This method sets them.

        Remove this class method when a drizpars reffile is delivered.
        """
        configspec = self.load_spec_file()
        config = ConfigObj(configspec=configspec)
        if config.validate(Validator()):
            kwargs = config.dict()

        # Force definition of good bits
        kwargs['good_bits'] = GOOD_BITS

        if kwargs['pixfrac'] is None:
            kwargs['pixfrac'] = 1.0
        if kwargs['kernel'] is None:
            kwargs['kernel'] = 'square'
        if kwargs['fillval'] is None:
            kwargs['fillval'] = 'INDEF'
        if kwargs['weight_type'] is None:
            kwargs['weight_type'] = 'ivm'
        kwargs['pscale_ratio'] = self.pixel_scale_ratio
        kwargs.pop('pixel_scale_ratio')

        for k, v in kwargs.items():
            if k in ['pixfrac', 'kernel', 'fillval', 'weight_type', 'pscale_ratio']:
                log.info('  setting: %s=%s', k, repr(v))

        return kwargs
