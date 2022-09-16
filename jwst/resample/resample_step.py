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

    class_alias = "resample"

    spec = """
        pixfrac = float(default=1.0) # change back to None when drizpar reference files are updated
        kernel = string(default='square') # change back to None when drizpar reference files are updated
        fillval = string(default='INDEF' ) # change back to None when drizpar reference files are updated
        weight_type = option('ivm', 'exptime', None, default='ivm')  # change back to None when drizpar ref update
        output_shape = int_list(min=2, max=2, default=None)  # [x, y] order
        crpix = float_list(min=2, max=2, default=None)
        crval = float_list(min=2, max=2, default=None)
        rotation = float(default=None)
        pixel_scale_ratio = float(default=1.0) # Ratio of input to output pixel scale
        pixel_scale = float(default=None) # Absolute pixel scale in arcsec
        single = boolean(default=False)
        blendheaders = boolean(default=True)
        allowed_memory = float(default=None)  # Fraction of memory to use for the combined image.
        in_memory = boolean(default=True)
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
            input_models.asn_pool_name = input.meta.asn.pool_name
            input_models.asn_table_name = input.meta.asn.table_name
            output = input.meta.filename
            self.blendheaders = False

        # Check that input models are 2D images
        if len(input_models[0].data.shape) != 2:
            # resample can only handle 2D images, not 3D cubes, etc
            raise RuntimeError("Input {} is not a 2D image.".format(input_models[0]))

        #  Get drizzle parameters reference file, if there is one
        if 'drizpars' in self.reference_file_types:
            ref_filename = self.get_reference_file(input_models[0], 'drizpars')
            self.log.info('Using drizpars reference file: {}'.format(ref_filename))
            kwargs = self.get_drizpars(ref_filename, input_models)
        else:  # no drizpars reference file found
            ref_filename = 'N/A'

        if ref_filename == 'N/A':
            self.log.info('No drizpars reference file found.')
            kwargs = self._set_spec_defaults()

        self.wht_type = self.weight_type

        kwargs['allowed_memory'] = self.allowed_memory

        # Issue a warning about the use of exptime weighting
        if self.wht_type == 'exptime':
            self.log.warning("Use of EXPTIME weighting will result in incorrect")
            self.log.warning("propagated errors in the resampled product")

        # Custom output WCS parameters.
        # Modify get_drizpars if any of these get into reference files:
        kwargs['ref_wcs'] = None  # TODO: add mechanism of specifying a ref WCS
        kwargs['out_shape'] = _check_list_pars(self.output_shape, 'output_shape',
                                               min_vals=[1, 1])
        kwargs['crpix'] = _check_list_pars(self.crpix, 'crpix')
        kwargs['crval'] = _check_list_pars(self.crval, 'crval')
        kwargs['rotation'] = self.rotation
        kwargs['pscale'] = self.pixel_scale
        kwargs['pscale_ratio'] = self.pixel_scale_ratio
        kwargs['in_memory'] = self.in_memory

        # Call the resampling routine
        resamp = resample.ResampleData(input_models, output=output, **kwargs)
        result = resamp.do_drizzle()

        for model in result:
            model.meta.cal_step.resample = 'COMPLETE'
            util.update_s_region_imaging(model)
            model.meta.asn.pool_name = input_models.asn_pool_name
            model.meta.asn.table_name = input_models.asn_table_name

            # if pixel_scale exists, it will override pixel_scale_ratio.
            # calculate the actual value of pixel_scale_ratio based on pixel_scale
            # because source_catalog uses this value from the header.
            if not self.pixel_scale:
                model.meta.resample.pixel_scale_ratio = self.pixel_scale_ratio
            else:
                model.meta.resample.pixel_scale_ratio = (self.pixel_scale ** 2) / model.meta.photometry.pixelarea_arcsecsq
            model.meta.resample.pixfrac = kwargs['pixfrac']
            self.update_phot_keywords(model)

        if len(result) == 1:
            result = result[0]

        input_models.close()
        return result

    def update_phot_keywords(self, model):
        """Update pixel scale keywords"""
        if model.meta.photometry.pixelarea_steradians is not None:
            model.meta.photometry.pixelarea_steradians *= model.meta.resample.pixel_scale_ratio**2
        if model.meta.photometry.pixelarea_arcsecsq is not None:
            model.meta.photometry.pixelarea_arcsecsq *= model.meta.resample.pixel_scale_ratio**2

    def get_drizpars(self, ref_filename, input_models):
        """
        Extract drizzle parameters from reference file.

        This method extracts parameters from the drizpars reference file and
        uses those to set defaults on the following ResampleStep configuration
        parameters:

        pixfrac = float(default=None)
        kernel = string(default=None)
        fillval = string(default=None)
        wht_type = option('ivm', 'exptime', None, default=None)

        Once the defaults are set from the reference file, if the user has
        used a resample.cfg file or run ResampleStep using command line args,
        then these will overwrite the defaults pulled from the reference file.
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
            wht_type=self.weight_type
            # pscale_ratio=self.pixel_scale_ratio, # I think this can be removed JEM (??)
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

        if self.pixfrac is None:
            self.pixfrac = 1.0
        if self.kernel is None:
            self.kernel = 'square'
        if self.fillval is None:
            self.fillval = 'INDEF'
        # Force definition of good bits
        kwargs['good_bits'] = GOOD_BITS
        kwargs['pixfrac'] = self.pixfrac
        kwargs['kernel'] = str(self.kernel)
        kwargs['fillval'] = str(self.fillval)
        #  self.weight_type has a default value of None
        # The other instruments read this parameter from a reference file
        if self.wht_type is None:
            self.wht_type = 'ivm'

        kwargs['wht_type'] = str(self.wht_type)
        kwargs['pscale_ratio'] = self.pixel_scale_ratio
        kwargs.pop('pixel_scale_ratio')

        for k, v in kwargs.items():
            if k in ['pixfrac', 'kernel', 'fillval', 'wht_type', 'pscale_ratio']:
                log.info('  using: %s=%s', k, repr(v))

        return kwargs


def _check_list_pars(vals, name, min_vals=None):
    if vals is None:
        return None
    if len(vals) != 2:
        raise ValueError(f"List '{name}' must have exactly two elements.")
    n = sum(x is None for x in vals)
    if n == 2:
        return None
    elif n == 0:
        if min_vals and sum(x >= y for x, y in zip(vals, min_vals)) != 2:
            raise ValueError(f"'{name}' values must be larger or equal to {list(min_vals)}")
        return list(vals)
    else:
        raise ValueError(f"Both '{name}' values must be either None or not None.")
