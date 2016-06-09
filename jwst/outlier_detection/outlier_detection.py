import time
import numpy as np
from collections import OrderedDict

from .. import datamodels
from .. import assign_wcs
from .. import resample

from . import flag_cr
from . import blot_median
from . import create_median

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

DEFAULT_DOMAIN = {'lower':None,'upper':None,'includes_lower':True, 'includes_upper':False}

class OutlierDetection(object):
    """
    This is the controlling routine for the outlier detection process.
    It loads and sets the various input data and parameters needed by
    the various functions and then controls the operation of this process
    through all the steps used for the detection.

    Notes
    -----
    This routine performs the following operations::

      1. Extracts parameter settings from input model and merges
         them with any user-provided values
      2. Resamples all input images into grouped observation mosaics.
      3. Creates a median image from all grouped observation mosaics.
      4. Blot median image to match each original input image.
      5. Perform statistical comparison between blotted image and original image
         to identify outliers.
      6. Updates input data model DQ arrays with mask of detected outliers.

    """
    outlierpars = {'kernel':'square','pixfrac':1.0,'resample_bits':None,
                        'fillval':'INDEF','wht_type':'exptime',
                    'nlow': 0, 'nhigh':1,
                        'hthresh':None, 'lthresh':None,
                        'nsigma': '4 3', 'maskpt':0.7,
                    'grow': 1, 'ctegrow':0, 'snr': "4.0 3.0", 
                        'scale': "0.5 0.4", 'backg': 0
                }

    def __init__(self, input_models, ref_filename=None, to_file=False, **pars):
        """
        Parameters
        ----------
        input_models : list of objects or str
            list of data models as ModelContainer or ASN file,
            one data model for each input image

        drizzled_models : list of objects
            ModelContainer containing drizzled grouped input images
      
        to_file : bool
            Control whether or not to write out intermediate results to files

        """
        self.input_models = input_models
        self.ref_filename = ref_filename
        self.to_file = to_file

        self.num_groups = len(self.input_models.group_names)

        self.outlierpars = {}
        if 'outlierpars' in ref_filename:
            self._get_outlier_pars()
        self.outlierpars.update(pars)


    def _get_outlier_pars(self):
        """ Extract outlier detection parameters from reference file
        """
        # start by interpreting input data models to define selection criteria
        input_dm = self.input_models[0]
        filtname = input_dm.meta.instrument.filter

        ref_model = datamodels.OutlierParsModel(self.ref_filename['outlierpars'])

        # look for row that applies to this set of input data models
        # NOTE:
        #  This logic could be replaced by a method added to the DrizParsModel object
        #  to select the correct row based on a set of selection parameters
        row = None
        outlierpars = ref_model.outlierpars_table

        filter_match = False # flag to support wild-card rows in outlierpars table
        for n,filt,num in zip(range(1,outlierpars.numimages.shape[0]+1),outlierpars.filter,
                            outlierpars.numimages):
            # only remember this row if no exact match has already been made for
            # the filter. This allows the wild-card row to be anywhere in the
            # table; since it may be placed at beginning or end of table.

            if filt == "ANY" and not filter_match and self.num_groups >= num:
                row = n
            # always go for an exact match if present, though...
            if filtname == filt and self.num_groups>= num:
                row = n
                filter_match = True

        # With presence of wild-card rows, code should never trigger this logic
        if row is None:
            log.error("No row found in %s that matches input data.",self.ref_filename)
            raise ValueError

        # read in values from that row for each parameter
        for kw in list(self.outlierpars.keys()):
            self.outlierpars[kw] = ref_model['outlierpars_table.{0}'.format(kw)]

    def do_detection (self):
        """ Perform drizzling operation on input images's to create a new output

        """
        pars = self.outlierpars
        
        # Start by creating resampled/mosaic images for each group of exposures
        sdriz = resample.resample.ResampleData(self.input_models, single=True, **pars)
        sdriz.do_drizzle(**pars)
        drizzled_models = sdriz.output_models
        if self.to_file:
            log.info("Saving resampled grouped exposures to disk...")
            drizzled_models.save(None)

        # Initialize intermediate products used in the outlier detection
        median_model = datamodels.ImageModel(init=drizzled_models[0].data.shape)
        median_model.meta = drizzled_models[0].meta # provide median with initial metadata
        base_filename = self.input_models[0].meta.filename
        median_filename = '_'.join(base_filename.split('_')[:2]+['median.fits'])
        median_model.meta.filename = median_filename

        # Perform median combination on set of drizzled mosaics
        drizzle_groups_sci = [i.data for i in drizzled_models]
        drizzle_groups_wht = [i.wht for i in drizzled_models]
        median_model.data = create_median.do_median(drizzle_groups_sci,
                                        drizzle_groups_wht,
                                        **pars)
        if self.to_file:
            log.info("Writing out MEDIAN image to: {}".format(median_model.meta.filename))
            median_model.save(median_model.meta.filename)
        # Blot the median image back to recreate each input image specified in
        # the original input list/ASN/ModelContainer
        blot_models = blot_median.do_blot(median_model, self.input_models,
                     **pars)
        if self.to_file:
            log.info("Writing out BLOT input images...")
            blot_models.save(None)
            
        # Perform outlier detection using statistical comparisons between
        # original input images and their blotted-median images
        flag_cr.do_detection(self.input_models, blot_models, 
            self.ref_filename,**pars)

        # clean-up (just to be explicit about being finished with these results)
        del median_model, blot_models
