"""Class definition for performing outlier detection on IFU data."""

import numpy as np
from stsci.image import median
from astropy.stats import sigma_clipped_stats

from stdatamodels.jwst import datamodels
from scipy.signal import medfilt, medfilt2d

from .outlier_detection import OutlierDetection
from jwst.datamodels import ModelContainer
from stdatamodels.jwst.datamodels import dqflags
import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["OutlierDetectionIFU"]


class OutlierDetectionIFU(OutlierDetection):
    """Sub-class defined for performing outlier detection on IFU data.

    This is the controlling routine for the outlier detection process.
    It loads and sets the various input data and parameters needed to flag
    outliers.  Pixel are flagged as outliers based on the MINIMUM difference 
    a pixel has with its neighbor across all the input cal files. 

    Notes
    -----
    This routine performs the following operations::

      1. Extracts parameter settings from input ModelContainer and merges
         them with any user-provided values
      2. Loop over cal files
         a. read in science data 
         b. Store computed neighbor differences for all the pixels.
            The neighbor pixel  differences are defined by the dispersion axis.
            For MIRI (disp axis = 1) the neighbors to find differences  are to the left and right of pixel
            For NIRSpec (disp axis = 0) the neighbors to find the differences are above and below the pixel 

      3. For each input file store the  minimum of the two diferences
      4. Comparing all the differences from all the input data find the minimum difference
      5. Normalize minimum difference to local median of difference array
      6. select outliers by flagging those normailzed minimum values > thershold_percent

      7. Updates input ImageModel DQ arrays with mask of detected outliers.

    """

    def __init__(self, input_models, reffiles=None, **pars):
        """Initialize class for IFU data processing.

        Parameters
        ----------
        input_models : ~jwst.datamodels.ModelContainer, str
            list of data models as ModelContainer or ASN file,
            one data model for each input 2-D ImageModel

        reffiles : dict of `~stdatamodels.jwst.datamodels.JwstDataModel`
            Dictionary of datamodels.  Keys are reffile_types.


        """
        OutlierDetection.__init__(self, input_models,reffiles=reffiles, **pars)

    def create_optional_results_model(self, opt_info):
        """
        Creates an OutlierOutputModel from the computed arrays from outlier detection on IFU data.

        Parameter
        ---------
        input_model: ~jwst.datamodels.RampModel

        opt_info: tuple
        The output arrays needed for the OultierOutputModel.

        Return
        ---------
        opt_model: OultierOutputModel
        The optional OutlierOutputModel to be returned from the outlier_detection_ifu step.
        """
        (diffarr, minarr, normarr, minnorm) = opt_info
        
        opt_model = datamodels.OutlierIFUOutputModel(
            diffarr=diffarr,
            minarr=minarr,
            normarr=normarr,
            minnorm=minnorm)
        return opt_model


    def _find_detector_parameters(self):
        if self.input_models[0].meta.instrument.name.upper() == 'MIRI':
            diffaxis = 1
        elif self.input_models[0].meta.instrument.name.upper() == 'NIRSPEC':
            diffaxis = 0
            
        ny,nx = self.inputs[0].data.shape
        return (diffaxis, ny, nx)

    def convert_inputs(self):
        self.input_models = self.inputs
        self.converted = False
                        
    def do_detection(self):
        """Split data by detector to find outliers."""
        self.convert_inputs()

        self.build_suffix(**self.outlierpars)
        save_intermediate_results = \
            self.outlierpars['save_intermediate_results']

        kernel_size = self.outlierpars['kernel_size']
        sizex, sizey = [int(val) for val in kernel_size.split()]
        kern_size = np.zeros(2, dtype = int)
        kern_size[0]=  sizex
        kern_size[1] = sizey

        threshold_percent = self.outlierpars['threshold_percent']
        (diffaxis, ny, nx) = self._find_detector_parameters()

        nfiles = len(self.input_models)
        detector=np.empty(nfiles,dtype='<U15')
        for i, model in enumerate(self.input_models):
            detector[i] = model.meta.instrument.detector.lower()

        exptype = self.input_models[0].meta.exposure.type
        log.info("Performing IFU outlier_detection for exptype {}".format(
                 exptype))
        # How many unique values of detector?
        uq_det=np.unique(detector)
        ndet=len(uq_det)
        for idet in range(ndet):
            indx=(np.where(detector == uq_det[idet]))[0]
            ndet_files = int(len(indx))
            self.flag_outliers(idet, uq_det, ndet_files,
                               diffaxis, nx, ny,
                               kern_size, threshold_percent,
                               save_intermediate_results)

        # send input_models back - that is what is returned from outlier_detection.py
        self.detect_outliers_ifu(self.input_models)            
        
    def flag_outliers(self, idet, uq_det, ndet_files,
                      diffaxis, nx, ny,
                      kern_size, threshold_percent,
                      save_intermediate_results):
        """Flag outlier pixels in DQ of input images."""
        
        # set up array to hold group differences
        diffarr = np.zeros([ndet_files, ny, nx])
        print('working on detector', uq_det[idet])

        j = 0 
        for i, model in enumerate(self.input_models):
            sci = model.data
            dq = model.dq
            detector = model.meta.instrument.detector.lower()
            # only use data from the same detector
            if detector == uq_det[idet]:
                bad = np.where(np.bitwise_and(dq, dqflags.pixel['DO_NOT_USE']).astype(bool))
                sci[bad] = np.nan

                # Compute left and right differences (MIRI dispersion axis = 1)
                # For NIRSpec dispersion axis = 0, these differences are top, bottom
                # prepend = 0 has the effect of keeping the same shape as sci and
                # for MIRI data (disp axis = 1) the first column = sci data
                # OR
                # for NIRSpec data (disp axis = 0) the first row = sci data
                leftdiff=np.diff(sci,axis=diffaxis,prepend=0)

                flip=np.flip(sci,axis=diffaxis)
                rightdiff=np.diff(flip,axis=diffaxis,prepend=0)
                rightdiff=np.flip(rightdiff,axis=diffaxis)

                # Combine left and right differences with minimum of the abs value
                # to avoid artifacts from bright edges
                comb=np.zeros([2,ny,nx])
                comb[0,:,:]=np.abs(leftdiff)
                comb[1,:,:]=np.abs(rightdiff)
                combdiff=np.nanmin(comb,axis=0)

                diffarr[j,:,:]=combdiff
                j = j + 1

        # minarr final minimum combined differences, size: ny X nx
        minarr=np.nanmin(diffarr,axis=0)
        print(minarr[502,478])
        print(minarr[500:504, 477:481])
        # Normalise the differences to a local median image to deal with ultra-bright sources
        print(type(minarr))
        normarr=medfilt(minarr,kernel_size=kern_size)
        print('normarr', normarr[502,478])
        print(normarr[500:504, 477:481])
        nfloor=np.nanmedian(minarr)/3
        print(nfloor)
        normarr[normarr < nfloor] = nfloor # Ensure we never divide by a tiny number
        print('normarr', normarr[502,478])
        minarr_norm=minarr/normarr
        
        # Percentile cut of the central region (cutting out weird detector edge effects)
        pctmin=np.nanpercentile(minarr_norm[4:ny-4,4:nx-4],threshold_percent)
        print('pctmin', pctmin)
        log.info("Flag pixels with values above {} {}: ".format(threshold_percent,pctmin))
        # Flag everything above this percentile value
        indx=np.where(minarr_norm > pctmin)
        print(minarr_norm[500:504, 477:481])
        log.info("Number of outlier pixels flagged: {}".format(
                len(indx[0])))
            
        if save_intermediate_results:
            detector_name = uq_det[idet]
            opt_info = (diffarr, minarr, normarr, minarr_norm)        
            opt_model = self.create_optional_results_model(opt_info)
            opt_model.meta.filename = self.make_output_path(
                basepath=self.input_models.meta.asn_table.products[0].name,
                suffix=detector_name + '_outlier_output')
            log.info("Writing out intermediate outlier file {}".format(opt_model.meta.filename))
            opt_model.save(opt_model.meta.filename)
              
        del diffarr

        # Update DQ flag
        for i in range(len(self.input_models)):
            model = datamodels.open(self.input_models[i])
            sci = model.data
            dq = model.dq
            detector = model.meta.instrument.detector.lower()
            # only use data from the same detector
            if detector == uq_det[idet]:
                count_existing = np.count_nonzero(dq & dqflags.pixel['DO_NOT_USE'])

                sci[indx]=np.nan
                dq[indx] = np.bitwise_or(dq[indx], dqflags.pixel['DO_NOT_USE'])
                dq[indx] = np.bitwise_or(dq[indx], dqflags.pixel['OUTLIER'])
           
                # Report number (and percent) of new DO_NOT_USE pixels found
                count_outlier = np.count_nonzero(dq & dqflags.pixel['DO_NOT_USE'])
                count_added = count_outlier - count_existing
                percent_cr = count_added / (model.data.shape[0] * model.data.shape[1]) * 100
                log.info(f"New pixels flagged as outliers: {count_added} ({percent_cr:.2f}%)")
                # update model
                model.dq = dq
                model.data = sci
                self.input_models[i] = model
            model.close()
