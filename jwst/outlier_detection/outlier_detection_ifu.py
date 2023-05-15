"""Class definition for performing outlier detection on IFU data."""

import numpy as np
from stsci.image import median
from astropy.stats import sigma_clipped_stats

from stdatamodels.jwst import datamodels
from scipy.signal import medfilt

from .outlier_detection import OutlierDetection
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
        OutlierDetection.__init__(self, input_models,
                                  reffiles=reffiles, **pars)

    def _find_detector_parameters(self):
        print('Instrument', self.inputs[0].meta.instrument.name.upper())
        
        if self.inputs[0].meta.instrument.name.upper() == 'MIRI':
            diffaxis = 1
        elif self.inputs[0].meta.instrument.name.upper() == 'NIRSPEC':
            diffaxis = 0
            
        ny,nx = self.inputs[0].data.shape
        print(' Shape of array', ny, nx)
        return (diffaxis, ny, nx)
        
    def do_detection(self):

        log.info("Flagging outliers")

        self.build_suffix(**self.outlierpars)
        print('Suffix ', self.resample_suffix)
              
        save_intermediate_results = \
            self.outlierpars['save_intermediate_results']

        kernel_size = self.outlierpars['kernel_size']
        sizex, sizey = [int(val) for val in kernel_size.split()]
        kern_size = np.zeros(2, dtype = int)
        kern_size[0]=  sizex
        kern_size[1] = sizey
        print(kern_size)
        
        print(type(kern_size))
        
        threshold_percent = self.outlierpars['threshold_percent']

        (diffaxis, ny, nx) = self._find_detector_parameters()
        
        # set up array to hold group differences 
        n = len(self.inputs)
        diffarr = np.zeros([n,ny,nx])

        for i, model in enumerate(self.inputs):
            sci = model.data
            dq = model.dq
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

            diffarr[i,:,:]=combdiff

        # minarr final minimum combined differences, size: ny X nx
        minarr=np.nanmin(diffarr,axis=0)

        # Normalise the differences to a local median image to deal with ultra-bright sources
        norm=medfilt(minarr,kernel_size=kern_size)
        nfloor=np.nanmedian(minarr)/3
        norm[norm < nfloor] = nfloor # Ensure we never divide by a tiny number
        minarr_norm=minarr/norm
        # Percentile cut of the central region (cutting out weird detector edge effects)
        pctmin=np.nanpercentile(minarr_norm[4:ny-4,4:nx-4],threshold_percent)
        print('Percentile min: ',threshold_percent,pctmin)

        if save_intermediate_results:
              #model.meta.filename = self.make_output_path(
              #    basepath=model.meta.filename,
              #    suffix='diffarr')
              #log.info("Writing out (single) IFU cube {}".format(model.meta.filename))
              #model.save(model.meta.filename)
              
              diffarr_filename = self.make_output_path(
                  basepath=model.meta.filename,
                  suffix='diffarr')
              hdu=fits.PrimaryHDU(diffarr)
              hdu.writeto(diffarr_filename,overwrite=True)

              minarr_filename = self.make_output_path(
                  basepath=model.meta.filename,
                  suffix='minarr')
              hdu=fits.PrimaryHDU(minarr)
              hdu.writeto(minarr_filename,overwrite=True)

              norm_filename = self.make_output_path(
                  basepath=model.meta.filename,
                  suffix='norm')
              hdu=fits.PrimaryHDU(norm)
              hdu.writeto(norm_filename,overwrite=True)
              
              minarr_norm_filename = self.make_output_path(
                  basepath=model.meta.filename,
                  suffix='minarr_norm')
              hdu=fits.PrimaryHDU(minarr_norm)
              hdu.writeto(minarr_norm_filename,overwrite=True)
              
        
        # Flag everything above this percentile value
        indx=np.where(minarr_norm > pctmin)
        print('number of flagged pixels', len(indx[0]))
        print(indx)
        # Update in place dq flag
        for i, model in enumerate(self.inputs):
              count_existing = np.count_nonzero(model.dq & dqflags.pixel['DO_NOT_USE'])
              sci = model.data
              dq = model.dq
              sci[indx]=np.nan
              dq[indx] = np.bitwise_or(dq[indx], dqflags.pixel['DO_NOT_USE'])
              dq[indx] = np.bitwise_or(dq[indx], dqflags.pixel['OUTLIER'])
              model.data = sci
              model.dq = dq

              #hdu.writeto(files[ii].replace(calstr,outstr),overwrite=True)
              # Update the DQ array in the input image.
              #sci_image.dq = np.bitwise_or(sci_image.dq, cr_mask * (DO_NOT_USE | OUTLIER))

              # Report number (and percent) of new DO_NOT_USE pixels found
              count_outlier = np.count_nonzero(dq & dqflags.pixel['DO_NOT_USE'])
              count_added = count_outlier - count_existing
              percent_cr = count_added / (model.data.shape[0] * model.data.shape[1]) * 100
              log.info(f"New pixels flagged as outliers: {count_added} ({percent_cr:.2f}%)")
              

        del diffarr
        

class ErrorWrongInstrument(Exception):
    """ Raises an exception if the instrument is not MIRI or NIRSPEC
    """
    pass
