#
#  Module for applying fringe correction
#
from .. import datamodels
import numpy as np
import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

from . import residual_fringe_tools

class ResidualFringeCorrection():

    def __init__(self,
                 input_model,
                 residual_fringe_reference_file,
                 regions_reference_file,
                 transmission_level):
        
        self.input_model = input_model
        self.residual_fringe_reference_file = residual_fringe_reference_file
        self.regions_reference_file = regions_reference_file
        self.transmission_level = int(transmission_level)
        self.output_data = None
        self.rfc_factors = None
        self.fit_mask = None
        self.weighted_pix_num = None
        self.rejected_fit = None
        self.weights_feat = None
        self.input_weights = None
        self.max_amp = None
        self.freq_table = None
        self.slice_map = None

        self.band = None
        self.channel = None
        
    def do_correction(self):
        
        """
        Short Summary
        -------------
        Residual Fringe-correct a JWST data model using a residual fringe model

        Parameters
        ----------
        input_model: JWST data model
        input science data model to be residual fringe-corrected

        residual_fringe_model: JWST data model
        data model containing residual fringe correction

        Returns
        -------
        output_model: JWST data model
        residual fringe-corrected science data model

        """

        # Load the fringe reference file
        residual_fringe_model = datamodels.ResidualFringeModel(self.residual_fringe_reference_file)

        # read in the band
        band = self.input_model.meta.instrument.band.lower()
        if band == 'short':
            residual_fringe_table = residual_fringe_model.rfc_freq_short_table
        elif band == 'medium':
            residual_fringe_table = residual_fringe_model.rfc_freq_medium_table
        else:
            residual_fringe_table = residual_fringe_model.rfc_freq_long_table

        self.max_amp = residual_fringe_model.max_amp_table
        residual_fringe_model.close()

        self.freq_table = residual_fringe_table
        # Read in the regions reference file
        # Use throughput array defined by self.transmission_level
        allregions = datamodels.RegionsModel(self.regions_reference_file)
        self.transmission_level= int(self.transmission_level/10)
        print(allregions.regions.shape)
        slice_map = (allregions.regions)[self.transmission_level-1, :, :].copy()
        log.info(" Using {} throughput threshhold.".format(self.transmission_level))

        self.slice_map = slice_map

        # set up the channels for the detector
        detector = self.input_model.meta.instrument.detector.lower()        

        if 'short' in detector:
            self.channels = [1, 2]
        elif 'long' in detector:
            self.channels = [3, 4]
        print('detector',detector,self.channels)
        
        self.input_weights = self.calc_weights()
        self.output_data = self.input_model.data.copy()
        self.weights_feat = self.input_weights.copy()

        self.rfc_factors = np.zeros(self.output_data.shape)
        self.fit_mask = np.zeros(self.input_model.data.shape)
        self.weighted_pix_num = np.zeros(self.input_model.data.shape)
        self.rejected_fit =  np.zeros(self.input_model.data.shape)
        allregions.close()

        for c in self.channels:
            log.info("Processing channel {}".format(c))
            result = residual_fringe_tools.slice_info(slice_map,c)
            (slices_in_band, xrange_channel, all_slice_masks) = result
            
            ysize = self.input_model.data.shape[0]
            xsize = self.input_model.data.shape[1]
            #y, x = np.mgrid[:ysize, xrange_channel[0]: xrange_channel[1]]
            y, x = np.mgrid[:ysize, :xsize]
            _, _, wave_map = self.input_model.meta.wcs(x,y)

            for n, ss in enumerate(slices_in_band):

                log.info(" Processing slice {} =================================".format(ss))

                # initialise the list to store correction quality numbers for slice
                correction_quality = []

                # use the mask to set all out-of-slice pixels to 0 in wmap and data
                # set out-of-slice pixels to 0 in arrays
                ss_data = all_slice_masks[n] * self.input_model.data.copy()
                ss_wmap = all_slice_masks[n] * wave_map
                ss_weight = all_slice_masks[n] * self.input_weights.copy()
                ss_mask = all_slice_masks[n]

                # get the freq_table info for this slice
                slice_row = self.freq_table[(self.freq_table['Slice'] == str(ss))]
                print(' Using slice freq info:\n {}'.format(slice_row))
                ffreq = slice_row['ffreq'][0][0]
                dffreq = slice_row['dffreq'][0][0]
                min_nfringes = slice_row['min_nfringes'][0][0]
                max_nfringes = slice_row['max_nfringes'][0][0]
                min_snr = slice_row['min_snr'][0][0]
                pgram_res = slice_row['pgram_res'][0][0]


            
        result = self.apply_residual_fringe()


    #output_model.meta.cal_step.residual_fringe = 'COMPLETE'
    #return output_model


    def apply_residual_fringe(self):
        """
        Short Summary
        -------------
        Residual Fringe Correction:
        
        Parameters
        ----------
        input_model: JWST data model
        input science data model to be fringe-corrected

        residual_fringe_table:  

        Returns
        -------
        output_model: JWST data model
        input science data model which has been fringe-corrected
        """

        # Initialize the output model as a copy of the input
        #output_model = self.input_model.copy()

        # Fringe-correct data and error arrays, applying correction only
        # to pixels having good calibration values
        #fringe_data = fringe.data
        #fringe_data[np.isnan(fringe_data)] = 1.0
        #output_model.data /= fringe_data
        #output_model.err /= fringe_data
        print('stopping')
        exit(1)


        return output_model


    def calc_weights(self):

        """Make a weights array based on flux. This is a placeholder function, for now just returns a normalised
        flux array to use as weights array. This is because any smoothing results in incorrect rfc around emission lines.
        This can be changed in the future if need be.

        :Parameters:

        flux:  numpy array, required
        the 1D array of fluxes

        :Returns:

        weights array which is just a copy of the flux array, normalised by the mean

        """
         
        weights = np.zeros(self.input_model.data.shape)
        for c in np.arange(weights.shape[1]):
            flux_1d = self.input_model.data[:,c]
            w =  flux_1d/np.mean(flux_1d)
            weights[:, c] = w
            #print(c,w)

        # replace infs and nans in weights with 0
        weights[weights == np.inf] = 0
        weights[np.isnan(weights)] = 0
        return weights
