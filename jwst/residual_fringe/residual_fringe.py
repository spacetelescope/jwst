#
#  Module for applying fringe correction
#
from .. import datamodels
import numpy as np
import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

from . import residual_fringe_tools
import astropy.units as u
from specutils import Spectrum1D
from specutils.analysis import snr_derived

class ResidualFringeCorrection():

    def __init__(self,
                 input_model,
                 residual_fringe_reference_file,
                 regions_reference_file,
                 transmission_level):
        
        self.input_model = input_model
        self.output_model =None
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
        
        # normalise the output_data to remove units

        self.output_data  = self.input_model.data.copy()
        self.output_data /= np.median(self.input_model.data)
        
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
            (slices_in_band, xrange_channel, all_slice_masks) = residual_fringe_tools.slice_info(slice_map,c)
            print(all_slice_masks.shape)
            print(xrange_channel)
            
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
                ss_data = all_slice_masks[n] * self.output_data.copy()
                print(' columns 450 to 500 in data*mask',ss_data[0,450:500])
                
                ss_wmap = all_slice_masks[n] * wave_map
                ss_weight = all_slice_masks[n] * self.input_weights.copy()
                ss_mask = all_slice_masks[n]

                
                print('weights',ss_weight[0,450:500])
                # get the freq_table info for this slice
                this_row = np.where(self.freq_table['slice'] == float(ss))[0][0]
                print('this_row',this_row)
                
                ffreq = self.freq_table['ffreq'][this_row][0][0]
                dffreq = self.freq_table['dffreq'][this_row][0][0]
                print(' freq',ffreq, dffreq)
                min_nfringes = self.freq_table['min_nfringes'][0][0]
                max_nfringes = self.freq_table['max_nfringes'][0][0]
                min_snr = self.freq_table['min_snr'][0][0]
                pgram_res = self.freq_table['pgram_res'][0][0]
                # cycle through the cols and fit the fringes
                for col in np.arange(ss_data.shape[1]):
                    # set a root name for outputs
                    #out_root = self.input_detector + '_' + self.input_band + '_CH' + str(c) + \
                    #           '_SL' + str(ss) + '_COL' + str(col)

                    print('ss_data',ss_data.shape)
                    print('col',col)
                    # get the column data
                    col_data = ss_data[:, col]
                    col_wmap = ss_wmap[:, col]
                    col_weight = ss_weight[:, col]
                    col_mask = ss_mask[:, col]
                    col_max_amp = np.interp(col_wmap, self.max_amp['Wavelength'], self.max_amp['Amplitude'])
                    #print(col_max_amp.shape)
                    

                    # Transform wavelength in micron to wavenumber in cm^-1.
                    col_wnum = 10000.0 / col_wmap

                    #print(col_data.shape)
                    #print(col_wnum.shape)

                    
                    # do some checks on column to make sure there is reasonable signal. If the SNR < min_snr (CDP), pass
                    # to use the method http://www.stecf.org/software/ASTROsoft/DER_SNR/ in specutils, need to
                    # use a Spectrum1D option. Use arbitraty units
                    test_wnum = col_wnum.copy()
                    test_flux = col_data.copy()
                    check_spectrum = Spectrum1D(flux=test_flux[::-1] * u.Jy,
                                               spectral_axis=test_wnum[::-1] / u.cm)

                    
                    # determine SNR
                    snr = snr_derived(check_spectrum, region=None)
                    if snr.value > min_snr[0]:
                        print('check spectrum shape',check_spectrum.shape)
                        print(col_wnum[0:100])
                        print(col_data[0:100])
                        print(check_spectrum[0:100])
                        exit(-1)

                    # Sometimes can return nan, inf for bad data so include this in check
                    if (snr.value < min_snr[0]) or (str(snr.value) == 'inf') or (str(snr.value) == 'nan'):
                        print('not fitting colum')
                        pass
                    else:
                        log.info("Fitting column{}".format(col))



            
        #result = self.apply_residual_fringe()


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
