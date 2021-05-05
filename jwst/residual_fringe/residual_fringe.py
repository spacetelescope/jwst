#
#  Module for applying fringe correction
#

import numpy as np
from functools import partial

from .. import datamodels

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

from astropy.table import Table
from astropy.io import ascii

from . import utils
import astropy.units as u
from .diagnostic_plot import diagnostic_plot
from ..stpipe import Step

# remove matplotlib after testing 
import matplotlib.pyplot as plt


class ResidualFringeCorrection():

    def __init__(self,
                 input_model,
                 residual_fringe_reference_file,
                 regions_reference_file,
                 **pars):
        
        self.input_model = input_model.copy()
        self.model = input_model.copy()
        
        self.residual_fringe_reference_file = residual_fringe_reference_file
        self.regions_reference_file = regions_reference_file

        self.save_intermediate_results = pars['save_intermediate_results']
        self.transmission_level = int(pars['transmission_level'])
        # define how filenames are created
        self.make_output_path = pars.get('make_output_path',
                                        partial(Step._make_output_path,None))
                                                
        print('make output path',self.make_output_path)
        
        self.rfc_factors = None
        self.fit_mask = None
        self.weighted_pix_num = None
        self.rejected_fit = None
        self.weights_feat = None
        self.input_weights = None
        self.max_amp = None
        self.freq_table = None
        self.slice_map = None

        self.background_fit = None
        self.knot_locations = None
        
        self.band = None
        self.channel = None

        # remove after tested code - tries to make some plots and  store some output
        # output can be binary output - TODO need organize what data is saved.
        
        self.diagnostic_mode = True
        
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

        output_data = self.model.data.copy()
        output_data /= np.median(self.input_model.data)
        #print('normalizing the data', np.median(self.input_model.data))
        
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
        log.info("Detector {} ".format(detector,self.channels))
        
        self.input_weights = self.calc_weights()
        self.weights_feat = self.input_weights.copy()

        self.rfc_factors = np.zeros(self.input_model.data.shape)
        self.fit_mask = np.zeros(self.input_model.data.shape)
        self.weighted_pix_num = np.zeros(self.input_model.data.shape)
        self.rejected_fit =  np.zeros(self.input_model.data.shape)

        self.background_fit = np.zeros(self.input_model.data.shape)
        self.knot_locations = np.zeros(self.input_model.data.shape)
                
        allregions.close()

        # intermediate output product - Table 
        stat_table = Table(names=('Slice', 'mean', 'median', 'stddev', 'max', 'pmean', 'pmedian', 'pstddev', 'pmax'),
                           dtype=('i4', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8'))

        out_table = Table(names=('Slice', 'col', 'fringe', 'sn', 'pre_contrast', 'post_contrast',
                                 'periodogram_res', 'opt_fringes', 'peak_freq', 'freq_min', 'freq_max' ),
                           dtype=('i4', 'i4', 'i4', 'f8', 'f8', 'f8',
                                  'f8', 'f8', 'f8', 'f8', 'f8'))

        # intermediate output
        n_wav_samples = 1024
        qual_table = Table(names=('column', 'quality'),
                   units=('', ''),
                   dtype=('i4', '({},3)f8'.format(n_wav_samples)))

        
        for c in self.channels:
            log.info("Processing channel {}".format(c))
            (slices_in_band, xrange_channel, slice_x_ranges, all_slice_masks) = \
                utils.slice_info(slice_map,c)

            print(xrange_channel)
            
            ysize = self.input_model.data.shape[0]
            xsize = self.input_model.data.shape[1]
            y, x = np.mgrid[:ysize, :xsize]
            _, _, wave_map = self.input_model.meta.wcs(x,y)
            
            for n, ss in enumerate(slices_in_band):

                log.info(" Processing slice {} =================================".format(ss))
                log.info(" X ranges of slice {} {}".format(slice_x_ranges[n,1], slice_x_ranges[n,2]))
                # initialise the list to store correction quality numbers for slice
                correction_quality = []

                # use the mask to set all out-of-slice pixels to 0 in wmap and data
                # set out-of-slice pixels to 0 in arrays
                ss_data = all_slice_masks[n] * output_data.copy()
                
                ss_wmap = all_slice_masks[n] * wave_map
                ss_weight = all_slice_masks[n] * self.input_weights.copy()
                ss_mask = all_slice_masks[n]

                # get the freq_table info for this slice
                this_row = np.where(self.freq_table['slice'] == float(ss))[0][0]
                log.debug('Row in reference file for slice {}'.format(this_row))

                slice_row = self.freq_table[(self.freq_table['slice'] == float(ss))]
                ffreq = slice_row['ffreq'][0][0]
                dffreq = slice_row['dffreq'][0][0]
                min_nfringes = slice_row['min_nfringes'][0][0]
                max_nfringes = slice_row['max_nfringes'][0][0]
                min_snr = slice_row['min_snr'][0][0]
                pgram_res = slice_row['pgram_res'][0][0]
                
                #print(' freq',ffreq, dffreq, min_nfringes,max_nfringes)
                # cycle through the cols and fit the fringes
                for col in np.arange(slice_x_ranges[n,1], slice_x_ranges[n,2]):
                    col_data = ss_data[:, col]
                    col_wmap = ss_wmap[:, col]

                    valid = np.logical_and( (col_wmap>0), ~np.isnan(col_wmap))
                    # because of the curvature of the slices there can be large regions not falling on a column
                    num_good  = len(np.where(valid)[0])
                    #print('num of valid values in wavelength array',num_good)
                    # Need at least 50 pixels in column to proceed
                    
                    if num_good > 50:

                        test_flux = col_data[valid]
                        # Transform wavelength in micron to wavenumber in cm^-1.
                        col_wnum = 10000.0 / col_wmap

                        # do some checks on column to make sure there is reasonable signal. If the SNR < min_snr (CDP), pass
                        # determine SNR for this column of data
                        n = len(test_flux)
                        signal = np.median(test_flux)
                        noise = 0.6052697 * np.median(np.abs(2.0 * test_flux[2:n-2] - test_flux[0:n-4] - test_flux[4:n]))
                        
                        snr2  = 0.0 # initialize
                        if noise !=0:
                            snr2 = signal/noise
                        
                        # Sometimes can return nan, inf for bad data so include this in check
                        if (snr2 < min_snr[0]):
                                print('not fitting column',col)
                                pass
                        else:
                            log.debug("Fitting column{} ".format(col))
                            log.debug("SNR > {} ".format(min_snr[0]))
                            col_weight = ss_weight[:, col]
                            col_mask = ss_mask[:, col]
                            col_max_amp = np.interp(col_wmap, self.max_amp['Wavelength'], self.max_amp['Amplitude'])

                            # get the in-slice pixel indices for replacing in output later
                            idx = np.where(col_data > 0)

                            # BayesicFitting doesn't like 0s at data or weight array edges so set to small value
                            # replacing array 0s with arbitrarily low number
                            col_data[col_data <= 0] = 1e-08
                            col_weight[col_weight <= 0] = 1e-08

                            # check for off-slice pixels and send to be filled with interpolated/extrapolated wnums
                            # to stop BayesicFitting crashing, will not be fitted anyway
                            # finding out-of-slice pixels in column and filling

                            found_bad = np.logical_or( np.isnan(col_wnum), np.isinf(col_wnum))
                            num_bad  = len(np.where(found_bad)[0])

                            if num_bad> 0:
                                col_wnum[found_bad] = 0
                                col_wnum = utils.fill_wavenumbers(col_wnum)

                            # do feature finding on slice now column-by-column
                            log.debug(" starting feature finding")
                            
                            # find spectral features (env is spline fit of troughs and peaks) 
                            env, l_x, l_y, _, _, _ = utils.fit_envelope(np.arange(col_data.shape[0]), col_data)
                            mod = np.abs(col_data / env) -1

                            # given signal in mod find location of lines > col_max_amp * 2
                            weight_factors = utils.find_lines(mod, col_max_amp * 2)
                            weights_feat = col_weight * weight_factors
                            #print('size of weights_feat',weights_feat.shape)# - 1024
                            # iterate over the fringe components to fit, initialise pre-contrast, other output arrays
                            # in case fit fails
                            proc_data = col_data.copy()
                            proc_factors = np.ones(col_data.shape)
                            pre_contrast = 0.0
                            bg_fit = col_data.copy()
                            res_fringes = np.zeros(col_data.shape)
                            res_fringe_fit = np.zeros(col_data.shape)
                            res_fringe_fit_flag = np.zeros(col_data.shape)
                            wpix_num = 1024
                            # the reference file is set up with 2 values for ffreq but currently one one is used. The other value
                            # is a place holder and set to a small value
                            try:
                                for fn, ff in enumerate(ffreq):
                                    # ignore place holder fringes
                                    if ff > 1e-03:
                                        log.debug(' starting ffreq = {}'.format(ff))
                                        
                                        # check if snr criteria is met for fringe component, should always be true for fringe 1
                                        
                                        if snr2> min_snr[fn]:
                                            log.debug(' SNR > {}'.format(min_snr[fn]))
                                            log.debug(" fitting spectral baseline")
                                            bg_fit, bgindx, fitter = \
                                                utils.fit_1d_background_complex(proc_data, weights_feat,
                                                                                col_wnum, ffreq=ffreq[fn])

                                           # print('size of bg_fit and bgindx', bg_fit.shape, bgindx.shape)# 1024 and varying for bgindx
                                            
                                            # get the residual fringes as fraction of signal
                                            res_fringes = np.divide(proc_data, bg_fit, out=np.zeros_like(proc_data),
                                                                    where=bg_fit != 0)
                                            res_fringes = np.subtract(res_fringes, 1, where=res_fringes != 0)

                                            res_fringes *= np.where(col_weight > 1e-07, 1, 1e-08)
                                            
                                            # get the pre-correction contrast using fringe component 1
                                            if fn == 0:
                                                log.debug(" get pre-correction contrast")
                                                pre_contrast, quality  = utils.fit_quality(col_wnum,
                                                                                 res_fringes,
                                                                                 weights_feat,
                                                                                 ffreq[0],
                                                                                 dffreq[0])

                                                log.debug(" pre-correction contrast = {}".format(pre_contrast))

                                            # fit the residual fringes
                                            log.debug(" set up bayes ")
                                            res_fringe_fit, wpix_num, opt_nfringe, peak_freq, freq_min, freq_max= \
                                                utils.fit_1d_fringes_bayes_evidence(res_fringes,
                                                                                    weights_feat,
                                                                                    col_wnum,
                                                                                    ffreq[fn],
                                                                                    dffreq[fn],
                                                                                    min_nfringes=min_nfringes[fn],
                                                                                    max_nfringes=max_nfringes[fn],
                                                                                    pgram_res=pgram_res[fn])

                                            
                                            # check for fit blowing up, reset rfc fit to 0, raise a flag
                                            log.debug("check residual fringe fit for bad fit regions")
                                            res_fringe_fit, res_fringe_fit_flag = utils.check_res_fringes(res_fringe_fit,
                                                                                                    col_max_amp)
                                            
                                            # correct for residual fringes
                                            log.debug(" divide out residual fringe fit, get fringe corrected column")
                                            _, _, _, env, u_x, u_y = utils.fit_envelope(np.arange(res_fringe_fit.shape[0]),
                                                                                        res_fringe_fit)
                                            rfc_factors = 1 / (res_fringe_fit * (col_weight > 1e-05).astype(int) + 1)
                                            proc_data *= rfc_factors
                                            proc_factors *= rfc_factors

                                            # handle nans or infs that may exist
                                            proc_data[proc_data == np.inf] = 0
                                            proc_data = np.nan_to_num(proc_data)

                                            out_table.add_row((ss, col, fn ,snr2, pre_contrast, pre_contrast, pgram_res[fn],
                                                               opt_nfringe,peak_freq, freq_min, freq_max))


                                # define fringe sub after all fringe components corrections
                                fringe_sub = proc_data.copy()
                                rfc_factors = proc_factors.copy()

                                # get the new fringe contrast
                                log.debug(" analysing fit quality")
                                pbg_fit, pbgindx, pfitter = utils.fit_1d_background_complex(fringe_sub,
                                                                                            weights_feat,
                                                                                            col_wnum,
                                                                                            ffreq=ffreq[0])

                                # get the residual fringes as fraction of signal
                                fit_res = np.divide(fringe_sub, pbg_fit, out=np.zeros_like(fringe_sub),
                                                    where=pbg_fit != 0)
                                fit_res = np.subtract(fit_res, 1, where=fit_res != 0)
                                fit_res *= np.where(col_weight > 1e-07, 1, 1e-08)

                                
                                contrast, quality = utils.fit_quality(col_wnum,
                                                                      fit_res,
                                                                      weights_feat,
                                                                      ffreq,
                                                                      dffreq,
                                                                      save_results=self.save_intermediate_results)
                                print(type(quality), quality.shape)
                                
                                qual_table.add_row(col, quality)
                                correction_quality.append([contrast, pre_contrast])
                                log.debug(" residual contrast = {}".format(contrast))
                                
                               # if self.diagnostic_mode:

                                    # MASK1D plot
                                    #out_name = os.path.join(self.plot_dir, out_root + '_mask.pdf')
                                    #out_name = 'mask.pdf'
                                    #xdata, ydata = [[col, col]], [[1, 1022]]
                                    #diagnostic_plot('mask1d', 'columns', out_name,
                                    #                xdata=xdata, ydata=ydata, idata=[all_slice_masks[n]])

                                    # FITS1D plot
                                    #out_name = os.path.join(self.plot_dir, out_root + '_res-fringe-fits.pdf')
                                #    out_name = 'res-friinge-fits.pdf'
                                #    xdata = [col_wnum]
                                 #   ydata = [col_data, bg_fit, res_fringes, col_weight, res_fringe_fit, weights_feat,
                                 #            fringe_sub]
                                 #   diagnostic_plot('fits1d', 'columns', out_name, xdata=xdata, ydata=ydata,
                                 #                   save_data=True)

                                # replace the corrected in-slice column pixels in the data_cor array
                                log.debug(" updating the trace pixels in the output")
                                output_data[idx, col] = fringe_sub[idx]
                                self.rfc_factors[idx, col] = rfc_factors[idx]
                                self.fit_mask[idx, col] = np.ones(1024)[idx]
                                self.weights_feat[idx, col] = weights_feat[idx]
                                self.weighted_pix_num[idx, col] = np.ones(1024)[idx] * (wpix_num/1024)
                                self.rejected_fit[idx, col] = res_fringe_fit_flag[idx]
                                self.background_fit[idx,col] = bg_fit[idx]
                                #self.knot_locations[idx,col] = bgindx

                            except Exception as e:
                                log.warning(" Skipping col={}".format(col, ss))
                                log.warning(' %s' % (str(e)))

                del ss_data, ss_wmap, ss_weight  # end on column

                # asses the fit quality statistics and make plot if in diagnostic mode
                log.debug(" analysing fit statistics")
                if len(correction_quality) > 0:
                    contrasts = np.asarray(correction_quality)[:,0]
                    pre_contrasts = np.asarray(correction_quality)[:, 1]
                    mean, median, stddev, fmax = utils.fit_quality_stats(contrasts)
                    pmean, pmedian, pstddev, pfmax = utils.fit_quality_stats(pre_contrasts)

                    # add the statistics to the Table
                    stat_table.add_row((ss, mean, median, stddev, fmax, pmean, pmedian, pstddev, pfmax))

            del slice_x_ranges, all_slice_masks, slices_in_band, wave_map # end of channel

        log.info("Processing complete")

        # add units back to output data
        log.info(" adding units back to output array")
        output_data *= np.median(self.input_model.data)
        self.model.data = output_data
        del output_data
        
        if self.save_intermediate_results:
            self.logger.debug(" saving the output data")
            stat_table_name = self.make_output_path(
                basepath=self.input_model.meta.filename,
                suffix='stat_table',ext='.ecvs')
            print('Stat table name',stat_table_name)
            ascii.write(stat_table,stat_table_name, format='ecsv',fast_writer=False)

            out_table_name = self.make_output_path(
                basepath=self.input_model.meta.filename,
                suffix='out_table',ext='.ecvs')
            print('out table name',out_table_name)
            ascii.write(out_table,out_table_name, format='ecsv',fast_writer=False)


        #save_data_astropy(self.input, self.input_type, self.output_dir, self.output_data,
        #                  self.rfc_factors, self.fit_mask, self.weights_feat,
        #                  self.weighted_pix_num, stat_table, self.rejected_fit)

        #self.logger.info("Output saved")
            

        # garbage collect when finished
        #gc.collect()

        return self.model



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

        # replace infs and nans in weights with 0
        weights[weights == np.inf] = 0
        weights[np.isnan(weights)] = 0
        print('value of weights',weights[0:100,473])
        return weights
