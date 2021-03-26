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

        # remove after tested code
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
            (slices_in_band, xrange_channel, slice_x_ranges, all_slice_masks) = \
                residual_fringe_tools.slice_info(slice_map,c)
            print(all_slice_masks.shape)
            print(xrange_channel)
            
            ysize = self.input_model.data.shape[0]
            xsize = self.input_model.data.shape[1]
            y, x = np.mgrid[:ysize, :xsize]
            _, _, wave_map = self.input_model.meta.wcs(x,y)
            
            for n, ss in enumerate(slices_in_band):

                log.info(" Processing slice {} =================================".format(ss))
                print('x range slice', slice_x_ranges[n,:])
                # initialise the list to store correction quality numbers for slice
                correction_quality = []

                # use the mask to set all out-of-slice pixels to 0 in wmap and data
                # set out-of-slice pixels to 0 in arrays
                ss_data = all_slice_masks[n] * self.output_data.copy()
                
                ss_wmap = all_slice_masks[n] * wave_map
                ss_weight = all_slice_masks[n] * self.input_weights.copy()
                ss_mask = all_slice_masks[n]

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
                for col in np.arange(slice_x_ranges[n,1], slice_x_ranges[n,2]):
                    col_data = ss_data[:, col]
                    col_wmap = ss_wmap[:, col]

                    valid = np.logical_and( (col_wmap>0), ~np.isnan(col_wmap))
                    num_good  = len(np.where(valid)[0])
                    #num_good = len(good_index[0])
                    print('num of valid values in wavelength array',num_good)

                    if num_good > 25:

                        test_flux = col_data[valid]
                        # Transform wavelength in micron to wavenumber in cm^-1.
                        col_wnum = 10000.0 / col_wmap
                        test_wnum = col_wnum[valid]

                        # do some checks on column to make sure there is reasonable signal. If the SNR < min_snr (CDP), pass
                        # to use the method http://www.stecf.org/software/ASTROsoft/DER_SNR/ in specutils, need to
                        # use a Spectrum1D option. Use arbitraty units
                        check_spectrum = Spectrum1D(flux=test_flux[::-1] * u.Jy,
                                                    spectral_axis=test_wnum[::-1] / u.cm)

                    
                        # determine SNR for this column of data
                        snr = snr_derived(check_spectrum, region=None)
                        print('SNR value',snr.value)
                        
                        # Sometimes can return nan, inf for bad data so include this in check
                        if (snr.value < min_snr[0]) or (str(snr.value) == 'inf') or (str(snr.value) == 'nan'):
                                print('not fitting colum')
                                pass
                        else:
                            log.info("Fitting column{} ".format(col))

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
                                col_wnum = residual_fringe_tools.fill_wavenumbers(col_wnum)

                            # do feature finding on slice now column-by-column
                            log.info(" starting feature finding")
                            # require at least 50 pixels
                            if (col_mask > 0).sum() > 50:
                                # find spectral features
                                env, l_x, l_y, _, _, _ = residual_fringe_tools.fit_envelope(np.arange(col_data.shape[0]), col_data)
                                mod = np.abs(col_data / env) - 1
                                weight_factors = residual_fringe_tools.find_lines(mod, col_max_amp * 2)
                                weights_feat = col_weight * weight_factors

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
                            try:
                                for fn, ff in enumerate(ffreq):
                                    # ignore place holder fringes
                                    if ff > 1e-03:
                                        log.info(' starting ffreq = {}'.format(ff))
                                        # check if snr criteria is met for fringe component, should always be true for fringe 1
                                        if snr.value > min_snr[fn]:
                                            log.info(' SNR > {}'.format(min_snr[fn]))
                                            log.info(" fitting spectral baseline")
                                            bg_fit, bgindx, fitter = \
                                                residual_fringe_tools.fit_1d_background_complex(proc_data, weights_feat,
                                                                                                col_wnum, ffreq=ffreq[fn])

                                            # get the residual fringes as fraction of signal
                                            res_fringes = np.divide(proc_data, bg_fit, out=np.zeros_like(proc_data),
                                                                    where=bg_fit != 0)
                                            res_fringes = np.subtract(res_fringes, 1, where=res_fringes != 0)
                                            res_fringes *= np.where(col_weight > 1e-07, 1, 1e-08)

                                            # get the pre-correction contrast using fringe component 1)
                                            if fn == 0:
                                                log.info(" get pre-correction contrast")
                                                pre_contrast = residual_fringe_tools.fit_quality(col_wnum,
                                                                                                 res_fringes,
                                                                                                 weights_feat,
                                                                                                 ffreq[0],
                                                                                                 dffreq[0])
                                                log.info(" pre-correction contrast = {}".format(pre_contrast))

                                            # fit the residual fringes
                                            #pdgm_name = os.path.join(self.plot_dir, out_root + '_' + \
                                            #                         '_periodogram.pdf')

                                            res_fringe_fit, wpix_num = residual_fringe_tools.fit_1d_fringes_bayes_evidence(res_fringes, weights_feat,
                                                                                                     col_wnum, ffreq[fn], dffreq[fn], min_nfringes=min_nfringes[fn],
                                                                                                     max_nfringes=max_nfringes[fn], pgram_res=pgram_res[fn],
                                                                                                     plot_pdgm=self.diagnostic_mode,
                                                                                                     pdgm_name=pdgm_name)

                                            # check for fit blowing up, reset rfc fit to 0, raise a flag
                                            log.info("check residual fringe fit for bad fit regions")
                                            res_fringe_fit, res_fringe_fit_flag = check_res_fringes(res_fringe_fit,
                                                                                                    col_max_amp)

                                            # correct for residual fringes
                                            self.logger.debug(" divide out residual fringe fit, get fringe corrected column")
                                            _, _, _, env, u_x, u_y = fit_envelope(np.arange(res_fringe_fit.shape[0]),
                                                                                  res_fringe_fit)
                                            #rfc_factors = 1 / ((res_fringe_fit-env) * (col_weight > 1e-05).astype(int) + 1)
                                            rfc_factors = 1 / (res_fringe_fit * (col_weight > 1e-05).astype(int) + 1)
                                            proc_data *= rfc_factors
                                            proc_factors *= rfc_factors

                                            # handle nans or infs that may exist
                                            proc_data[proc_data == np.inf] = 0
                                            proc_data = np.nan_to_num(proc_data)

                                # define fringe sub after all fringe components corrections
                                fringe_sub = proc_data.copy()
                                rfc_factors = proc_factors.copy()

                                # get the new fringe contrast
                                self.logger.debug(" analysing fit quality")
                                pbg_fit, pbgindx, pfitter = fit_1d_background_complex(fringe_sub, weights_feat,
                                                                                      col_wnum, ffreq=ffreq[0])

                                # get the residual fringes as fraction of signal
                                fit_res = np.divide(fringe_sub, pbg_fit, out=np.zeros_like(fringe_sub),
                                                    where=pbg_fit != 0)
                                fit_res = np.subtract(fit_res, 1, where=fit_res != 0)
                                fit_res *= np.where(col_weight > 1e-07, 1, 1e-08)

                                if self.diagnostic_mode:
                                    out_name = os.path.join(self.plot_dir,
                                                            out_root + '_fit_quality.pdf')
                                    contrast = fit_quality(col_wnum, fit_res, weights_feat, ffreq, dffreq,
                                                           save_fig=True, plot_name=out_name)
                                else:
                                    contrast = fit_quality(col_wnum, fit_res, weights_feat, ffreq, dffreq)

                                correction_quality.append([contrast, pre_contrast])
                                self.logger.debug(" residual contrast = {}".format(contrast))
                                
                                if self.diagnostic_mode:

                                    # MASK1D plot
                                    out_name = os.path.join(self.plot_dir, out_root + '_mask.pdf')
                                    xdata, ydata = [[col, col]], [[1, 1022]]
                                    diagnostic_plot('mask1d', 'columns', out_name,
                                                    xdata=xdata, ydata=ydata, idata=[all_slice_masks[n]])

                                    # FITS1D plot
                                    out_name = os.path.join(self.plot_dir, out_root + '_res-fringe-fits.pdf')
                                    xdata = [col_wnum]
                                    ydata = [col_data, bg_fit, res_fringes, col_weight, res_fringe_fit, weights_feat,
                                             fringe_sub]
                                    diagnostic_plot('fits1d', 'columns', out_name, xdata=xdata, ydata=ydata,
                                                    save_data=True)

                                # replace the corrected in-slice column pixels in the data_cor array
                                self.logger.debug(" updating the trace pixels in the output")
                                self.output_data[idx, col] = fringe_sub[idx]
                                self.rfc_factors[idx, col] = rfc_factors[idx]
                                self.fit_mask[idx, col] = np.ones(1024)[idx]
                                self.weights_feat[idx, col] = weights_feat[idx]
                                self.weighted_pix_num[idx, col] = np.ones(1024)[idx] * (wpix_num/1024)
                                self.rejected_fit[idx, col] = res_fringe_fit_flag[idx]

                            except Exception as e:
                                self.logger.warning(" Skipping col={}".format(col, ss))
                                self.logger.warning(' %s' % (str(e)))

                                if self.diagnostic_mode:
                                    # save data for testing
                                    self.logger.debug(" saving failed fit data for testing")
                                    out_data = np.array([col_data, col_wnum, col_weight])
                                    out_data_name = self.plot_dir + '/' + out_root + '_data_ERROR_{}'.format(type(e))
                                    np.save(out_data_name, out_data)
                                    del out_data_name, out_data
                                    self.logger.debug(" fit npy file saved")

                del ss_data, ss_wmap, ss_amap, ss_weight  # end on column
                
                # garbage collect after each slice
                gc.collect()

                # asses the fit quality statistics and make plot if in diagnostic mode
                self.logger.debug(" analysing fit statistics")
                if len(correction_quality) > 0:
                    contrasts = np.asarray(correction_quality)[:,0]
                    pre_contrasts = np.asarray(correction_quality)[:, 1]
                    mean, median, stddev, fmax = fit_quality_stats(contrasts)
                    pmean, pmedian, pstddev, pfmax = fit_quality_stats(pre_contrasts)

                    # add the statistics to the Table
                    stat_table.add_row((ss, mean, median, stddev, fmax, pmean, pmedian, pstddev, pfmax))

            del xranges, all_slice_masks, slice_map, slices_in_band, lambda_map, alpha_map

            # garbage collect after each channel
            gc.collect()

        self.logger.info("Processing complete")

        # add units back to output data
        self.logger.debug(" adding units back to output array")
        self.output_data *= np.median(self.input_data)

        # save output to file
        self.logger.debug(" saving the output to file")
        save_data_astropy(self.input, self.input_type, self.output_dir, self.output_data,
                          self.rfc_factors, self.fit_mask, self.weights_feat,
                          self.weighted_pix_num, stat_table, self.rejected_fit)

        self.logger.info("Output saved")
            


        # plot the before and after and the fit statistics histogram
        if self.diagnostic_mode:

            self.logger.debug(" creating fringe residual plot")
            out_name = os.path.join(self.plot_dir, 'fit-stats.pdf')
            fig, axs = plt.subplots(1, 1, figsize=(15, 7))

            x_pos = np.arange(stat_table['Slice'].shape[0]) + 1
            axs.bar(x_pos, stat_table['pmedian'] * 100, align='center', alpha=0.4, color='r', label='before')
            axs.bar(x_pos, stat_table['median'] * 100, align='center', alpha=0.4, color='b', label='after')
            axs.set_xticks(x_pos)
            axs.set_xticklabels(stat_table['Slice'].tolist())
            axs.grid()
            axs.set_ylim(0, 10)
            axs.set_ylabel('Median fringe contrast (%)', fontsize=16)
            axs.set_xlabel('Slice', fontsize=16)
            axs.legend(prop={'size': 10}, loc=0)
            plt.tight_layout()
            fig.savefig(out_name, format='pdf', dpi=150)
            fig.clear()
            plt.close(fig)
            del out_name, fig

        # garbage collect when finished
        gc.collect()

        return None

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

        # replace infs and nans in weights with 0
        weights[weights == np.inf] = 0
        weights[np.isnan(weights)] = 0
        print('value of weights',weights[0:100,473])
        return weights
