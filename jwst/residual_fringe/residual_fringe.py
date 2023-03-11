#
#  Module for applying fringe correction
#

import numpy as np
from functools import partial

from stdatamodels.jwst import datamodels

from ..stpipe import Step
from astropy.table import Table
from astropy.io import ascii
from astropy.io import fits
from . import utils

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class ResidualFringeCorrection():

    def __init__(self,
                 input_model,
                 residual_fringe_reference_file,
                 regions_reference_file,
                 ignore_regions,
                 **pars):

        self.input_model = input_model.copy()
        self.model = input_model.copy()
        self.residual_fringe_reference_file = residual_fringe_reference_file
        self.regions_reference_file = regions_reference_file
        self.ignore_regions = ignore_regions
        self.save_intermediate_results = pars['save_intermediate_results']
        self.transmission_level = int(pars['transmission_level'])
        # define how filenames are created
        self.make_output_path = pars.get('make_output_path',
                                         partial(Step._make_output_path, None))

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

        # used to create additional data that can be plotted outside of step
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

        # Check that the fringe flat has been applied
        if self.input_model.meta.cal_step.fringe != 'COMPLETE':
            raise ErrorNoFringeFlat("The fringe flat step has not been run on file %s",
                                    self.input_model.meta.filename)

        # Remove any NaN values from the data prior to processing
        # Set them to 0 for the residual fringe routine
        # They will be re-added at the end
        output_data = self.model.data.copy()
        nanval_indx = np.where(~np.isfinite(output_data))
        output_data[nanval_indx] = 0

        # normalise the output_data to remove units
        pos_data = self.input_model.data[self.input_model.data > 0]
        normalization_factor = np.median(pos_data)
        output_data /= normalization_factor

        # Load the fringe reference file
        residual_fringe_model = datamodels.FringeFreqModel(self.residual_fringe_reference_file)

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
        self.transmission_level = int(self.transmission_level / 10)

        slice_map = (allregions.regions)[self.transmission_level - 1, :, :].copy()
        log.info(" Using {} throughput threshold.".format(self.transmission_level))

        self.slice_map = slice_map

        # set up the channels for the detector
        detector = self.input_model.meta.instrument.detector.lower()

        if 'short' in detector:
            self.channels = [1, 2]
        elif 'long' in detector:
            self.channels = [3, 4]
        log.info("Detector {} {} ".format(detector, self.channels))

        self.input_weights = self.calc_weights()
        self.weights_feat = self.input_weights.copy()

        self.rfc_factors = np.zeros(self.input_model.data.shape)
        self.fit_mask = np.zeros(self.input_model.data.shape)
        self.weighted_pix_num = np.zeros(self.input_model.data.shape)
        self.rejected_fit = np.zeros(self.input_model.data.shape)

        self.background_fit = np.zeros(self.input_model.data.shape)
        self.knot_locations = np.full_like(self.input_model.data, np.NaN)
        allregions.close()

        # intermediate output product - Table
        stat_table = Table(names=('Slice', 'mean', 'median', 'stddev', 'max', 'pmean', 'pmedian', 'pstddev', 'pmax'),
                           dtype=('i4', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8'))

        out_table = Table(names=('Slice', 'col', 'fringe', 'sn', 'pre_contrast', 'post_contrast',
                                 'periodogram_res', 'opt_fringes', 'peak_freq', 'freq_min', 'freq_max'),
                          dtype=('i4', 'i4', 'i4', 'f8', 'f8', 'f8',
                                 'f8', 'f8', 'f8', 'f8', 'f8'))

        # intermediate output
        n_wav_samples = 1024
        qual_table = Table(names=('column', 'quality'),
                           units=('', ''),
                           dtype=('i4', '(3,{})f8'.format(n_wav_samples)))

        ysize = self.input_model.data.shape[0]
        xsize = self.input_model.data.shape[1]

        for c in self.channels:
            num_corrected = 0
            log.info("Processing channel {}".format(c))
            (slices_in_band, xrange_channel, slice_x_ranges, all_slice_masks) = \
                utils.slice_info(slice_map, c)

            log.debug(' Slice Ranges {}'.format(slice_x_ranges))

            y, x = np.mgrid[:ysize, :xsize]
            _, _, wave_map = self.input_model.meta.wcs(x, y)

            # if the user wants to ignore some values use the wave_map array to set the corresponding
            # weight values to 0
            if self.ignore_regions['num'] > 0:
                for r in range(self.ignore_regions['num']):
                    min_wave = self.ignore_regions['min'][r]
                    max_wave = self.ignore_regions['max'][r]
                    self.input_weights[((wave_map > min_wave) & (wave_map < max_wave))] = 0

            for n, ss in enumerate(slices_in_band):
                log.info(" Processing slice {} =================================".format(ss))
                log.debug(" X ranges of slice {} {}".format(slice_x_ranges[n, 1], slice_x_ranges[n, 2]))
                # initialise the list to store correction quality numbers for slice
                correction_quality = []

                # use the mask to set all out-of-slice pixels to 0 in wmap and data
                # set out-of-slice pixels to 0 in arrays
                ss_data = all_slice_masks[n] * output_data.copy()
                ss_wmap = all_slice_masks[n] * wave_map
                ss_weight = all_slice_masks[n] * self.input_weights.copy()
                #  ss_mask = all_slice_masks[n]

                # get the freq_table info for this slice
                this_row = np.where(self.freq_table['slice'] == float(ss))[0][0]
                log.debug('Row in reference file for slice {}'.format(this_row))

                slice_row = self.freq_table[(self.freq_table['slice'] == float(ss))]
                ffreq = slice_row['ffreq'][0]
                dffreq = slice_row['dffreq'][0]
                min_nfringes = slice_row['min_nfringes'][0]
                max_nfringes = slice_row['max_nfringes'][0]
                min_snr = slice_row['min_snr'][0]
                pgram_res = slice_row['pgram_res'][0]

                # cycle through the cols and fit the fringes
                for col in np.arange(slice_x_ranges[n, 1], slice_x_ranges[n, 2]):
                    col_data = ss_data[:, col]
                    col_wmap = ss_wmap[:, col]

                    valid = np.logical_and((col_wmap > 0), ~np.isnan(col_wmap))
                    # because of the curvature of the slices there can be large regions not falling on a column
                    num_good = len(np.where(valid)[0])
                    # Need at least 50 pixels in column to proceed

                    if num_good > 50:
                        test_flux = col_data[valid]
                        test_flux[test_flux < 0] = 1e-08
                        # Transform wavelength in micron to wavenumber in cm^-1.
                        col_wnum = 10000.0 / col_wmap

                        # use the error array to get col snr, used to remove noisey pixels
                        col_snr = self.model.data.copy()[:, col] / self.model.err.copy()[:, col]

                        # do some checks on column to make sure there is reasonable signal. If the SNR < min_snr (CDP), pass
                        # determine SNR for this column of data
                        n = len(test_flux)
                        signal = np.nanmean(test_flux)
                        noise = 0.6052697 * np.nanmedian(np.abs(2.0 * test_flux[2:n - 2] - test_flux[0:n - 4] - test_flux[4:n]))

                        snr2 = 0.0  # initialize
                        if noise != 0:
                            snr2 = signal / noise

                        # Sometimes can return nan, inf for bad data so include this in check
                        if snr2 < min_snr[0]:
                            log.debug('SNR too high not fitting column {}, {}, {}'.format(col, snr2, min_snr[0]))
                            pass
                        else:
                            log.debug("Fitting column{} ".format(col))
                            log.debug("SNR > {} ".format(min_snr[0]))

                            col_weight = ss_weight[:, col]
                            col_max_amp = np.interp(col_wmap, self.max_amp['Wavelength'], self.max_amp['Amplitude'])
                            col_snr2 = np.where(col_snr > 10, 1, 0)  # hardcoded at snr > 10 for now

                            # get the in-slice pixel indices for replacing in output later
                            idx = np.where(col_data > 0)

                            # BayesicFitting doesn't like 0s at data or weight array edges so set to small value
                            # replacing array 0s with arbitrarily low number
                            col_data[col_data <= 0] = 1e-08
                            col_weight[col_weight <= 0] = 1e-08

                            # check for off-slice pixels and send to be filled with interpolated/extrapolated wnums
                            # to stop BayesicFitting crashing, will not be fitted anyway
                            # finding out-of-slice pixels in column and filling

                            found_bad = np.logical_or(np.isnan(col_wnum), np.isinf(col_wnum))
                            num_bad = len(np.where(found_bad)[0])

                            if num_bad > 0:
                                col_wnum[found_bad] = 0
                                col_wnum = utils.fill_wavenumbers(col_wnum)

                            # do feature finding on slice now column-by-column
                            log.debug(" starting feature finding")

                            # narrow features (similar or less than fringe #1 period)
                            # find spectral features (env is spline fit of troughs and peaks)
                            env, l_x, l_y, _, _, _ = utils.fit_envelope(np.arange(col_data.shape[0]), col_data)
                            mod = np.abs(col_data / env) - 1

                            # given signal in mod find location of lines > col_max_amp * 2 (fringe contrast)
                            # use col_snr to ignore noisy pixels
                            weight_factors = utils.find_lines(mod * col_snr2, col_max_amp * 2)
                            weights_feat = col_weight * weight_factors

                            # account for fringe 2 on broad features in channels 3 and 4
                            # need to smooth out the dichroic fringe as it breaks the feature finding method
                            if c in [3, 4]:
                                win = 7  # smooting window hardcoded to 7 for now (based on testing)
                                cumsum = np.cumsum(np.insert(col_data, 0, 0))
                                sm_col_data = (cumsum[win:] - cumsum[:-win]) / float(win)

                                # find spectral features (env is spline fit of troughs and peaks)
                                env, l_x, l_y, _, _, _ = utils.fit_envelope(np.arange(col_data.shape[0]), sm_col_data)
                                mod = np.abs(col_data / env) - 1

                                # given signal in mod find location of lines > col_max_amp * 2
                                weight_factors = utils.find_lines(mod, col_max_amp * 2)
                                weights_feat *= weight_factors

                            # iterate over the fringe components to fit, initialize pre-contrast, other output arrays
                            # in case fit fails
                            proc_data = col_data.copy()
                            proc_factors = np.ones(col_data.shape)
                            pre_contrast = 0.0
                            bg_fit = col_data.copy()
                            res_fringes = np.zeros(col_data.shape)
                            res_fringe_fit = np.zeros(col_data.shape)
                            res_fringe_fit_flag = np.zeros(col_data.shape)
                            wpix_num = 1024

                            # check the end points. A single value followed by gap of zero can cause
                            # problems in the fitting.
                            index = np.where(weights_feat != 0.0)
                            length = np.diff(index[0])

                            if weights_feat[0] != 0 and length[0] > 1:
                                weights_feat[0] = 1e-08

                            if weights_feat[-1] != 0 and length[-1] > 1:
                                weights_feat[-1] = 1e-08

                            # jane addded this - fit can fail in evidence function.
                            # once we replace evidence funtion with astropy routine - we can test
                            # removing setting weights < 0.003 to zero (1e-08)

                            weights_feat[weights_feat <= 0.003] = 1e-08

                            # currently the reference file fits one fringe originating in the detector pixels, and a
                            # second high frequency, low amplitude fringe in channels 3 and 4 which has been
                            # attributed to the dichroics.
                            try:
                                for fn, ff in enumerate(ffreq):
                                    # ignore place holder fringes
                                    if ff > 1e-03:
                                        log.debug(' starting ffreq = {}'.format(ff))

                                        # check if snr criteria is met for fringe component, should always be true for fringe 1
                                        if snr2 > min_snr[fn]:
                                            log.debug(" fitting spectral baseline")

                                            bg_fit, bgindx = \
                                                utils.fit_1d_background_complex(proc_data, weights_feat,
                                                                                col_wnum, ffreq=ffreq[fn], channel=c)

                                            # get the residual fringes as fraction of signal
                                            res_fringes = np.divide(proc_data, bg_fit, out=np.zeros_like(proc_data),
                                                                    where=bg_fit != 0)
                                            res_fringes = np.subtract(res_fringes, 1, where=res_fringes != 0)
                                            res_fringes *= np.where(col_weight > 1e-07, 1, 1e-08)
                                            # get the pre-correction contrast using fringe component 1
                                            # TODO: REMOVE CONTRAST CHECKS
                                            # set dummy values for contrast check parameters until removed
                                            pre_contrast = 0.0
                                            quality = np.array([np.zeros(col_data.shape), np.zeros(col_data.shape),
                                                                np.zeros(col_data.shape)])
                                            # if fn == 0:
                                            #    pre_contrast, quality = utils.fit_quality(col_wnum,
                                            #                                              res_fringes,
                                            #                                              weights_feat,
                                            #                                              ffreq[0],
                                            #                                              dffreq[0])
                                            #
                                            #   log.debug(" pre-correction contrast = {}".format(pre_contrast))
                                            #
                                            # fit the residual fringes
                                            log.debug(" set up bayes ")
                                            res_fringe_fit, wpix_num, opt_nfringe, peak_freq, freq_min, freq_max = \
                                                utils.new_fit_1d_fringes_bayes_evidence(res_fringes,
                                                                                        weights_feat,
                                                                                        col_wnum,
                                                                                        ffreq[fn],
                                                                                        dffreq[fn],
                                                                                        min_nfringes=min_nfringes[fn],
                                                                                        max_nfringes=max_nfringes[fn],
                                                                                        pgram_res=pgram_res[fn],
                                                                                        col_snr2=col_snr2)

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
                                            proc_data = np.nan_to_num(proc_data, posinf=1e-08, neginf=1e-08)
                                            proc_data[proc_data < 0] = 1e-08

                                            out_table.add_row((ss, col, fn, snr2, pre_contrast, pre_contrast, pgram_res[fn],
                                                               opt_nfringe, peak_freq, freq_min, freq_max))

                                # define fringe sub after all fringe components corrections
                                fringe_sub = proc_data.copy()
                                rfc_factors = proc_factors.copy()

                                # get the new fringe contrast
                                log.debug(" analysing fit quality")

                                pbg_fit, pbgindx = utils.fit_1d_background_complex(fringe_sub,
                                                                                   weights_feat,
                                                                                   col_wnum,
                                                                                   ffreq=ffreq[0], channel=c)

                                # get the residual fringes as fraction of signal
                                fit_res = np.divide(fringe_sub, pbg_fit, out=np.zeros_like(fringe_sub),
                                                    where=pbg_fit != 0)
                                fit_res = np.subtract(fit_res, 1, where=fit_res != 0)
                                fit_res *= np.where(col_weight > 1e-07, 1, 1e-08)

                                # TODO: REMOVE CONTRAST CHECKS
                                # set dummy values for contrast check parameters until removed
                                contrast = 0.0
                                quality = np.array([np.zeros(col_data.shape), np.zeros(col_data.shape), np.zeros(col_data.shape)])
                                # contrast, quality = utils.fit_quality(col_wnum,
                                #                                      fit_res,
                                #                                      weights_feat,
                                #                                      ffreq,
                                #                                      dffreq,
                                #                                      save_results=self.save_intermediate_results)

                                out_table.add_row((ss, col, fn, snr2, pre_contrast, contrast, pgram_res[0],
                                                   opt_nfringe, peak_freq, freq_min, freq_max))

                                qual_table.add_row((col, quality))
                                correction_quality.append([contrast, pre_contrast])
                                log.debug(" residual contrast = {}".format(contrast))

                                # replace the corrected in-slice column pixels in the data_cor array
                                log.debug(" updating the trace pixels in the output")
                                output_data[idx, col] = fringe_sub[idx]
                                self.rfc_factors[idx, col] = rfc_factors[idx]
                                self.fit_mask[idx, col] = np.ones(1024)[idx]
                                self.weights_feat[idx, col] = weights_feat[idx]
                                self.weighted_pix_num[idx, col] = np.ones(1024)[idx] * (wpix_num / 1024)
                                self.rejected_fit[idx, col] = res_fringe_fit_flag[idx]
                                self.background_fit[idx, col] = bg_fit[idx]
                                self.knot_locations[:bgindx.shape[0], col] = bgindx
                                num_corrected = num_corrected + 1

                            except Exception as e:
                                log.warning(" Skipping col={} {} ".format(col, ss))
                                log.warning(' %s' % (str(e)))

                del ss_data, ss_wmap, ss_weight  # end on column

                # asses the fit quality statistics and set up data to make plot outside of step
                log.debug(" analysing fit statistics")
                if len(correction_quality) > 0:
                    contrasts = np.asarray(correction_quality)[:, 0]
                    pre_contrasts = np.asarray(correction_quality)[:, 1]
                    mean, median, stddev, fmax = utils.fit_quality_stats(contrasts)
                    pmean, pmedian, pstddev, pfmax = utils.fit_quality_stats(pre_contrasts)

                    # add the statistics to the Table
                    stat_table.add_row((ss, mean, median, stddev, fmax, pmean, pmedian, pstddev, pfmax))

            del slice_x_ranges, all_slice_masks, slices_in_band, wave_map  # end of channel
            log.info('Number of columns corrected for channel {}'.format(num_corrected))
        log.info("Processing complete")

        # add units back to output data
        log.info(" adding units back to output array")
        output_data *= normalization_factor
        # Add NaNs back to output data
        output_data[nanval_indx] = np.nan
        self.model.data = output_data
        del output_data

        if self.save_intermediate_results:
            stat_table_name = self.make_output_path(
                basepath=self.input_model.meta.filename,
                suffix='stat_table', ext='.ecsv')
            log.info(' Saving intermediate Stat table {}'.format(stat_table_name))
            ascii.write(stat_table, stat_table_name, format='ecsv', fast_writer=False, overwrite=True)

            out_table_name = self.make_output_path(
                basepath=self.input_model.meta.filename,
                suffix='out_table', ext='.ecsv')
            log.info(' Saving intermediate Output table {}'.format(out_table_name))
            ascii.write(out_table, out_table_name, format='ecsv', fast_writer=False, overwrite=True)

            t = fits.BinTableHDU(data=qual_table, name='FIT_QUAL')

            fit_results_name = self.make_output_path(
                basepath=self.input_model.meta.filename,
                suffix='fit_results', ext='.fits')
            log.info('Saving intermediate fit results output {}'.format(fit_results_name))
            h = fits.open(self.input_model.meta.filename)
            hdr = h[0].header
            h.close()

            hdu0 = fits.PrimaryHDU(header=hdr)
            hdu1 = fits.ImageHDU(self.rfc_factors, name='RFC_FACTORS')
            hdu2 = fits.ImageHDU(self.fit_mask, name='FIT_MASK')
            hdu3 = fits.ImageHDU(self.weights_feat, name='WEIGHTS_FEATURES')
            hdu4 = fits.ImageHDU(self.weighted_pix_num, name='WEIGHTED_PIXEL_FRACTION')
            hdu5 = fits.ImageHDU(self.background_fit, name='BACKGROUND_FIT')
            hdu6 = fits.ImageHDU(self.knot_locations, name='KNOT_LOCATIONS')

            hdu = fits.HDUList([hdu0, hdu1, hdu2, hdu3, hdu4, hdu5, hdu6, t])
            hdu.writeto(fit_results_name, overwrite=True)
            hdu.close()

        return self.model

    def calc_weights(self):

        """Make a weights array based on flux. This is a placeholder function,
        for now just returns a normalised flux array to use as weights array.
        This is because any smoothing results in incorrect rfc around emission lines.
        This can be changed in the future if need be.

        :Parameters:

        flux:  numpy array, required
        the 1D array of fluxes

        :Returns:

        weights array which is just a copy of the flux array, normalised by the mean

        """

        weights = np.zeros(self.input_model.data.shape)
        for c in np.arange(weights.shape[1]):
            flux_1d = self.input_model.data[:, c]
            w = flux_1d / np.nanmean(flux_1d)
            weights[:, c] = w

        # replace infs and nans in weights with 0
        weights[weights == np.inf] = 0
        weights[np.isnan(weights)] = 0
        return weights


class ErrorNoFringeFlat(Exception):
    pass
