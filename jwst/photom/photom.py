import logging
import numpy as np
from astropy import units as u

from .. import datamodels
from .. datamodels import dqflags
from .. master_background import expand_to_2d

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

PHOT_TOL = 0.001  # relative tolerance between PIXAR_* keys

# Conversion factor from MJy/sr to uJy/arcsec^2
MJSR_TO_UJA2 = (u.megajansky/u.steradian).to(u.microjansky/u.arcsecond/u.arcsecond)


class DataSet():
    """
    Input dataset to which the photom information will be applied

    Parameters
    ----------

    """
    def __init__(self, model):
        """
        Short Summary
        -------------
        Store vital params in DataSet object, such as
        instrument, detector, filter, pupil, and exposure type.

        Parameters
        ----------
        model : `~jwst.datamodels.DataModel`
            input Data Model object

        """
        # Create a copy of the input model
        self.input = model.copy()

        self.instrument = model.meta.instrument.name.upper()
        self.detector = model.meta.instrument.detector.upper()
        self.exptype = model.meta.exposure.type.upper()
        self.filter = None
        if model.meta.instrument.filter is not None:
            self.filter = model.meta.instrument.filter.upper()
        self.pupil = None
        if model.meta.instrument.pupil is not None:
            self.pupil = model.meta.instrument.pupil.upper()
        self.grating = None
        if model.meta.instrument.grating is not None:
            self.grating = model.meta.instrument.grating.upper()
        self.band = None
        if model.meta.instrument.band is not None:
            self.band = model.meta.instrument.band.upper()
        self.slitnum = -1

        # Let the user know what we're working with
        log.info('Using instrument: %s', self.instrument)
        log.info(' detector: %s', self.detector)
        log.info(' exp_type: %s', self.exptype)
        if self.filter is not None:
            log.info(' filter: %s', self.filter)
        if self.pupil is not None:
            log.info(' pupil: %s', self.pupil)
        if self.grating is not None:
            log.info(' grating: %s', self.grating)
        if self.band is not None:
            log.info(' band: %s', self.band)

    def calc_nirspec(self, ftab, area_fname):
        """
        Extended Summary
        -------------
        For the NIRSPEC instrument, reference file matching is based on
        FILTER and GRATING, as well as SLIT name for the fixed-slits mode.
        The routine will find the corresponding information in the reference
        file, apply it to the data, and write the scalar conversion
        factor to the output model. All NIRSpec modes use wavelength-dependent
        flux calibration factors.

        Parameters
        ----------
        ftab : `~jwst.datamodels.NirspecPhotomModel` or `~jwst.datamodels.NirspecFSPhotomModel`
            NIRSpec photom reference file data model

        area_fname : str
            Pixel area map reference file name

        Returns
        -------

        """

        # Get the GRATING value from the input data model
        grating = self.input.meta.instrument.grating.upper()

        # Normal fixed-slit exposures get handled as a MultiSlitModel
        if self.exptype == 'NRS_FIXEDSLIT':

            # We have to find and apply a separate set of flux cal
            # data for each of the fixed slits in the input
            for slit in self.input.slits:

                log.info('Working on slit %s' % slit.name)
                match = False
                self.slitnum += 1

                # Loop through reference table to find matching row
                for tabdata in ftab.phot_table:
                    ref_filter = tabdata['filter'].strip().upper()
                    ref_grating = tabdata['grating'].strip().upper()
                    ref_slit = tabdata['slit'].strip().upper()
                    log.debug(' Ref table data: %s %s %s' %
                              (ref_filter, ref_grating, ref_slit))

                    # Match on filter, grating, and slit name
                    if (self.filter == ref_filter and
                            grating == ref_grating and slit.name == ref_slit):
                        self.photom_io(tabdata)
                        match = True
                        break

                if not match:
                    log.warning('No match in reference file')

        # Bright object fixed-slit exposures use a CubeModel
        elif self.exptype == 'NRS_BRIGHTOBJ':

            match = False

            # Bright object always uses S1600A1 slit
            slit_name = 'S1600A1'
            log.info('Working on slit %s' % slit_name)

            # Loop through reference table to find matching row
            for tabdata in ftab.phot_table:
                ref_filter = tabdata['filter'].strip().upper()
                ref_grating = tabdata['grating'].strip().upper()
                ref_slit = tabdata['slit'].strip().upper()
                log.debug(' Ref table data: %s %s %s' %
                          (ref_filter, ref_grating, ref_slit))

                # Match on filter, grating, and slit name
                if (self.filter == ref_filter and
                        grating == ref_grating and slit_name == ref_slit):
                    self.photom_io(tabdata)
                    match = True
                    break

            if not match:
                log.warning('No match in reference file')

        # IFU and MSA exposures use one set of flux cal data
        else:

            match = False

            # Loop through the reference table to find matching row
            for tabdata in ftab.phot_table:

                ref_filter = tabdata['filter'].strip().upper()
                ref_grating = tabdata['grating'].strip().upper()

                # Match on filter and grating only
                if self.filter == ref_filter and grating == ref_grating:

                    match = True

                    # MSA data
                    if (isinstance(self.input, datamodels.MultiSlitModel) and
                            self.exptype == 'NRS_MSASPEC'):

                        # Loop over the MSA slits, applying the same photom
                        # ref data to all slits
                        for slit in self.input.slits:
                            log.info('Working on slit %s' % slit.name)
                            self.slitnum += 1
                            self.photom_io(tabdata)

                    # IFU data
                    else:

                        # Get the conversion factor from the PHOTMJSR column
                        conv_factor = tabdata['photmjsr']

                        # Populate the photometry keywords
                        self.input.meta.photometry.conversion_megajanskys = \
                            conv_factor
                        self.input.meta.photometry.conversion_microjanskys = \
                            conv_factor * MJSR_TO_UJA2

                        # Get the length of the relative response arrays in
                        # this table row
                        nelem = tabdata['nelem']

                        # If the relative response arrays have length > 0,
                        # load them for use in creating a 2-d array of
                        # flux conversion factors
                        if nelem > 0:
                            waves = tabdata['wavelength'][:nelem]
                            relresps = tabdata['relresponse'][:nelem]

                            # Convert wavelengths from meters to microns,
                            # if necessary
                            microns_100 = 1.e-4    # 100 microns, in meters
                            if waves.max() > 0. and waves.max() < microns_100:
                                waves *= 1.e+6

                        # Load the pixel area table for the IFU slices
                        area_model = datamodels.NirspecIfuAreaModel(area_fname)
                        area_data = area_model.area_table

                        # Compute 2D wavelength and pixel area arrays for the
                        # whole image
                        wave2d, area2d, dqmap = self.calc_nrs_ifu_sens2d(area_data)

                        # Compute relative sensitivity for each pixel based
                        # on its wavelength
                        sens2d = np.interp(wave2d, waves, relresps)

                        # Include the scalar conversion factor
                        sens2d *= conv_factor

                        # Divide by pixel area
                        sens2d /= area2d

                        # Reset NON_SCIENCE pixels to 1 in sens2d array and flag
                        # them in the science data DQ array
                        where_dq = \
                            np.bitwise_and(dqmap, dqflags.pixel['NON_SCIENCE'])
                        sens2d[where_dq > 0] = 1.
                        self.input.dq = np.bitwise_or(self.input.dq, dqmap)

                        # Divide the science data and uncertainty arrays by the
                        # conversion factors
                        self.input.data /= sens2d
                        self.input.err /= sens2d
                        self.input.var_poisson /= sens2d**2
                        self.input.var_rnoise /= sens2d**2

                        # Update BUNIT values for the science data and err
                        self.input.meta.bunit_data = 'mJy/arcsec^2'
                        self.input.meta.bunit_err = 'mJy/arcsec^2'

                        area_model.close()

                    break

            if not match:
                log.warning('No match in reference file')

        return

    def calc_niriss(self, ftab):
        """
        Extended Summary
        -------------
        For NIRISS matching is based on FILTER and PUPIL, as well as ORDER
        for spectroscopic modes.
        There may be multiple entries for a given FILTER+PUPIL combination,
        corresponding to different spectral orders. Data for all orders will
        be retrieved.

        The routine will find the corresponding information in the reference
        file and apply the conversions to the science arrays.
        If wavelength-dependent conversion factors exist, which will be the
        case for spectroscopic modes, they will be loaded and applied along
        with the scalar conversion.

        Parameters
        ----------
        ftab : `~jwst.datamodels.NirissPhotomModel`
            NIRISS photom reference file data model

        Returns
        -------
        """

        # Handle MultiSlit models separately, which are used for NIRISS WFSS
        if isinstance(self.input, datamodels.MultiSlitModel):

            # We have to find and apply a separate set of flux cal
            # data for each of the slits/orders in the input
            for slit in self.input.slits:

                # Initialize the output conversion factor and increment slit number
                match = False
                self.slitnum += 1

                # Get the spectral order number for this slit
                order = slit.meta.wcsinfo.spectral_order

                log.info("Working on slit: {} order: {}".format(slit.name, order))

                # Locate matching row in reference file
                for tabdata in ftab.phot_table:

                    ref_filter = tabdata['filter'].strip().upper()
                    ref_pupil = tabdata['pupil'].strip().upper()
                    ref_order = tabdata['order']

                    # Find matching values of FILTER, PUPIL, ORDER
                    if (self.filter == ref_filter and self.pupil == ref_pupil and order == ref_order):
                        self.photom_io(tabdata)
                        match = True
                        break

                if not match:
                    log.warning('No match in reference file')

        # NIRISS imaging and SOSS modes
        else:

            # Hardwire the science data order number to 1 for now
            order = 1
            match = False

            # Locate matching row in reference file
            for tabdata in ftab.phot_table:
                ref_filter = tabdata['filter'].strip().upper()
                ref_pupil = tabdata['pupil'].strip().upper()
                ref_order = tabdata['order']

                # SOSS mode
                if self.exptype in ['NIS_SOSS']:

                    # Find matching values of FILTER, PUPIL, and ORDER
                    if (self.filter == ref_filter and self.pupil == ref_pupil and order == ref_order):
                        self.photom_io(tabdata, order)
                        match = True
                        break

                # Imaging mode
                else:

                    # Find matching values of FILTER and PUPIL
                    if (self.filter == ref_filter and self.pupil == ref_pupil):
                        self.photom_io(tabdata)
                        match = True
                        break

            if not match:
                log.warning('No match in reference file')

        return

    def calc_miri(self, ftab):
        """
        Extended Summary
        -------------
        For MIRI imaging and LRS modes, matching is based on FILTER and SUBARRAY.
        MIRI MRS uses dedicated photom reference files per CHANNEL+BAND.

        For Imaging and LRS, the routine will find the corresponding row of
        information in the reference file, apply it, and store the scalar
        conversion factor in the output model PHOTMJSR keyword. If
        wavelength-dependent conversion values exist, which is the case for LRS
        mode, they will be included in the applied conversion.

        Parameters
        ----------
        ftab : `~jwst.datamodels.MiriImgPhotomModel` or `~jwst.datamodels.MiriMrsPhotomModel`
            MIRI photom reference file data model

        Returns
        -------
        """

        match = False

        # Imaging detector
        if self.detector == 'MIRIMAGE':

            # Get the subarray value of the input data model
            subarray = self.input.meta.subarray.name
            log.info(' subarray: %s', subarray)

            # Find the matching row in the reference file
            for tabdata in ftab.phot_table:

                ref_filter = tabdata['filter'].strip().upper()
                ref_subarray = tabdata['subarray'].strip().upper()

                # If the ref file subarray entry is GENERIC, it's an automatic
                # match to the science data subarray value
                if ref_subarray == 'GENERIC':
                    ref_subarray = subarray

                # Find matching FILTER and SUBARRAY values
                if self.filter == ref_filter and subarray == ref_subarray:
                    self.photom_io(tabdata)
                    match = True
                    break

            if not match:
                log.warning('No match in reference file')

        # MRS detectors
        elif self.detector == 'MIRIFUSHORT' or self.detector == 'MIRIFULONG':

            # Reset conversion and pixel size values with DQ=NON_SCIENCE to 1,
            # so no conversion is applied
            where_dq = np.bitwise_and(ftab.dq, dqflags.pixel['NON_SCIENCE'])
            ftab.data[where_dq > 0] = 1.0
            ftab.pixsiz[where_dq > 0] = 1.0

            # Reset NaN's in conversion array to 1
            where_nan = np.isnan(ftab.data)
            ftab.data[where_nan] = 1.0

            # Reset zeros in pixsiz array to 1
            where_zero = np.where(ftab.pixsiz == 0.0)
            ftab.pixsiz[where_zero] = 1.0

            # Make sure all NaN's and zeros have DQ flags set
            ftab.dq[where_nan] = np.bitwise_or(ftab.dq[where_nan],
                                               dqflags.pixel['NON_SCIENCE'])
            ftab.dq[where_zero] = np.bitwise_or(ftab.dq[where_zero],
                                                dqflags.pixel['NON_SCIENCE'])

            # Compute the combined 2D sensitivity factors
            sens2d = ftab.data * ftab.pixsiz

            # Divide the science data and uncertainty arrays by the 2D sensitivity factors
            self.input.data /= sens2d
            self.input.err /= sens2d
            self.input.var_poisson /= sens2d**2
            self.input.var_rnoise /= sens2d**2

            # Update the science dq
            self.input.dq = np.bitwise_or(self.input.dq, ftab.dq)

            # Retrieve the scalar conversion factor from the reference data
            conv_factor = ftab.meta.photometry.conversion_megajanskys

            # Store the conversion factors in the meta data
            self.input.meta.photometry.conversion_megajanskys = \
                conv_factor
            self.input.meta.photometry.conversion_microjanskys = \
                conv_factor * MJSR_TO_UJA2

            # Update BUNIT values for the science data and err
            self.input.meta.bunit_data = 'mJy/arcsec^2'
            self.input.meta.bunit_err = 'mJy/arcsec^2'

        return

    def calc_nircam(self, ftab):
        """
        Extended Summary
        -------------
        For NIRCAM, matching is based on FILTER and PUPIL.
        The routine will find the corresponding information in the reference
        file, apply the conversion factors, and store the scalar conversion
        factor in the output model. If wavelength-dependent conversion factors
        exist, they will be included in the calibration.
        For WFSS (grism) mode, the calibration information extracted from the
        reference file is applied to each slit instance in the science data.

        Parameters
        ----------
        ftab : `~jwst.datamodels.NircamPhotomModel`
            NIRCam photom reference file data model

        Returns
        -------
        """

        match = False

        # Locate relevant information in reference file
        for tabdata in ftab.phot_table:

            ref_filter = tabdata['filter'].strip().upper()
            ref_pupil = tabdata['pupil'].strip().upper()

            # Finding matching FILTER and PUPIL values
            if self.filter == ref_filter and self.pupil == ref_pupil:
                match = True

                # Handle WFSS data separately from regular imaging
                if (isinstance(self.input, datamodels.MultiSlitModel) and self.exptype == 'NRC_WFSS'):

                    # Loop over the WFSS slits, applying the same photom
                    # ref data to all slits
                    for slit in self.input.slits:
                        log.info('Working on slit %s' % slit.name)
                        self.slitnum += 1
                        self.photom_io(tabdata)

                else:

                    # Regular imaging data only requires 1 call
                    self.photom_io(tabdata)

                break

        if not match:
            log.warning('No match in reference file')

        return

    def calc_fgs(self, ftab):
        """
        Extended Summary
        -------------
        For FGS, there is no matching required, because the instrument does
        not contain any filters or pupil wheel elements. The only mode is CLEAR.

        The routine will find the corresponding information in the reference
        file (which should have only a single row), apply it to the data, and
        write the conversion factor to the output model.

        Parameters
        ----------
        ftab : `~jwst.datamodels.FgsPhotomModel`
            FGS photom reference file data model

        Returns
        -------
        """

        # Read the first (and only) row in the reference file
        for tabdata in ftab.phot_table:
            self.photom_io(tabdata)
            break

        return

    def calc_nrs_ifu_sens2d(self, area_data):
        """Create the 2-D wavelength and pixel area arrays needed for
        constructing a NIRSpec IFU sensitivity map.

        Parameters
        ----------
        area_data : 1-D ndarray
            Array of pixel area values for the IFU slices

        Returns
        -------
        wave2d : 2-D ndarray
            Array of wavelengths per pixel

        area2d : 2-D ndarray
            Array of pixel area values

        dqmap : 2-D ndarray
            Array of DQ flags per pixel
        """

        import numpy as np
        from .. assign_wcs import nirspec       # for NIRSpec IFU data
        import gwcs

        microns_100 = 1.e-4                     # 100 microns, in meters

        # Create empty 2D arrays for the wavelengths and pixel areas
        wave2d = np.zeros_like(self.input.data)
        area2d = np.ones_like(self.input.data)

        # Create and initialize an array for the 2D dq map to be returned.
        # initialize all pixels to NON_SCIENCE, because operations below
        # only touch pixels within the bounding_box of each slice
        dqmap = np.zeros_like(self.input.dq) + dqflags.pixel['NON_SCIENCE']

        # Get the list of wcs's for the IFU slices
        list_of_wcs = nirspec.nrs_ifu_wcs(self.input)

        # Loop over the slices
        for (k, ifu_wcs) in enumerate(list_of_wcs):

            # Construct array indexes for pixels in this slice
            x, y = gwcs.wcstools.grid_from_bounding_box(ifu_wcs.bounding_box,
                                                        step=(1, 1),
                                                        center=True)

            log.debug("Slice %d: %g %g %g %g" %
                      (k, x[0][0], x[-1][-1], y[0][0], y[-1][-1]))

            # Get the world coords for all pixels in this slice
            coords = ifu_wcs(x, y)

            # Pull out the wavelengths only
            wl = coords[2]
            nan_flag = np.isnan(wl)
            good_flag = np.logical_not(nan_flag)
            if wl[good_flag].max() < microns_100:
                log.info("Wavelengths in WCS table appear to be in meters")

            # Set NaNs to a harmless value, but don't modify nan_flag.
            wl[nan_flag] = 0.

            # Mark pixels with no wavelength as non-science
            dq = np.zeros_like(wl)
            dq[nan_flag] = dqflags.pixel['NON_SCIENCE']
            dqmap[y.astype(int), x.astype(int)] = dq

            # Insert the wavelength values for this slice into the
            # whole image array
            wave2d[y.astype(int), x.astype(int)] = wl

            # Insert the pixel area value for this slice into the
            # whole image array
            ar = np.ones_like(wl)
            ar[:, :] = area_data[np.where(area_data['slice_id'] == k)]['pixarea'][0]
            ar[nan_flag] = 1.
            area2d[y.astype(int), x.astype(int)] = ar

        return wave2d, area2d, dqmap

    def photom_io(self, tabdata, order=None):
        """
        Short Summary
        -------------
        Combine photometric scalar and wavelength-dependent conversion factors
        and apply to the science dataset.

        Parameters
        ----------
        tabdata : FITS record
            Single row of data from reference table

        order : int
            Spectral order number

        Returns
        -------

        """
        # Get the scalar conversion factor from the PHOTMJSR column of the table row
        conversion = tabdata['photmjsr']

        # Store the conversion factor in the meta data
        log.info('PHOTMJSR value: %g', conversion)
        if isinstance(self.input, datamodels.MultiSlitModel):
            self.input.slits[self.slitnum].meta.photometry.conversion_megajanskys = \
                conversion
            self.input.slits[self.slitnum].meta.photometry.conversion_microjanskys = \
                conversion * MJSR_TO_UJA2
        else:
            self.input.meta.photometry.conversion_megajanskys = conversion
            self.input.meta.photometry.conversion_microjanskys = conversion * MJSR_TO_UJA2

        # Get the length of the relative response arrays in this row
        nelem = tabdata['nelem']

        # If the relative response arrays have length > 0, load and include them in
        # the flux conversion
        if nelem > 0:
            waves = tabdata['wavelength'][:nelem]
            relresps = tabdata['relresponse'][:nelem]

            # Make sure waves and relresps are in increasing wavelength order
            if not np.all(np.diff(waves) > 0):
                index = np.argsort(waves)
                waves = waves[index].copy()
                relresps = relresps[index].copy()

            # Convert wavelengths from meters to microns, if necessary
            microns_100 = 1.e-4         # 100 microns, in meters
            if waves.max() > 0. and waves.max() < microns_100:
                waves *= 1.e+6

            # Compute a 2-D grid of conversion factors, as a function of wavelength
            if isinstance(self.input, datamodels.MultiSlitModel):
                wl_array = expand_to_2d.get_wavelengths(
                                self.input.slits[self.slitnum],
                                self.input.meta.exposure.type,
                                order)
            else:
                wl_array = expand_to_2d.get_wavelengths(
                                self.input,
                                self.input.meta.exposure.type,
                                order)

            wl_array[np.isnan(wl_array)] = -1.
            conv_2d = np.interp(wl_array, waves, relresps, left=1., right=1.)

            # Combine the scalar and 2-D conversions
            # NOTE: the 2-D conversion is divided into the data for now, until the
            # instrument teams deliver multiplicative conversions in photom ref files
            conversion = conversion / conv_2d

        # Apply the conversion to the data and all uncertainty arrays
        if isinstance(self.input, datamodels.MultiSlitModel):
            self.input.slits[self.slitnum].data *= conversion
            self.input.slits[self.slitnum].err *= conversion
            if (self.input.slits[self.slitnum].var_poisson is not None and
                np.size(self.input.slits[self.slitnum].var_poisson) > 0):
                    self.input.slits[self.slitnum].var_poisson *= conversion**2
            if (self.input.slits[self.slitnum].var_rnoise is not None and
                np.size(self.input.slits[self.slitnum].var_rnoise) > 0):
                    self.input.slits[self.slitnum].var_rnoise *= conversion**2
            self.input.slits[self.slitnum].meta.bunit_data = 'MJy/sr'
            self.input.slits[self.slitnum].meta.bunit_err = 'MJy/sr'
        else:
            self.input.data *= conversion
            self.input.err *= conversion
            if (self.input.var_poisson is not None and
                np.size(self.input.var_poisson) > 0):
                    self.input.var_poisson *= conversion**2
            if (self.input.var_rnoise is not None and
                np.size(self.input.var_rnoise) > 0):
                    self.input.var_rnoise *= conversion**2
            self.input.meta.bunit_data = 'MJy/sr'
            self.input.meta.bunit_err = 'MJy/sr'

        return

    def save_area_info(self, ftab, area_fname):
        """
        Short Summary
        -------------
        Read the pixel area values in the PIXAR_A2 and PIXAR_SR keys from the
        meta data in the photom reference file and the pixel area reference
        file. Copy the values from the pixel area reference file header
        keywords to the output product. If the difference between the values
        of the pixel area (in units of arc seconds) between the two reference
        files exceeds a defined threshold, issue a warning.

        Also copy the pixel area data array from the pixel area reference file
        to the area extension of the output product.

        Parameters
        ----------
        ftab : `~jwst.datamodels.DataModel`
            Photom reference file data model

        area_fname : str
            Pixel area reference file name

        Returns
        -------

        """

        # Load the pixel area reference file
        pix_area = datamodels.PixelAreaModel(area_fname)

        # Copy the pixel area data array to the appropriate attribute
        # of the science data model
        if isinstance(self.input, datamodels.MultiSlitModel):
            self.input.slits[0].area = pix_area.data
        else:
            self.input.area = pix_area.data
        log.info('Pixel area map copied to output.')

        # Load the average pixel area values from the photom reference file
        try:
            #tab_ster = None
            tab_a2 = None
            #tab_ster = ftab.meta.photometry.pixelarea_steradians
            tab_a2 = ftab.meta.photometry.pixelarea_arcsecsq
        except AttributeError:
            # If one or both of them are missing, issue a warning, but carry on
            log.warning('At least one of the PIXAR_nn keyword values is')
            log.warning('missing from the photom reference file')

        # Load the average pixel area values from the pixel area reference file
        try:
            area_ster = None
            area_a2 = None
            area_ster = pix_area.meta.photometry.pixelarea_steradians
            area_a2 = pix_area.meta.photometry.pixelarea_arcsecsq
        except AttributeError:
            # If one or both of them are missing, issue a warning
            log.warning('At least one of the PIXAR_nn keyword values is')
            log.warning('missing from the reference file %s', area_fname)
            log.warning('Pixel area keyword values will not be set in output')

        # Compute the relative difference between the pixel area values from
        # the two different sources, if they exist
        if (tab_a2 is not None) and (area_a2 is not None):
            a2_tol = abs(tab_a2 - area_a2) / (tab_a2 + area_a2)

            # If the difference is greater than the defined tolerance,
            # issue a warning
            if (a2_tol > PHOT_TOL):
                log.warning('The relative difference between the values for')
                log.warning('the pixel area in sq arcsec (%s)', a2_tol)
                log.warning('exceeds the defined tolerance (%s)', PHOT_TOL)

        # Copy the pixel area values to the output
        log.debug('The values of the pixel areas (PIXAR_A2 and PIXAR_SR)')
        log.debug('will be copied to the output.')
        if area_a2 is not None:
            self.input.meta.photometry.pixelarea_arcsecsq = float(area_a2)
        if area_ster is not None:
            self.input.meta.photometry.pixelarea_steradians = float(area_ster)

        return

    def apply_photom(self, photom_fname, area_fname):
        """
        Short Summary
        -------------
        Open the reference file, retrieve the conversion factors from the reference
        file that are appropriate to the instrument mode. This can consist of both
        a scalar factor and an array of wavelength-dependent factors. The combination
        of both factors are applied to the input model. The scalar factor is also
        written to the PHOTMJSR and PHOTUJA2 keywords in the model. If a pixel area
        map reference file exists for the instrument mode, it is attached to
        the input model.

        Parameters
        ----------
        photom_fname : str
            photom reference file name

        area_fname: str
            pixel area map reference file name

        Returns
        -------
        output_model : ~jwst.datamodels.DataModel
            output data model with the flux calibrations applied

        """

        # Load the photom reference file into the appropriate type of datamodel
        # for the instrument mode in use and then call the calculation routine
        # for that instrument mode
        if self.instrument == 'NIRISS':
            ftab = datamodels.NirissPhotomModel(photom_fname)
            self.calc_niriss(ftab)

        if self.instrument == 'NIRSPEC':
            if self.exptype in ['NRS_FIXEDSLIT', 'NRS_BRIGHTOBJ']:
                ftab = datamodels.NirspecFSPhotomModel(photom_fname)
            else:
                ftab = datamodels.NirspecPhotomModel(photom_fname)
            self.calc_nirspec(ftab, area_fname)

        if self.instrument == 'NIRCAM':
            ftab = datamodels.NircamPhotomModel(photom_fname)
            self.calc_nircam(ftab)

        if self.instrument == 'MIRI':
            if self.detector == 'MIRIMAGE':
                ftab = datamodels.MiriImgPhotomModel(photom_fname)
            else:
                ftab = datamodels.MiriMrsPhotomModel(photom_fname)
            self.calc_miri(ftab)

        if self.instrument == 'FGS':
            ftab = datamodels.FgsPhotomModel(photom_fname)
            self.calc_fgs(ftab)

        # Load the pixel area reference file, if it exists, and attach the
        # the reference data to the science model
        if area_fname:
            if 'IMAGE' in self.exptype and area_fname != 'N/A':
                self.save_area_info(ftab, area_fname)

        ftab.close()  # Close the photom reference table

        return self.input
