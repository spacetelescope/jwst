import logging
import numpy as np
from astropy import units as u

from .. import datamodels
from .. datamodels import dqflags
from .. lib.wcs_utils import get_wavelengths

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

PHOT_TOL = 0.001  # relative tolerance between PIXAR_* keys

# Conversion factor from MJy/sr to uJy/arcsec^2
MJSR_TO_UJA2 = (u.megajansky/u.steradian).to(u.microjansky/u.arcsecond/u.arcsecond)

# Conversion factor from square arcseconds to steradians
A2_TO_SR = (np.pi / (180. * 3600.))**2


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
        ftab : `~jwst.datamodels.NrsFsPhotomModel` or `~jwst.datamodels.NrsMosPhotomModel`
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
                        # this table row.  If the nelem column is not present,
                        # we'll use the entire wavelength and relresponse
                        # arrays.
                        try:
                            nelem = tabdata['nelem']
                        except KeyError:
                            nelem = None

                        waves = tabdata['wavelength']
                        relresps = tabdata['relresponse']
                        if nelem is not None:
                            waves = waves[:nelem]
                            relresps = relresps[:nelem]

                        # Convert wavelengths from meters to microns,
                        # if necessary
                        microns_100 = 1.e-4    # 100 microns, in meters
                        if waves.max() > 0. and waves.max() < microns_100:
                            waves *= 1.e+6

                        # Load the pixel area table for the IFU slices
                        area_model = datamodels.open(area_fname)
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

                        # Multiply the science data and uncertainty arrays by
                        # the conversion factors
                        self.input.data *= sens2d
                        self.input.err *= sens2d
                        self.input.var_poisson *= sens2d**2
                        self.input.var_rnoise *= sens2d**2
                        if self.input.var_flat is not None and np.size(self.input.var_flat) > 0:
                            self.input.var_flat *= sens2d**2

                        # Update BUNIT values for the science data and err
                        self.input.meta.bunit_data = 'MJy/sr'
                        self.input.meta.bunit_err = 'MJy/sr'

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
        ftab : `~jwst.datamodels.NisSossPhotomModel` or
               `~jwst.datamodels.NisWfssPhotomModel` or
               `~jwst.datamodels.NisImgPhotomModel`
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
                try:
                    ref_order = tabdata['order']
                except KeyError:
                    ref_order = order

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
        ftab : `~jwst.datamodels.MirImgPhotomModel` or
               `~jwst.datamodels.MirMrsPhotomModel` or
               `~jwst.datamodels.MirLrsPhotomModel`
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

            # Reset NaN's in conversion array to 1
            where_nan = np.isnan(ftab.data)
            ftab.data[where_nan] = 1.0

            # Make sure all NaN's and zeros have DQ flags set
            ftab.dq[where_nan] = np.bitwise_or(ftab.dq[where_nan],
                                               dqflags.pixel['NON_SCIENCE'])

            # Compute the combined 2D sensitivity factors
            sens2d = ftab.data

            # Multiply the science data and uncertainty arrays by the 2D
            # sensitivity factors
            self.input.data *= sens2d
            self.input.err *= sens2d
            self.input.var_poisson *= sens2d**2
            self.input.var_rnoise *= sens2d**2
            if self.input.var_flat is not None and np.size(self.input.var_flat) > 0:
                self.input.var_flat *= sens2d**2

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
            self.input.meta.bunit_data = 'MJy/sr'
            self.input.meta.bunit_err = 'MJy/sr'

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
        ftab : `~jwst.datamodels.NrcImgPhotomModel` or
               `~jwst.datamodels.NrcWfssPhotomModel`
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
        ftab : `~jwst.datamodels.FgsImgPhotomModel`
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
        # First get the scalar conversion factor.
        # For most modes, the scalar conversion factor in the photom reference
        # file is in units of (MJy / sr) / (DN / s), and the output from
        # the photom step will be in units of surface brightness, specifically
        # MJy / sr.  For NIRSpec and NIRISS SOSS, however, the scalar
        # conversion factor is in units of MJy / (DN / s); for point-source
        # targets, the output will be in units of flux density, MJy.  For an
        # extended source (or if the source type is unknown), the conversion
        # factor will be divided by the solid angle of a pixel, so the output
        # of the photom step will be in units of surface brightness, as for
        # other types of data.
        unit_is_surface_brightness = True               # default
        try:
            conversion = tabdata['photmjsr']            # unit is MJy / sr
        except KeyError:
            conversion = tabdata['photmj']              # unit is MJy
            if isinstance(self.input, datamodels.MultiSlitModel):
                slit = self.input.slits[self.slitnum]
                if self.input.meta.exposure.type == 'NRS_MSASPEC':
                    srctype = slit.source_type
                else:
                    srctype = self.input.meta.target.source_type
                if srctype is None or srctype.upper() != 'POINT':
                    log.debug('Converting conversion factor from flux '
                              'to surface brightness')
                    conversion /= slit.meta.photometry.pixelarea_steradians
                else:
                    unit_is_surface_brightness = False
            else:
                srctype = self.input.meta.target.source_type
                if srctype is None or srctype.upper() != 'POINT':
                    log.debug('Converting conversion factor from flux '
                              'to surface brightness')
                    conversion /= self.input.meta.photometry.pixelarea_steradians
                else:
                    unit_is_surface_brightness = False

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

        # If the photom reference file is for spectroscopic data, the table
        # in the reference file should contain a 'wavelength' column (among
        # other columns).
        try:
            wl_test = tabdata['wavelength']
            is_spectroscopic = True
            del wl_test
        except KeyError:
            is_spectroscopic = False

        # Get the length of the relative response arrays in this row.  If the
        # nelem column is not present, we'll use the entire wavelength and
        # relresponse arrays.
        if is_spectroscopic:
            try:
                nelem = tabdata['nelem']
            except KeyError:
                nelem = None
        else:
            nelem = None

        # For spectroscopic data, include the relative response array in
        # the flux conversion.
        no_cal = None
        if is_spectroscopic:
            waves = tabdata['wavelength']
            relresps = tabdata['relresponse']
            if nelem is not None:
                waves = waves[:nelem]
                relresps = relresps[:nelem]

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
                wl_array = get_wavelengths(self.input.slits[self.slitnum],
                                           self.input.meta.exposure.type,
                                           order)
            else:
                wl_array = get_wavelengths(self.input,
                                           self.input.meta.exposure.type,
                                           order)

            wl_array[np.isnan(wl_array)] = -1.
            conv_2d = np.interp(wl_array, waves, relresps, left=np.NaN, right=np.NaN)
            no_cal = np.isnan(conv_2d)

            # Combine the scalar and 2-D conversions
            # NOTE: the data are now multiplied by the 2-D conversion
            conversion = conversion * conv_2d
            conversion[no_cal] = 0.

        # Apply the conversion to the data and all uncertainty arrays
        if isinstance(self.input, datamodels.MultiSlitModel):
            slit = self.input.slits[self.slitnum]
            slit.data *= conversion
            slit.err *= conversion
            if slit.var_poisson is not None and np.size(slit.var_poisson) > 0:
                slit.var_poisson *= conversion**2
            if slit.var_rnoise is not None and np.size(slit.var_rnoise) > 0:
                slit.var_rnoise *= conversion**2
            if slit.var_flat is not None and np.size(slit.var_flat) > 0:
                slit.var_flat *= conversion**2
            if no_cal is not None:
                slit.dq[..., no_cal] = np.bitwise_or(slit.dq[..., no_cal],
                                                     dqflags.pixel['DO_NOT_USE'])
            if unit_is_surface_brightness:
                slit.meta.bunit_data = 'MJy/sr'
                slit.meta.bunit_err = 'MJy/sr'
            else:
                slit.meta.bunit_data = 'MJy'
                slit.meta.bunit_err = 'MJy'
        else:
            self.input.data *= conversion
            self.input.err *= conversion
            if self.input.var_poisson is not None and np.size(self.input.var_poisson) > 0:
                self.input.var_poisson *= conversion**2
            if self.input.var_rnoise is not None and np.size(self.input.var_rnoise) > 0:
                self.input.var_rnoise *= conversion**2
            if self.input.var_flat is not None and np.size(self.input.var_flat) > 0:
                self.input.var_flat *= conversion**2
            if no_cal is not None:
                self.input.dq[..., no_cal] = np.bitwise_or(self.input.dq[..., no_cal],
                                                           dqflags.pixel['DO_NOT_USE'])
            if unit_is_surface_brightness:
                self.input.meta.bunit_data = 'MJy/sr'
                self.input.meta.bunit_err = 'MJy/sr'
            else:
                self.input.meta.bunit_data = 'MJy'
                self.input.meta.bunit_err = 'MJy'

        return

    def save_area_info(self, area_fname):
        """
        Short Summary
        -------------
        Read the pixel area values in the PIXAR_A2 and PIXAR_SR keys from the
        primary header of the pixel area reference file or (for NIRSpec data)
        from the PIXAREA column in the selected row of the AREA table in the
        area reference file.  Use that information to populate the pixel area
        keywords in the output product.

        Except for NIRSpec data, also copy the pixel area data array from the
        pixel area reference file to the area extension of the output product.

        Parameters
        ----------
        area_fname : str
            Pixel area reference file name
        """

        # We need the instrument name in order to check for NIRSpec.
        instrument = self.input.meta.instrument.name.upper()

        # Load the pixel area reference file
        pix_area = datamodels.open(area_fname)

        # Copy the pixel area data array to the appropriate attribute
        # of the science data model
        if instrument != 'NIRSPEC':
            if isinstance(self.input, datamodels.MultiSlitModel):
                # Note that this only copied to the first slit.
                self.input.slits[0].area = pix_area.data
            else:
                self.input.area = pix_area.data
            log.info('Pixel area map copied to output.')

        if instrument == 'NIRSPEC':
            self.save_area_nirspec(pix_area)
        else:
            # Load the average pixel area values from the pixel area reference
            try:
                area_ster = pix_area.meta.photometry.pixelarea_steradians
            except AttributeError:
                area_ster = None
                log.warning('The PIXAR_SR keyword is missing from %s',
                            area_fname)
            try:
                area_a2 = pix_area.meta.photometry.pixelarea_arcsecsq
            except AttributeError:
                area_a2 = None
                log.warning('The PIXAR_A2 keyword is missing from %s',
                            area_fname)

            # Copy the pixel area values to the output
            log.debug('PIXAR_SR = %s, PIXAR_A2 = %s', str(area_ster), str(area_a2))
            if area_a2 is not None:
                self.input.meta.photometry.pixelarea_arcsecsq = float(area_a2)
            if area_ster is not None:
                self.input.meta.photometry.pixelarea_steradians = float(area_ster)

        pix_area.close()

    def save_area_nirspec(self, pix_area):
        """
        Short Summary
        -------------
        Read the pixel area value from the PIXAREA column in the selected row
        of the AREA table in the area reference file.
        Use that information to populate the pixel area keywords in the output
        product.

        Parameters
        ----------
        pix_area : `~jwst.datamodels.DataModel`
            Pixel area reference file data model
        """

        exp_type = self.input.meta.exposure.type
        pixarea = pix_area.area_table['pixarea']
        if exp_type == 'NRS_MSASPEC':
            quadrant = pix_area.area_table['quadrant']
            shutter_x = pix_area.area_table['shutter_x']
            shutter_y = pix_area.area_table['shutter_y']
            n_failures = 0
            for slit in self.input.slits:
                match_q = (quadrant == slit.quadrant)
                match_x = (shutter_x == slit.xcen)
                match_y = (shutter_y == slit.ycen)
                match = np.logical_and(match_q,
                                       np.logical_and(match_x, match_y))
                n_matches = match.sum(dtype=np.int64)
                if n_matches != 1:
                    n_failures += 1
                    slit.meta.photometry.pixelarea_arcsecsq = 1.
                    slit.meta.photometry.pixelarea_steradians = 1.
                else:
                    slit.meta.photometry.pixelarea_arcsecsq = float(pixarea[match])
                    slit.meta.photometry.pixelarea_steradians = \
                        slit.meta.photometry.pixelarea_arcsecsq * A2_TO_SR
            if n_failures > 0:
                log.warning('%d failures out of %d, matching MSA data to '
                            'area reference file',
                            n_failures, len(self.input.slits))

        elif exp_type in ['NRS_LAMP', 'NRS_BRIGHTOBJ', 'NRS_FIXEDSLIT']:
            slit_id = pix_area.area_table['slit_id']
            foundit = False
            for k, slit in enumerate(self.input.slits):
                if slit_id[k] == slit.name:
                    foundit = True
                    slit.meta.photometry.pixelarea_arcsecsq = float(pixarea[k])
                    slit.meta.photometry.pixelarea_steradians = \
                        slit.meta.photometry.pixelarea_arcsecsq * A2_TO_SR
            if not foundit:
                slit.meta.photometry.pixelarea_arcsecsq = 1.
                slit.meta.photometry.pixelarea_steradians = 1.

        elif exp_type == 'NRS_IFU':
            # There is a slice_id column for selecting a matching slice, but
            # we're just going to average the pixel area for all slices.
            pixel_area = np.nanmean(pixarea)
            slit.meta.photometry.pixelarea_arcsecsq = pixel_area
            slit.meta.photometry.pixelarea_steradians = \
                slit.meta.photometry.pixelarea_arcsecsq * A2_TO_SR

        else:
            log.warning('EXP_TYPE of NIRSpec data is %s, which is not an '
                        'expected value; pixel area keywords will be set to 1.',
                        exp_type)
            self.input.meta.photometry.pixelarea_arcsecsq = 1.
            self.input.meta.photometry.pixelarea_steradians = 1.

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

        ftab = datamodels.open(photom_fname)

        # Load the pixel area reference file, if it exists, and attach the
        # reference data to the science model
        if area_fname != 'N/A':
            self.save_area_info(area_fname)

        if self.instrument == 'NIRISS':
            self.calc_niriss(ftab)

        elif self.instrument == 'NIRSPEC':
            self.calc_nirspec(ftab, area_fname)

        elif self.instrument == 'NIRCAM':
            self.calc_nircam(ftab)

        elif self.instrument == 'MIRI':
            self.calc_miri(ftab)

        elif self.instrument == 'FGS':
            self.calc_fgs(ftab)

        else:
            raise RuntimeError('Instrument {} is not recognized'
                               .format(self.instrument))

        ftab.close()  # Close the photom reference table

        return self.input
