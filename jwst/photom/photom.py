import logging
import functools
import warnings

import numpy as np
from astropy import units as u

from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels import dqflags

from jwst.lib.pipe_utils import match_nans_and_flags
from jwst.lib.wcs_utils import get_wavelengths
from jwst.lib.dispaxis import get_dispersion_direction
from . import miri_mrs
from . import miri_imager

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

PHOT_TOL = 0.001  # relative tolerance between PIXAR_* keys

# Conversion factor from MJy/sr to uJy/arcsec^2
MJSR_TO_UJA2 = (u.megajansky / u.steradian).to(u.microjansky / u.arcsecond / u.arcsecond)

# Conversion factor from square arcseconds to steradians
A2_TO_SR = (np.pi / (180.0 * 3600.0)) ** 2


class MatchFitsTableRowError(Exception):
    """Class to handle error when matching FITS Table rows."""

    def __init__(self, message):
        if message is None:
            message = "Expected to match one row in FITS table."
        super().__init__(message)


class DataModelTypeError(Exception):
    """Class to handle case of unexpected datamodel type."""

    def __init__(self, message):
        if message is None:
            message = "Unexpected DataModel type."
        super().__init__(message)


def find_row(fits_table, match_fields):
    """
    Find a row in a FITS table matching fields.

    Parameters
    ----------
    fits_table : `~astropy.io.fits.fitsrec.FITS_rec`
        FITS table
    match_fields : dict
        {field_name: value} pair to use as a matching criteria.

    Returns
    -------
    row : int, or None
        FITS table row index, None if no match.

    Raises
    ------
    Warning
        When a field name is not in the table.
    MatchFitsTableRowError
        When more than one row matches.
    """

    def _normalize_strings(field):
        if isinstance(field[0], str):
            return np.array([s.upper() for s in field])
        return field

    # item[1] is always converted to upper case in the `DataSet` initializer.
    results = [
        _normalize_strings(fits_table.field(item[0])) == item[1] for item in match_fields.items()
    ]
    row = functools.reduce(np.logical_and, results).nonzero()[0]
    if len(row) > 1:
        raise MatchFitsTableRowError(
            f"Expected to find one matching row in table, found {len(row)}."
        )
    if len(row) == 0:
        log.warning("Expected to find one matching row in table, found 0.")
        return None
    return row[0]


class DataSet:
    """
    Input dataset to which the photom information will be applied.

    Store vital params, such as
    instrument, detector, filter, pupil, and exposure type.
    """

    def __init__(
        self,
        model,
        inverse=False,
        source_type=None,
        mrs_time_correction=False,
        correction_pars=None,
    ):
        """
        Instantiate a DataSet object.

        Parameters
        ----------
        model : `~jwst.datamodels.JwstDataModel`
            Input Data Model object.
        inverse : bool
            Invert the math operations used to apply the corrections.
        source_type : str or None
            Force processing using the specified source type.
        mrs_time_correction : bool
            Switch to apply/not apply the MRS time correction.
        correction_pars : dict
            Correction meta-data from a previous run.
        """
        # Set up attributes necessary for calculation.
        if correction_pars:
            self.update(correction_pars["dataset"])
        else:
            self.band = None
            if model.meta.instrument.band is not None:
                self.band = model.meta.instrument.band.upper()
            self.instrument = model.meta.instrument.name.upper()
            self.detector = model.meta.instrument.detector.upper()
            self.exptype = model.meta.exposure.type.upper()
            self.filter = None
            if model.meta.instrument.filter is not None:
                self.filter = model.meta.instrument.filter.upper()
            self.grating = None
            if model.meta.instrument.grating is not None:
                self.grating = model.meta.instrument.grating.upper()
            self.order = None
            if (
                model.meta.hasattr("wcsinfo")
                and model.meta.wcsinfo.hasattr("spectral_order")
                and model.meta.wcsinfo.spectral_order is not None
            ):
                self.order = model.meta.wcsinfo.spectral_order
            self.pupil = None
            if model.meta.instrument.pupil is not None:
                self.pupil = model.meta.instrument.pupil.upper()
            self.subarray = None
            if model.meta.subarray.name is not None:
                self.subarray = model.meta.subarray.name.upper()
            correction_pars = {}
        correction_pars["dataset"] = self.attributes
        self.correction_pars = correction_pars

        # Initialize other non-correction pars attributes.
        self.slitnum = -1
        self.specnum = -1
        self.integ_row = -1
        self.inverse = inverse
        self.source_type = None
        self.mrs_time_correction = mrs_time_correction

        # For MultiSlitModels, only set a generic source_type value for the
        # entire datamodel if the user has set the source_type parameter.
        # Otherwise leave the generic source_type set to None, which will
        # force the use of per-slit source_types later in processing.
        if isinstance(model, datamodels.MultiSlitModel):
            if source_type is not None:
                self.source_type = source_type

        # For non-MultiSlitModel inputs, where there's only 1 target and
        # one source_type, use the user-provided source_type value, if it
        # exists. Otherwise, use the generic source_type value provided in
        # the input model (if it exists).
        else:
            if source_type is not None:
                self.source_type = source_type
            else:
                if model.meta.target.source_type is not None:
                    self.source_type = model.meta.target.source_type.upper()

        # Create a copy of the input model
        self.input = model.copy()

        # Let the user know what we're working with
        log.info("Using instrument: %s", self.instrument)
        log.info(" detector: %s", self.detector)
        log.info(" exp_type: %s", self.exptype)
        if self.filter is not None:
            log.info(" filter: %s", self.filter)
        if self.pupil is not None:
            log.info(" pupil: %s", self.pupil)
        if self.grating is not None:
            log.info(" grating: %s", self.grating)
        if self.band is not None:
            log.info(" band: %s", self.band)

    @property
    def attributes(self):
        """
        Retrieve DataSet attributes.

        Returns
        -------
        attributes : dict
            A dict of `DataSet` attributes.
        """
        attributes = vars(self)

        # Remove some attributes
        for attribute in ["correction_pars", "input", "inverse", "slitnum", "source_type"]:
            if attribute in attributes:
                del attributes[attribute]

        return attributes

    def update(self, attributes):
        """
        Set DataSet attributes.

        Parameters
        ----------
        attributes : dict
            The attributes to be set on DataSet.
        """
        for key, value in attributes.items():
            setattr(self, key, value)

    def calc_nirspec(self, ftab, area_fname):
        """
        Apply photometric calibration data to dataset and update conversion factor.

        For the NIRSPEC instrument, reference file matching is based on
        FILTER and GRATING, as well as SLIT name for the fixed-slits mode.
        The routine will find the corresponding information in the reference
        file, apply it to the data, and write the scalar conversion
        factor to the output model. All NIRSpec modes use wavelength-dependent
        flux calibration factors.

        Parameters
        ----------
        ftab : `~jwst.datamodels.NrsFsPhotomModel` or `~jwst.datamodels.NrsMosPhotomModel`
            NIRSpec photom reference file data model.
        area_fname : str
            Pixel area map reference file name.
        """
        # Normal fixed-slit exposures get handled as a MultiSlitModel
        if self.exptype == "NRS_FIXEDSLIT":
            # We have to find and apply a separate set of flux cal
            # data for each of the fixed slits in the input
            for slit in self.input.slits:
                log.info(f"Working on slit {slit.name}")
                self.slitnum += 1

                fields_to_match = {
                    "filter": self.filter,
                    "grating": self.grating,
                    "slit": slit.name,
                }
                row = find_row(ftab.phot_table, fields_to_match)
                if row is None:
                    continue
                self.photom_io(ftab.phot_table[row])

        # Bright object fixed-slit exposures use a SlitModel
        elif self.exptype == "NRS_BRIGHTOBJ":
            # Bright object always uses S1600A1 slit
            slit_name = "S1600A1"
            log.info(f"Working on slit {slit_name}")
            fields_to_match = {"filter": self.filter, "grating": self.grating, "slit": slit_name}
            row = find_row(ftab.phot_table, fields_to_match)
            if row is None:
                return
            self.photom_io(ftab.phot_table[row])

        # IFU and MSA exposures use one set of flux cal data
        else:
            fields_to_match = {"filter": self.filter, "grating": self.grating}
            row = find_row(ftab.phot_table, fields_to_match)
            if row is None:
                return

            # MSA (MOS) data
            if isinstance(self.input, datamodels.MultiSlitModel) and self.exptype == "NRS_MSASPEC":
                # Loop over the MSA slits, applying the same photom
                # ref data to all slits
                for slit in self.input.slits:
                    log.info(f"Working on slit {slit.name}")
                    self.slitnum += 1
                    self.photom_io(ftab.phot_table[row])

            # IFU data
            else:
                tabdata = ftab.phot_table[row]

                # Get the scalar conversion factor from the PHOTMJ column;
                conv_factor = tabdata["photmj"]
                unit_is_surface_brightness = False

                # Convert to surface brightness units for all types of data
                if self.input.meta.photometry.pixelarea_steradians is None:
                    log.warning("Pixel area is None, so can't convert flux to surface brightness!")
                else:
                    log.debug("Converting conversion factor from flux to surface brightness")
                    conv_factor /= self.input.meta.photometry.pixelarea_steradians
                    unit_is_surface_brightness = True

                # Populate the photometry keywords
                log.info(f"PHOTMJSR value: {conv_factor:.6g}")
                self.input.meta.photometry.conversion_megajanskys = conv_factor
                self.input.meta.photometry.conversion_microjanskys = conv_factor * MJSR_TO_UJA2

                # Get the length of the relative response arrays in this table row.
                # If the nelem column is not present, we'll use the entire wavelength
                # and relresponse arrays.
                try:
                    nelem = tabdata["nelem"]
                except KeyError:
                    nelem = None

                waves = tabdata["wavelength"]
                relresps = tabdata["relresponse"]
                if nelem is not None:
                    waves = waves[:nelem]
                    relresps = relresps[:nelem]

                # Convert wavelengths from meters to microns, if necessary
                microns_100 = 1.0e-4  # 100 microns, in meters
                if waves.max() > 0.0 and waves.max() < microns_100:
                    waves *= 1.0e6

                # Load the pixel area table for the IFU slices
                area_model = datamodels.open(area_fname)
                area_data = area_model.area_table

                # Compute 2D wavelength and pixel area arrays for the whole image
                wave2d, area2d, dqmap = self.calc_nrs_ifu_sens2d(area_data)

                # Compute relative sensitivity for each pixel based
                # on its wavelength
                sens2d = np.interp(wave2d, waves, relresps)
                sens2d *= tabdata["photmj"]  # include the initial scalar conversion factor -> MJ
                # convert all data (both point source and extended) to surface brightness
                sens2d /= area2d * A2_TO_SR  # divide by pixel area * A2_TO_SR -> MJy/sr

                # Reset NON_SCIENCE pixels to 1 in sens2d array and flag
                # them in the science data DQ array
                where_dq = np.bitwise_and(dqmap, dqflags.pixel["NON_SCIENCE"])
                sens2d[where_dq > 0] = 1.0
                self.input.dq = np.bitwise_or(self.input.dq, dqmap)

                # Multiply the science data and uncertainty arrays by the conversion factors
                sens2d_squared = sens2d * sens2d
                if not self.inverse:
                    self.input.data *= sens2d
                    self.input.err *= sens2d
                    self.input.var_poisson *= sens2d_squared
                    self.input.var_rnoise *= sens2d_squared
                else:
                    self.input.data /= sens2d
                    self.input.err /= sens2d
                    self.input.var_poisson /= sens2d_squared
                    self.input.var_rnoise /= sens2d_squared

                if self.input.var_flat is not None and np.size(self.input.var_flat) > 0:
                    if not self.inverse:
                        self.input.var_flat *= sens2d_squared
                    else:
                        self.input.var_flat /= sens2d_squared

                # Update BUNIT values for the science data and err
                if not self.inverse:
                    if unit_is_surface_brightness:
                        self.input.meta.bunit_data = "MJy/sr"
                        self.input.meta.bunit_err = "MJy/sr"
                    else:
                        self.input.meta.bunit_data = "MJy"
                        self.input.meta.bunit_err = "MJy"
                else:
                    self.input.meta.bunit_data = "DN/s"
                    self.input.meta.bunit_err = "DN/s"

                area_model.close()

    def calc_niriss(self, ftab):
        """
        Apply photometric calibration data to dataset and update conversion factor.

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
        ftab : `~jwst.datamodels.NisSossPhotomModel` or `~jwst.datamodels.NisWfssPhotomModel` or
               `~jwst.datamodels.NisImgPhotomModel`
            NIRISS photom reference file data model.
        """
        # Handle MultiSlit models separately, which are used for NIRISS WFSS
        if isinstance(self.input, datamodels.MultiSlitModel):
            # We have to find and apply a separate set of flux cal
            # data for each of the slits/orders in the input
            for slit in self.input.slits:
                # Increment slit number
                self.slitnum += 1

                # Get the spectral order number for this slit
                order = slit.meta.wcsinfo.spectral_order
                log.info(f"Working on slit {slit.name}, order {order}")

                fields_to_match = {"filter": self.filter, "pupil": self.pupil, "order": order}
                row = find_row(ftab.phot_table, fields_to_match)
                if row is None:
                    continue
                self.photom_io(ftab.phot_table[row])

        elif isinstance(self.input, datamodels.CubeModel):
            raise DataModelTypeError(
                f"Unexpected input data model type for NIRISS: {str(self.input)}"
            )

        elif self.exptype in ["NIS_SOSS"]:
            if isinstance(self.input, datamodels.ImageModel):
                raise DataModelTypeError(
                    f"Unexpected input data model type for NIRISS: {str(self.input)}"
                )

            for spec in self.input.spec:
                self.specnum += 1
                self.order = spec.spectral_order
                fields_to_match = {"filter": self.filter, "pupil": self.pupil, "order": self.order}
                row = find_row(ftab.phot_table, fields_to_match)
                if row is None:
                    return
                # Correct each integration
                for integ_row in range(len(spec.spec_table)):
                    self.integ_row = integ_row
                    self.photom_io(ftab.phot_table[row], self.order)
        else:
            fields_to_match = {"filter": self.filter, "pupil": self.pupil}
            row = find_row(ftab.phot_table, fields_to_match)
            if row is None:
                return
            self.photom_io(ftab.phot_table[row])

    def calc_miri(self, ftab):
        """
        Apply photometric calibration data to dataset and update conversion factor.

        For MIRI imaging and LRS modes, matching is based on FILTER and SUBARRAY.
        MIRI MRS uses dedicated photom reference files per CHANNEL+BAND.

        For Imaging and LRS, the routine will find the corresponding row of
        information in the reference file, apply it, and store the scalar
        conversion factor in the output model PHOTMJSR keyword. If
        wavelength-dependent conversion values exist, which is the case for LRS
        mode, they will be included in the applied conversion.

        Parameters
        ----------
        ftab : `~jwst.datamodels.MirImgPhotomModel` or `~jwst.datamodels.MirMrsPhotomModel`
               or `~jwst.datamodels.MirLrsPhotomModel`
            MIRI photom reference file data model.
        """
        # Imaging detector
        if self.detector == "MIRIMAGE":
            # Get the subarray value of the input data model
            log.info(" subarray: %s", self.subarray)
            fields_to_match = {"subarray": self.subarray, "filter": self.filter}
            row = find_row(ftab.phot_table, fields_to_match)
            if row is None:
                # Search again using subarray="GENERIC" for old ref files
                fields_to_match = {"subarray": "GENERIC", "filter": self.filter}
                row = find_row(ftab.phot_table, fields_to_match)
                if row is None:
                    return

            # Check to see if the reference file contains the coefficients for the
            # time-dependent correction of the PHOTOM value
            try:
                ftab.getarray_noinit("timecoeff")
                log.info("Applying the time-dependent correction to the PHOTOM value.")

                mid_time = self.input.meta.exposure.mid_time
                photom_corr = miri_imager.time_corr_photom(ftab.timecoeff[row], mid_time)

                data = np.array(
                    [
                        (
                            self.filter,
                            self.subarray,
                            ftab.phot_table[row]["photmjsr"] + photom_corr,
                            ftab.phot_table[row]["uncertainty"],
                        )
                    ],
                    dtype=[
                        ("filter", "O"),
                        ("subarray", "O"),
                        ("photmjsr", "<f4"),
                        ("uncertainty", "<f4"),
                    ],
                )
                fftab = datamodels.MirImgPhotomModel(phot_table=data)
                self.photom_io(fftab.phot_table[0])
            except AttributeError:
                # No time-dependent correction is applied
                log.info(
                    " Skipping MIRI imager time correction."
                    " Extension not found in the reference file."
                )
                self.photom_io(ftab.phot_table[row])

        # MRS detectors
        elif self.detector == "MIRIFUSHORT" or self.detector == "MIRIFULONG":
            # Make sure all NaNs have DO_NOT_USE flag set
            where_nan = np.isnan(ftab.data)
            ftab.dq[where_nan] = np.bitwise_or(ftab.dq[where_nan], dqflags.pixel["DO_NOT_USE"])

            # Compute the combined 2D sensitivity factors
            sens2d = ftab.data

            # Multiply the science data and uncertainty arrays by the 2D
            # sensitivity factors
            sens2d_squared = sens2d * sens2d
            if not self.inverse:
                self.input.data *= sens2d
                self.input.err *= sens2d
                self.input.var_poisson *= sens2d_squared
                self.input.var_rnoise *= sens2d_squared
            else:
                self.input.data /= sens2d
                self.input.err /= sens2d
                self.input.var_poisson /= sens2d_squared
                self.input.var_rnoise /= sens2d_squared

            if self.input.var_flat is not None and np.size(self.input.var_flat) > 0:
                process_var_flat = True
                if not self.inverse:
                    self.input.var_flat *= sens2d_squared
                else:
                    self.input.var_flat /= sens2d_squared
            else:
                process_var_flat = False

            # Update the science dq
            self.input.dq = np.bitwise_or(self.input.dq, ftab.dq)

            # Check if reference file contains time dependent correction

            try:
                ftab.getarray_noinit("timecoeff_ch1")
            except AttributeError:
                # Old style ref file; skip the correction
                log.info(
                    "Skipping MRS MIRI time correction.  "
                    "Extensions not found in the reference file."
                )
                self.mrs_time_correction = False

            if self.mrs_time_correction:
                log.info("Applying MRS IFU time dependent correction.")
                mid_time = self.input.meta.exposure.mid_time
                correction = miri_mrs.time_correction(self.input, self.detector, ftab, mid_time)
                inv_correction_sq = 1.0 / (correction * correction)
                self.input.data /= correction
                self.input.err /= correction
                self.input.var_poisson *= inv_correction_sq
                self.input.var_rnoise *= inv_correction_sq
                if process_var_flat:
                    self.input.var_flat *= inv_correction_sq
            else:
                log.info("Not applying MRS IFU time dependent correction.")

            # Retrieve the scalar conversion factor from the reference data
            conv_factor = ftab.meta.photometry.conversion_megajanskys

            # Store the conversion factors in the meta data
            self.input.meta.photometry.conversion_megajanskys = conv_factor
            self.input.meta.photometry.conversion_microjanskys = conv_factor * MJSR_TO_UJA2

            # Update BUNIT values for the science data and err
            if not self.inverse:
                self.input.meta.bunit_data = "MJy/sr"
                self.input.meta.bunit_err = "MJy/sr"
            else:
                self.input.meta.bunit_data = "DN/s"
                self.input.meta.bunit_err = "DN/s"

    def calc_nircam(self, ftab):
        """
        Apply photometric calibration data to dataset and update conversion factor.

        For NIRCAM, matching is based on FILTER and PUPIL.
        The routine will find the corresponding information in the reference
        file, apply the conversion factors, and store the scalar conversion
        factor in the output model. If wavelength-dependent conversion factors
        exist, they will be included in the calibration.
        For WFSS (grism) mode, the calibration information extracted from the
        reference file is applied to each slit instance in the science data.

        Parameters
        ----------
        ftab : `~jwst.datamodels.NrcImgPhotomModel` or `~jwst.datamodels.NrcWfssPhotomModel`
            NIRCam photom reference file data model.
        """
        # Handle WFSS data separately from regular imaging
        if isinstance(self.input, datamodels.MultiSlitModel) and self.exptype == "NRC_WFSS":
            # Loop over the WFSS slits, applying the correct photom ref data
            for slit in self.input.slits:
                log.info(f"Working on slit {slit.name}")
                self.slitnum += 1
                order = slit.meta.wcsinfo.spectral_order
                # TODO: If it's reasonable to hardcode the list of orders for Nircam WFSS,
                # the code matching the two rows can be taken outside the loop.
                fields_to_match = {"filter": self.filter, "pupil": self.pupil, "order": order}
                row = find_row(ftab.phot_table, fields_to_match)
                if row is None:
                    continue
                self.photom_io(ftab.phot_table[row])
        elif self.exptype == "NRC_TSGRISM":
            fields_to_match = {"filter": self.filter, "pupil": self.pupil, "order": self.order}
            row = find_row(ftab.phot_table, fields_to_match)
            if row is None:
                return
            self.photom_io(ftab.phot_table[row])
        else:
            # check for subarray in the phot_table: older files do not have it
            fields_to_match = {"filter": self.filter, "pupil": self.pupil}
            if "subarray" in ftab.phot_table.columns.names:
                log.info("Matching to subarray: %s", self.subarray)
                fields_to_match["subarray"] = self.subarray

            row = find_row(ftab.phot_table, fields_to_match)
            if row is None:
                return
            self.photom_io(ftab.phot_table[row])

    def calc_fgs(self, ftab):
        """
        Apply photometric calibration data to dataset and update conversion factor.

        For FGS, there is no matching required, because the instrument does
        not contain any filters or pupil wheel elements. The only mode is CLEAR.

        The routine will find the corresponding information in the reference
        file (which should have only a single row), apply it to the data, and
        write the conversion factor to the output model.

        Parameters
        ----------
        ftab : `~jwst.datamodels.FgsImgPhotomModel`
            FGS photom reference file data model.
        """
        # Read the first (and only) row in the reference file
        self.photom_io(ftab.phot_table[0])

    def calc_nrs_ifu_sens2d(self, area_data):
        """
        Create the 2-D wavelength and pixel area arrays.

        These are needed for constructing a NIRSpec IFU sensitivity map.

        Parameters
        ----------
        area_data : 1-D numpy.ndarray
            Array of pixel area values for the IFU slices.

        Returns
        -------
        wave2d : 2-D numpy.ndarray
            Array of wavelengths per pixel.
        area2d : 2-D numpy.ndarray
            Array of pixel area values.
        dqmap : 2-D numpy.ndarray
            Array of DQ flags per pixel.
        """
        from jwst.assign_wcs import nirspec  # for NIRSpec IFU data
        import gwcs

        microns_100 = 1.0e-4  # 100 microns, in meters

        # Create empty 2D arrays for the wavelengths and pixel areas
        wave2d = np.zeros_like(self.input.data)
        area2d = np.ones_like(self.input.data)

        # Create and initialize an array for the 2D dq map to be returned.
        # initialize all pixels to NON_SCIENCE, because operations below
        # only touch pixels within the bounding_box of each slice
        dqmap = np.zeros_like(self.input.dq) + dqflags.pixel["NON_SCIENCE"]

        # Get the list of WCSs for the IFU slices
        list_of_wcs = nirspec.nrs_ifu_wcs(self.input)

        # Loop over the slices
        for k, ifu_wcs in enumerate(list_of_wcs):
            # Construct array indexes for pixels in this slice
            x, y = gwcs.wcstools.grid_from_bounding_box(
                ifu_wcs.bounding_box, step=(1, 1), center=True
            )

            log.debug(f"Slice {k}: {x[0][0]} {x[-1][-1]} {y[0][0]} {y[-1][-1]}")

            # Get the world coords for all pixels in this slice
            coords = ifu_wcs(x, y)
            dq = dqmap[y.astype(int), x.astype(int)]
            wl = coords[2]
            # pull out the valid wavelengths and reset other array to not include
            # nan values
            valid = ~np.isnan(wl)
            wl = wl[valid]
            dq = dq[valid]
            x = x[valid]
            y = y[valid]
            dq[:] = 0

            if wl.max() < microns_100:
                log.info("Wavelengths in WCS table appear to be in meters")

            dqmap[y.astype(int), x.astype(int)] = dq
            # Insert the wavelength values for this slice into the
            # whole image array
            wave2d[y.astype(int), x.astype(int)] = wl
            # Insert the pixel area value for this slice into the
            # whole image array
            ar = np.ones_like(wl)
            ar[:] = area_data[np.where(area_data["slice_id"] == k)]["pixarea"][0]
            area2d[y.astype(int), x.astype(int)] = ar

        return wave2d, area2d, dqmap

    def photom_io(self, tabdata, order=None):
        """
        Combine photometric conversion factors and apply to the science dataset.

        Parameters
        ----------
        tabdata : FITS record
            Single row of data from reference table.

        order : int
            Spectral order number.
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
        unit_is_surface_brightness = True  # default
        try:
            conversion = tabdata["photmjsr"]  # unit is MJy / sr
        except KeyError:
            conversion = tabdata["photmj"]  # unit is MJy
            if isinstance(self.input, datamodels.MultiSlitModel):
                slit = self.input.slits[self.slitnum]
                if self.exptype in ["NRS_MSASPEC", "NRS_FIXEDSLIT"]:
                    srctype = self.source_type if self.source_type else slit.source_type
                else:
                    srctype = self.source_type
                if srctype is None or srctype.upper() != "POINT":
                    if slit.meta.photometry.pixelarea_steradians is None:
                        log.warning(
                            "Pixel area is None, so can't convert flux to surface brightness!"
                        )
                    else:
                        log.debug("Converting conversion factor from flux to surface brightness")
                        conversion /= slit.meta.photometry.pixelarea_steradians
                else:
                    conversion_uniform = conversion / slit.meta.photometry.pixelarea_steradians
                    unit_is_surface_brightness = False
            elif isinstance(self.input, datamodels.TSOMultiSpecModel):
                # output from extract1d should not require this area conversion
                unit_is_surface_brightness = False
            else:
                if self.source_type is None or self.source_type != "POINT":
                    if self.input.meta.photometry.pixelarea_steradians is None:
                        log.warning(
                            "Pixel area is None, so can't convert flux to surface brightness!"
                        )
                    else:
                        log.debug("Converting conversion factor from flux to surface brightness")
                        conversion /= self.input.meta.photometry.pixelarea_steradians
                else:
                    pixel_area_steradians = self.input.meta.photometry.pixelarea_steradians
                    conversion_uniform = conversion / pixel_area_steradians
                    unit_is_surface_brightness = False

        # Store the conversion factor in the meta data
        log.info(f"PHOTMJSR value: {conversion:.6g}")
        if isinstance(self.input, datamodels.MultiSlitModel):
            self.input.slits[self.slitnum].meta.photometry.conversion_megajanskys = conversion
            self.input.slits[self.slitnum].meta.photometry.conversion_microjanskys = (
                conversion * MJSR_TO_UJA2
            )
        elif isinstance(self.input, datamodels.TSOMultiSpecModel):
            # No place in the schema to store photometry info
            pass
        else:
            self.input.meta.photometry.conversion_megajanskys = conversion
            self.input.meta.photometry.conversion_microjanskys = conversion * MJSR_TO_UJA2

        # If the photom reference file is for spectroscopic data, the table
        # in the reference file should contain a "wavelength" column (among
        # other columns).
        try:
            wl_test = tabdata["wavelength"]
            is_spectroscopic = True
            del wl_test
        except KeyError:
            is_spectroscopic = False

        # Get the length of the relative response arrays in this row.  If the
        # nelem column is not present, we'll use the entire wavelength and
        # relresponse arrays.
        if is_spectroscopic:
            try:
                nelem = tabdata["nelem"]
            except KeyError:
                nelem = None
        else:
            nelem = None

        # For spectroscopic data, include the relative response array in
        # the flux conversion.
        no_cal = None
        if is_spectroscopic:
            waves = tabdata["wavelength"]
            relresps = tabdata["relresponse"]
            if nelem is not None:
                waves = waves[:nelem]
                relresps = relresps[:nelem]

            # Make sure waves and relresps are in increasing wavelength order
            if not np.all(np.diff(waves) > 0):
                index = np.argsort(waves)
                waves = waves[index].copy()
                relresps = relresps[index].copy()

            # Convert wavelengths from meters to microns, if necessary
            microns_100 = 1.0e-4  # 100 microns, in meters
            if waves.max() > 0.0 and waves.max() < microns_100:
                waves *= 1.0e6

            # Compute a 2-D grid of conversion factors, as a function of wavelength
            if isinstance(self.input, datamodels.MultiSlitModel):
                slit = self.input.slits[self.slitnum]
                # The NIRSpec fixed-slit primary slit needs special handling if
                # it contains a point source
                if self.exptype.upper() == "NRS_FIXEDSLIT" and slit.source_type.upper() == "POINT":
                    # First, compute 2D array of photom correction values using
                    # uncorrected wavelengths, which is appropriate for a uniform source
                    conversion_2d_uniform, no_cal = self.create_2d_conversion(
                        slit,
                        self.exptype,
                        conversion_uniform,
                        waves,
                        relresps,
                        order,
                        use_wavecorr=False,
                    )
                    slit.photom_uniform = conversion_2d_uniform  # store the result

                    # Now repeat the process using corrected wavelength values,
                    # which is appropriate for a point source. This is the version of
                    # the correction that will actually get applied to the data below.
                    conversion, no_cal = self.create_2d_conversion(
                        slit, self.exptype, conversion, waves, relresps, order, use_wavecorr=True
                    )
                    slit.photom_point = conversion  # store the result

                elif self.exptype in ["NRC_WFSS", "NRC_TSGRISM"]:
                    log.info("Including spectral dispersion in 2-d flux calibration")
                    conversion, no_cal = self.create_2d_conversion(
                        slit,
                        self.exptype,
                        conversion,
                        waves,
                        relresps,
                        order,
                        include_dispersion=True,
                    )

                else:
                    conversion, no_cal = self.create_2d_conversion(
                        slit, self.exptype, conversion, waves, relresps, order
                    )

            elif isinstance(self.input, datamodels.TSOMultiSpecModel):
                # This input does not require a 2d conversion, but a 1d interpolation on the
                # input wavelength vector to find the relresponse.
                conversion, no_cal = self.create_1d_conversion(
                    self.input.spec[self.specnum], conversion, waves, relresps, self.integ_row
                )
            else:
                # NRC_TSGRISM data produces a SpecModel, which is handled here
                if self.exptype in ["NRC_WFSS", "NRC_TSGRISM"]:
                    log.info("Including spectral dispersion in 2-d flux calibration")
                    conversion, no_cal = self.create_2d_conversion(
                        self.input,
                        self.exptype,
                        conversion,
                        waves,
                        relresps,
                        order,
                        include_dispersion=True,
                    )

                else:
                    conversion, no_cal = self.create_2d_conversion(
                        self.input, self.exptype, conversion, waves, relresps, order
                    )
        # Apply the conversion to the data and all uncertainty arrays
        if isinstance(self.input, datamodels.MultiSlitModel):
            slit = self.input.slits[self.slitnum]
            conversion_squared = conversion * conversion
            if not self.inverse:
                slit.data *= conversion
                slit.err *= conversion
            else:
                slit.data /= conversion
                slit.err /= conversion
            if slit.var_poisson is not None and np.size(slit.var_poisson) > 0:
                if not self.inverse:
                    slit.var_poisson *= conversion_squared
                else:
                    slit.var_poisson /= conversion_squared
            if slit.var_rnoise is not None and np.size(slit.var_rnoise) > 0:
                if not self.inverse:
                    slit.var_rnoise *= conversion_squared
                else:
                    slit.var_rnoise /= conversion_squared
            if slit.var_flat is not None and np.size(slit.var_flat) > 0:
                if not self.inverse:
                    slit.var_flat *= conversion_squared
                else:
                    slit.var_flat /= conversion_squared
            if no_cal is not None:
                slit.dq[..., no_cal] = np.bitwise_or(
                    slit.dq[..., no_cal], dqflags.pixel["DO_NOT_USE"]
                )
            if not self.inverse:
                if unit_is_surface_brightness:
                    slit.meta.bunit_data = "MJy/sr"
                    slit.meta.bunit_err = "MJy/sr"
                else:
                    slit.meta.bunit_data = "MJy"
                    slit.meta.bunit_err = "MJy"
                # Setting top model to None so they will not be written to FITs File.
                # Information on the units should only come from the individual slits.
                self.input.meta.bunit_data = None
                self.input.meta.bunit_err = None
            else:
                self.input.meta.bunit_data = "DN/s"
                self.input.meta.bunit_err = "DN/s"

            # Make sure output model has consistent NaN and DO_NOT_USE values
            match_nans_and_flags(slit)

        elif isinstance(self.input, datamodels.TSOMultiSpecModel):
            # Does this block need to address SB columns as well, or will
            # they (presumably) never be populated for SOSS?
            # It appears flux_error is the only error column populated?
            spec = self.input.spec[self.specnum]
            spec.spec_table.FLUX[self.integ_row] *= conversion
            spec.spec_table.FLUX_ERROR[self.integ_row] *= conversion
            spec.spec_table.FLUX_VAR_POISSON[self.integ_row] *= conversion**2.0
            spec.spec_table.FLUX_VAR_RNOISE[self.integ_row] *= conversion**2.0
            spec.spec_table.FLUX_VAR_FLAT[self.integ_row] *= conversion**2.0
            spec.spec_table.BACKGROUND[self.integ_row] *= conversion
            spec.spec_table.BKGD_ERROR[self.integ_row] *= conversion
            spec.spec_table.BKGD_VAR_POISSON[self.integ_row] *= conversion**2.0
            spec.spec_table.BKGD_VAR_RNOISE[self.integ_row] *= conversion**2.0
            spec.spec_table.BKGD_VAR_FLAT[self.integ_row] *= conversion**2.0
            spec.spec_table.columns["FLUX"].unit = "MJy"
            spec.spec_table.columns["FLUX_ERROR"].unit = "MJy"
            spec.spec_table.columns["FLUX_VAR_POISSON"].unit = "MJy^2"
            spec.spec_table.columns["FLUX_VAR_RNOISE"].unit = "MJy^2"
            spec.spec_table.columns["FLUX_VAR_FLAT"].unit = "MJy^2"
            spec.spec_table.columns["BACKGROUND"].unit = "MJy"
            spec.spec_table.columns["BKGD_ERROR"].unit = "MJy"
            spec.spec_table.columns["BKGD_VAR_POISSON"].unit = "MJy^2"
            spec.spec_table.columns["BKGD_VAR_RNOISE"].unit = "MJy^2"
            spec.spec_table.columns["BKGD_VAR_FLAT"].unit = "MJy^2"

        else:
            conversion_squared = conversion * conversion
            if not self.inverse:
                self.input.data *= conversion
                self.input.err *= conversion

            else:
                self.input.data /= conversion
                self.input.err /= conversion
            if self.input.var_poisson is not None and np.size(self.input.var_poisson) > 0:
                if not self.inverse:
                    self.input.var_poisson *= conversion_squared
                else:
                    self.input.var_poisson /= conversion_squared
            if self.input.var_rnoise is not None and np.size(self.input.var_rnoise) > 0:
                if not self.inverse:
                    self.input.var_rnoise *= conversion_squared
                else:
                    self.input.var_rnoise /= conversion_squared
            if self.input.var_flat is not None and np.size(self.input.var_flat) > 0:
                if not self.inverse:
                    self.input.var_flat *= conversion_squared
                else:
                    self.input.var_flat /= conversion_squared
            if no_cal is not None:
                self.input.dq[..., no_cal] = np.bitwise_or(
                    self.input.dq[..., no_cal], dqflags.pixel["DO_NOT_USE"]
                )

            if not self.inverse:
                if unit_is_surface_brightness:
                    self.input.meta.bunit_data = "MJy/sr"
                    self.input.meta.bunit_err = "MJy/sr"
                else:
                    self.input.meta.bunit_data = "MJy"
                    self.input.meta.bunit_err = "MJy"
            else:
                self.input.meta.bunit_data = "DN/s"
                self.input.meta.bunit_err = "DN/s"

            # Make sure output model has consistent NaN and DO_NOT_USE values
            match_nans_and_flags(self.input)

    def create_2d_conversion(
        self,
        model,
        exptype,
        conversion,
        waves,
        relresps,
        order,
        use_wavecorr=None,
        include_dispersion=False,
    ):
        """
        Create a 2D array of photometric conversion values.

        This array is based on wavelengths per pixel and response
        as a function of wavelength.

        Parameters
        ----------
        model : `~jwst.datamodels.JwstDataModel`
            Input data model containing the necessary wavelength information.
        exptype : str
            Exposure type of the input.
        conversion : float
            Initial scalar photometric conversion value.
        waves : float numpy.ndarray
            1D wavelength vector on which relative response values are
            sampled.
        relresps : float numpy.ndarray
            1D photometric response values, as a function of waves.
        order : int
            Spectral order number.
        use_wavecorr : bool or None
            Flag indicating whether or not to use corrected wavelengths.
            Typically only used for NIRSpec fixed-slit data.
        include_dispersion : bool or None
            Flag indicating whether the dispersion needs to be incorporated
            into the 2-d conversion factors.

        Returns
        -------
        conversion : float numpy.ndarray
            2D array of computed photometric conversion values.
        no_cal : int numpy.ndarray
            2D mask indicating where no conversion is available.
        """
        # Get the 2D wavelength array corresponding to the input
        # image pixel values
        wl_array = get_wavelengths(model, exptype, order, use_wavecorr)
        wl_array[np.isnan(wl_array)] = -1.0

        # Interpolate the photometric response values onto the
        # 2D wavelength grid
        # waves is in microns, so wl_array must be in microns too
        conv_2d = np.interp(wl_array, waves, relresps, left=np.nan, right=np.nan)

        if include_dispersion:
            dispaxis = get_dispersion_direction(self.exptype, self.grating, self.filter, self.pupil)
            if dispaxis is not None:
                dispersion_array = self.get_dispersion_array(wl_array, dispaxis)
                # Convert dispersion from micron/pixel to angstrom/pixel
                dispersion_array *= 1.0e4
                conv_2d /= np.abs(dispersion_array)
            else:
                log.warning(
                    "Unable to get dispersion direction, so cannot calculate dispersion array"
                )
        # Combine the scalar and 2D conversion factors
        conversion = conversion * conv_2d
        no_cal = np.isnan(conv_2d)
        conversion[no_cal] = 0.0

        return conversion, no_cal

    def get_dispersion_array(self, wavelength_array, dispaxis):
        """
        Create an array of dispersion values from the wavelength array.

        Parameters
        ----------
        wavelength_array : float
            2-d array of wavelength values, assumed to be in microns.
        dispaxis : int
            Direction along which light is dispersed: 1 = along rows, 2 = along columns.

        Returns
        -------
        dispersion_array : float
            2-d array of dispersion values, in microns/pixel.
        """
        nrows, ncols = wavelength_array.shape
        dispersion_array = np.zeros(wavelength_array.shape)
        if dispaxis == 1:
            for row in range(nrows):
                dispersion_array[row] = np.gradient(wavelength_array[row])
        elif dispaxis == 2:
            for column in range(ncols):
                dispersion_array[:, column] = np.gradient(wavelength_array[:, column])
        else:
            log.warning(f"Can't process data with DISPAXIS={dispaxis}")
        return dispersion_array

    def create_1d_conversion(self, model, conversion, waves, relresps, integ_row):
        """
        Resample the photometric conversion array.

        Create a 1D array of photometric conversion values based on
        wavelength array of input spectrum and response as a function of wavelength.

        Parameters
        ----------
        model : `~jwst.datamodels.JwstDataModel`
            Input data model containing the necessary wavelength information.
        conversion : float
            Initial scalar photometric conversion value.
        waves : float numpy.ndarray
            1D wavelength vector on which relative response values are
            sampled.
        relresps : float numpy.ndarray
            1D photometric response values, as a function of waves.
        integ_row : int
            Table row number for the spectrum for the current integration.

        Returns
        -------
        conversion : float numpy.ndarray
            1D array of computed photometric conversion values.
        no_cal : int numpy.ndarray
            1D mask indicating where no conversion is available.
        """
        # Get the 2D wavelength array corresponding to the input
        # image pixel values
        wl_array = model.spec_table["WAVELENGTH"][integ_row]

        flip_wl = False
        if np.nanargmax(wl_array) - np.nanargmin(wl_array) < 0:
            # Need monotonically increasing wavelengths for interp
            # Bool flag to flip fit if True
            flip_wl = True
            wl_array = wl_array[::-1]

        wl_array[np.isnan(wl_array)] = -1.0

        if np.nanargmax(waves) - np.nanargmin(waves) < 0:
            # Need monotonically increasing wavelengths for interp
            # This shouldn't have effects external to method.
            waves = waves[::-1]
            relresps = relresps[::-1]

        # Interpolate the photometric response values onto the
        # 1D wavelength grid
        conv_1d = np.interp(wl_array, waves, relresps, left=np.nan, right=np.nan)

        if flip_wl:
            # If wl_array was flipped, flip the conversion before returning it.
            conv_1d = conv_1d[::-1]
        # Combine the scalar and 1D conversion factors
        conversion = conversion * conv_1d
        no_cal = np.isnan(conv_1d)
        conversion[no_cal] = 0.0

        return conversion, no_cal

    def pixarea_from_ftab(self, ftab):
        """
        Get pixel area in steradians and arsec^2 from photom reference file.

        Read the pixel area values in the PIXAR_A2 and PIXAR_SR keys from the
        primary header of the photom reference file.

        Parameters
        ----------
        ftab : `~jwst.datamodels.JwstDataModel`
            A photom reference file data model.

        Returns
        -------
        area_ster : float
            Pixel area in steradians
        area_a2 : float
            Pixel area in arcsec^2
        """
        area_ster, area_a2 = None, None
        area_ster = ftab.meta.photometry.pixelarea_steradians
        log.info("Attempting to obtain PIXAR_SR and PIXAR_A2 values from PHOTOM reference file.")
        if area_ster is None:
            log.warning("The PIXAR_SR keyword is missing from %s", ftab.meta.filename)
        area_a2 = ftab.meta.photometry.pixelarea_arcsecsq
        if area_a2 is None:
            log.warning("The PIXAR_A2 keyword is missing from %s", ftab.meta.filename)
        if area_ster is not None and area_a2 is not None:
            log.info("Values for PIXAR_SR and PIXAR_A2 obtained from PHOTOM reference file.")
        return area_ster, area_a2

    def save_area_info(self, ftab, area_fname):
        """
        Populate area attributes in science dataset.

        Read the pixel area values in the PIXAR_A2 and PIXAR_SR keys from the
        primary header of the pixel area reference file or (for NIRSpec data)
        from the PIXAREA column in the selected row of the AREA table in the
        area reference file.  Use that information to populate the pixel area
        keywords in the output product.

        Except for NIRSpec data, also copy the pixel area data array from the
        pixel area reference file to the area extension of the output product.

        Parameters
        ----------
        ftab : `~jwst.datamodels.JwstDataModel`
            A photom reference file data model.

        area_fname : str
            Pixel area reference file name.
        """
        use_pixarea_rfile = False
        area_ster, area_a2 = None, None
        if area_fname is not None and area_fname != "N/A":
            use_pixarea_rfile = True
            # Load the pixel area reference file
            pix_area = datamodels.open(area_fname)

        if self.instrument != "NIRSPEC":
            if use_pixarea_rfile:
                # Copy the pixel area data array to the appropriate attribute
                # of the science data model
                if isinstance(self.input, datamodels.MultiSlitModel):
                    # Note that this only copied to the first slit.
                    self.input.slits[0].area = pix_area.data
                else:
                    ystart = self.input.meta.subarray.ystart - 1
                    xstart = self.input.meta.subarray.xstart - 1
                    yend = ystart + self.input.meta.subarray.ysize
                    xend = xstart + self.input.meta.subarray.xsize
                    self.input.area = pix_area.data[ystart:yend, xstart:xend]
                log.info("Pixel area map copied to output.")

                # Load the average pixel area values from the AREA reference file header
                # Don't need to do this for NIRSpec, because pixel areas will be
                # copied using save_area_nirspec
                try:
                    area_ster = pix_area.meta.photometry.pixelarea_steradians
                    if area_ster is None:
                        log.warning("The PIXAR_SR keyword is missing from %s", area_fname)
                    area_a2 = pix_area.meta.photometry.pixelarea_arcsecsq
                    if area_a2 is None:
                        log.warning("The PIXAR_A2 keyword is missing from %s", area_fname)
                    if area_ster is not None and area_a2 is not None:
                        log.info(
                            "Values for PIXAR_SR and PIXAR_A2 obtained from AREA reference file."
                        )

                # The area reference file might be older, try the photom reference file
                except (AttributeError, KeyError):
                    area_ster, area_a2 = self.pixarea_from_ftab(ftab)

                pix_area.close()

            # The area reference file might be older, try the photom reference file
            else:
                area_ster, area_a2 = self.pixarea_from_ftab(ftab)

            # Copy the pixel area values to the output
            log.debug("PIXAR_SR = %s, PIXAR_A2 = %s", str(area_ster), str(area_a2))
            if area_a2 is not None:
                area_a2 = float(area_a2)
            if area_ster is not None:
                area_ster = float(area_ster)

            # For MultiSlitModels, copy to each slit attribute
            if isinstance(self.input, datamodels.MultiSlitModel):
                for slit in self.input.slits:
                    slit.meta.photometry.pixelarea_arcsecsq = area_a2
                    slit.meta.photometry.pixelarea_steradians = area_ster
            else:
                self.input.meta.photometry.pixelarea_arcsecsq = area_a2
                self.input.meta.photometry.pixelarea_steradians = area_ster
        else:
            self.save_area_nirspec(pix_area)
            pix_area.close()

    def save_area_nirspec(self, pix_area):
        """
        Populate the pixel area attributes of science dataset.

        Read the pixel area value from the PIXAREA column in the selected row
        of the AREA table in the area reference file.
        Use that information to populate the pixel area keywords in the output
        product.

        Parameters
        ----------
        pix_area : `~jwst.datamodels.JwstDataModel`
            Pixel area reference file data model.
        """
        exp_type = self.exptype
        pixarea = pix_area.area_table["pixarea"]
        if exp_type == "NRS_MSASPEC":
            quadrant = pix_area.area_table["quadrant"]
            shutter_x = pix_area.area_table["shutter_x"]
            shutter_y = pix_area.area_table["shutter_y"]
            n_failures = 0
            for slit in self.input.slits:
                match_q = quadrant == slit.quadrant
                match_x = shutter_x == slit.xcen
                match_y = shutter_y == slit.ycen
                match = np.logical_and(match_q, np.logical_and(match_x, match_y))
                n_matches = match.sum(dtype=np.int64)
                if n_matches != 1:
                    n_failures += 1
                    slit.meta.photometry.pixelarea_arcsecsq = 1.0
                    slit.meta.photometry.pixelarea_steradians = 1.0
                else:
                    slit.meta.photometry.pixelarea_arcsecsq = float(pixarea[match].item())
                    slit.meta.photometry.pixelarea_steradians = (
                        slit.meta.photometry.pixelarea_arcsecsq * A2_TO_SR
                    )
            if n_failures > 0:
                log.warning(
                    "%d failures out of %d, matching MSA data to area reference file",
                    n_failures,
                    len(self.input.slits),
                )

        elif exp_type == "NRS_BRIGHTOBJ":
            slit_id = pix_area.area_table["slit_id"]
            nrows = len(slit_id)
            slit_name = self.input.name  # "S1600A1"
            foundit = False
            for k in range(nrows):
                if slit_id[k] == slit_name:
                    foundit = True
                    self.input.meta.photometry.pixelarea_arcsecsq = float(pixarea[k])
                    self.input.meta.photometry.pixelarea_steradians = (
                        self.input.meta.photometry.pixelarea_arcsecsq * A2_TO_SR
                    )
                    break
            if not foundit:
                log.warning("%s not found in pixel area table", slit_name)
                self.input.meta.photometry.pixelarea_arcsecsq = 1.0
                self.input.meta.photometry.pixelarea_steradians = 1.0

        elif exp_type in ["NRS_LAMP", "NRS_FIXEDSLIT"]:
            slit_id = pix_area.area_table["slit_id"]
            nrows = len(slit_id)
            for slit in self.input.slits:
                foundit = False
                for k in range(nrows):
                    if slit_id[k] == slit.name:
                        foundit = True
                        slit.meta.photometry.pixelarea_arcsecsq = float(pixarea[k])
                        slit.meta.photometry.pixelarea_steradians = (
                            slit.meta.photometry.pixelarea_arcsecsq * A2_TO_SR
                        )
                        break
                if not foundit:
                    log.warning("%s not found in pixel area table", slit.name)
                    slit.meta.photometry.pixelarea_arcsecsq = 1.0
                    slit.meta.photometry.pixelarea_steradians = 1.0

        elif exp_type == "NRS_IFU":
            # There is a slice_id column for selecting a matching slice, but
            # we're just going to average the pixel area for all slices.
            pixel_area = np.nanmean(pixarea)
            self.input.meta.photometry.pixelarea_arcsecsq = pixel_area
            self.input.meta.photometry.pixelarea_steradians = (
                self.input.meta.photometry.pixelarea_arcsecsq * A2_TO_SR
            )

        else:
            log.warning(
                "EXP_TYPE of NIRSpec data is %s, which is not an "
                "expected value; pixel area keywords will be set to 1.0",
                exp_type,
            )
            self.input.meta.photometry.pixelarea_arcsecsq = 1.0
            self.input.meta.photometry.pixelarea_steradians = 1.0

    def apply_photom(self, photom_fname, area_fname):
        """
        Apply the photom calibration step.

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
            Photom reference file name.
        area_fname : str
            Pixel area map reference file name.

        Returns
        -------
        output_model : ~jwst.datamodels.JwstDataModel
            Output data model with the flux calibrations applied.
        """
        with datamodels.open(photom_fname) as ftab:
            # Load the pixel area reference file, if it exists, and attach the
            # reference data to the science model
            # SOSS data are in a TSOMultiSpecModel, which will not allow for
            # saving the area info.
            if self.exptype != "NIS_SOSS":
                self.save_area_info(ftab, area_fname)

            if self.instrument == "NIRISS":
                self.calc_niriss(ftab)

            elif self.instrument == "NIRSPEC":
                # Some NIRSpec flat variance values can overflow when multiplied by the
                # flux conversion factor.  Photom values for some regions are also set to zero.
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", "overflow encountered", RuntimeWarning)
                    warnings.filterwarnings("ignore", "divide by zero", RuntimeWarning)
                    warnings.filterwarnings("ignore", "invalid value", RuntimeWarning)
                    self.calc_nirspec(ftab, area_fname)

            elif self.instrument == "NIRCAM":
                self.calc_nircam(ftab)

            elif self.instrument == "MIRI":
                # MIRI data may have NaNs in data or error arrays
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", "invalid value", RuntimeWarning)
                    self.calc_miri(ftab)

            elif self.instrument == "FGS":
                self.calc_fgs(ftab)

            else:
                raise RuntimeError(f"Instrument {self.instrument} is not recognized")

        return self.input
