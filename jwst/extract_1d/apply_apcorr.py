import abc

from scipy.interpolate import RectBivariateSpline, interp1d
from stdatamodels.jwst.datamodels import MultiSlitModel

from jwst.assign_wcs.util import compute_scale
import numpy as np


class ApCorrBase(abc.ABC):
    """Base class for aperture correction classes."""

    match_pars = {
        "MIRI": {"LRS": {"subarray": ["name"]}},
        "NIRSPEC": {
            "MSASPEC": {"instrument": ["filter", "grating"]},
            "FIXEDSLIT": {
                "instrument": ["filter", "grating"]
            },  # Slit is also required; passed in as init arg
            "BRIGHTOBJ": {"instrument": ["filter", "grating"]},
        },
        "NIRCAM": {
            "WFSS": {"instrument": ["filter", "pupil"]},
            "NRC_GRISM": {"instrument": ["filter", "pupil"]},
        },
        "NIRISS": {"WFSS": {"instrument": ["filter", "pupil"]}},
    }

    def __init__(
        self, input_model, apcorr_table, sizeunit, location=None, slit_name=None, **match_kwargs
    ):
        """
        Create and apply an aperture correction.

        Parameters
        ----------
        input_model : `~/jwst.datamodels.JwstDataModel`
            Input data model used to determine matching parameters.
        apcorr_table : `~/astropy.io.fits.FITS_rec`
            Aperture correction table data from APCORR reference file.
        sizeunit : str
            Units for the aperture correction data "row" values (assuming the
            aperture correction is a 2D array).
        location : tuple, optional
            Reference location (RA, DEC) used to calculate the pixel scale
            used in converting values in arcec to pixels.
            Default is None, however, if the input reference data contains
            size/radius data in units of arcsecs, a location is required.
        slit_name : str, optional
            For MultiSlitModels, the name of the slit being processed.
        **match_kwargs : dict, optional
            Additional keywords for matching data to reference table entries.

        Raises
        ------
        ValueError :
            If apcorr_row_units are not supplied and units are undefined in apcorr_table ColDefs.
        ValueError :
            If the input apcorr_table cannot be reduced to a single row
            based on match criteria from input_model.
        """
        self.correction = None

        self.model = input_model
        self._reference_table = apcorr_table
        self.location = location
        self.apcorr_sizeunits = sizeunit
        self.slit_name = slit_name

        self.match_keys = self._get_match_keys()
        self.match_pars = self._get_match_pars()
        self.match_pars.update(match_kwargs)
        self.reference = self._reduce_reftable()
        self._convert_size_units()
        self.apcorr_func = self.approximate()
        self.tabulated_correction = None

    @property
    def size_key(self):
        """
        Size key for the reference.

        The value is intended to index into self.reference.  Implementations
        of this abstract class should define this property as appropriate.
        """
        return None

    def _convert_size_units(self):
        """If the SIZE or Radius column is in units of arcseconds, convert to pixels."""
        if self.apcorr_sizeunits.startswith("arcsec"):
            # compute_scale returns scale in degrees
            if self.location is not None:
                if isinstance(self.model, MultiSlitModel):
                    idx = [slit.name for slit in self.model.slits].index(self.slit_name)
                    scale_degrees = compute_scale(
                        self.model.slits[idx].meta.wcs,
                        self.location,
                        disp_axis=self.model.slits[idx].meta.wcsinfo.dispersion_direction,
                    )
                    scale_arcsec = scale_degrees * 3600.00
                    self.reference[self.size_key] /= scale_arcsec
                else:
                    scale_degrees = compute_scale(
                        self.model.meta.wcs,
                        self.location,
                        disp_axis=self.model.meta.wcsinfo.dispersion_direction,
                    )
                    scale_arcsec = scale_degrees * 3600.00
                    self.reference[self.size_key] /= scale_arcsec
            else:
                raise ValueError(
                    "If the size column for the input APCORR reference "
                    "file is in units with arcseconds, a location "
                    "(RA, DEC, wavelength) must be provided in order to "
                    "compute a pixel scale to convert arcseconds to "
                    "pixels."
                )

    def _get_match_keys(self):
        """
        Get column keys needed for reducing the reference table based on input.

        Returns
        -------
        dict
            Match keys relevant to the current instrument and exposure type.
        """
        instrument = self.model.meta.instrument.name.upper()
        exptype = self.model.meta.exposure.type.upper()

        relevant_pars = self.match_pars[instrument]

        for key in relevant_pars.keys():
            if key in exptype:
                return relevant_pars[key]

    def _get_match_pars(self):
        """
        Get meta parameters required for reference table row-selection.

        Returns
        -------
        dict
            Match meta-parameters for the current match keys.
        """
        match_pars = {}

        for node, keys in self.match_keys.items():
            meta_node = getattr(self.model.meta, node)

            for key in keys:
                match_pars[key if key != "name" else node] = getattr(meta_node, key)

        return match_pars

    def _reduce_reftable(self):
        """
        Reduce full reference table to a single matched row.

        Returns
        -------
        `~fits.FITS_record`
            A single row from the reference table.
        """
        table = self._reference_table.copy()

        for key, value in self.match_pars.items():
            if isinstance(value, str):
                # Not all files will have the same format as input model metadata values.
                table = table[table[key].upper() == value.upper()]
            else:
                table = table[table[key] == value]

        if len(table) != 1:
            raise ValueError("Could not resolve APCORR reference for input.")
        return table[0]

    @abc.abstractmethod
    def approximate(self):
        """Generate an approximate aperture correction function based on input data."""
        pass

    def apply(self, spec_table):
        """
        Apply interpolated aperture correction to extraction results in-place.

        Parameters
        ----------
        spec_table : `~fits.FITS_rec`
            Table of aperture corrections values from apcorr reference file.
        """
        flux_cols_to_correct = ("flux", "flux_error", "surf_bright", "sb_error")
        var_cols_to_correct = (
            "flux_var_poisson",
            "flux_var_rnoise",
            "flux_var_flat",
            "sb_var_poisson",
            "sb_var_rnoise",
            "sb_var_flat",
        )

        for row in spec_table:
            correction = self.apcorr_func(row["npixels"], row["wavelength"])

            for col in flux_cols_to_correct:
                row[col] *= correction.item()
            for col in var_cols_to_correct:
                row[col] *= correction.item() * correction.item()


class ApCorrPhase(ApCorrBase):
    """Produce and apply aperture correction for input data with pixel phase."""

    size_key = "size"

    def __init__(self, *args, pixphase=0.5, **kwargs):
        """
        Create an aperture correction for data with pixel phase.

        Parameters
        ----------
        *args : list
            Input arguments.  See `ApCorrBase` for more.
        pixphase : float, optional
            Pixel phase of the input data.
        **kwargs : dict, optional
            Additional parameters. See `ApCorrBase` for more.
        """
        self.phase = pixphase

        super().__init__(*args, **kwargs)

    def approximate(self):
        """
        Generate an approximate function for interpolating apcorr values.

        Returns
        -------
        callable
            The approximation function.
        """

        def _approx_func(wavelength, size, pixel_phase):
            """
            Create a function to approximate the aperture correction in two stages.

            Parameters
            ----------
            wavelength : float
                Input wavelength
            size : float
                Input size (In extract.py this would be `n_pixels`)
            pixel_phase : float
                Input pixel phase

            Returns
            -------
            callable
                Aperture correction approximation function that takes
                wavelength, size, and pixel_phase as inputs.
            """
            # apcorr column data has shape (pixphase, wavelength, size)
            # Reduce apcorr dimensionality by interpolating in the pixphase
            # dimension first, then size & wavelength
            apcorr_pixphase_func = interp1d(self.reference["pixphase"], self.reference["apcorr"])
            size_wl_func = interp1d(self.reference["wavelength"], self.reference["size"])

            apcorr_pixphase = apcorr_pixphase_func(pixel_phase)
            size_wl = size_wl_func(wavelength)

            # by default RectBivariateSpline is 3rd order,
            # fails for size_wl=3 as in e.g. the test data
            wl_sortidx = np.argsort(self.reference["wavelength"])
            apcorr_pixphase = apcorr_pixphase[:, wl_sortidx]
            wl_ref = self.reference["wavelength"][wl_sortidx]
            pixphase_size_func = RectBivariateSpline(wl_ref, size_wl, apcorr_pixphase.T, ky=1, kx=1)
            size_func = pixphase_size_func(wavelength, size).T

            return size_func

        return _approx_func

    def tabulate_correction(self, spec_table):
        """
        Tabulate the interpolated aperture correction value.

        Storing the values saves time when applying it later, especially if
        it is to be applied to multiple integrations.

        Modifies self.tabulated_correction.

        Parameters
        ----------
        spec_table : `~fits.FITS_rec`
            Table of aperture corrections values from apcorr reference file.
        """
        coefs = []
        for row in spec_table:
            try:
                correction = self.apcorr_func(row["wavelength"], row["npixels"], self.phase)
            except ValueError:
                # Some input wavelengths might not be supported
                # (especially at the ends of the range)
                correction = None

            if correction:
                coefs += [correction.item()]
            else:
                coefs += [1]

        self.tabulated_correction = np.asarray(coefs)

    def apply(self, spec_table, use_tabulated=False):
        """
        Apply interpolated aperture correction value to source-related extraction results in-place.

        Parameters
        ----------
        spec_table : `~fits.FITS_rec`
            Table of aperture corrections values from apcorr reference file.
        use_tabulated : bool, optional
            Use self.tabulated_correction to perform the aperture correction?
            Default False (recompute correction from scratch).
        """
        flux_cols_to_correct = ("flux", "flux_error", "surf_bright", "sb_error")
        var_cols_to_correct = (
            "flux_var_poisson",
            "flux_var_rnoise",
            "flux_var_flat",
            "sb_var_poisson",
            "sb_var_rnoise",
            "sb_var_flat",
        )

        if use_tabulated:
            if self.tabulated_correction is None:
                raise ValueError(
                    "Cannot call apply_tabulated_correction without first "
                    "calling tabulate_correction"
                )

            for col in flux_cols_to_correct:
                spec_table[col] *= self.tabulated_correction
            for col in var_cols_to_correct:
                spec_table[col] *= self.tabulated_correction**2
        else:
            for row in spec_table:
                try:
                    correction = self.apcorr_func(row["wavelength"], row["npixels"], self.phase)
                except ValueError:
                    # Some input wavelengths might not be supported
                    # (especially at the ends of the range)
                    correction = None

                if correction:
                    for col in flux_cols_to_correct:
                        row[col] *= correction.item()
                    for col in var_cols_to_correct:
                        row[col] *= correction.item() * correction.item()


class ApCorrRadial(ApCorrBase):
    """Aperture correction class for spectra produced from an extraction aperture radius."""

    def __init__(self, input_model, apcorr_table, location=None):
        self.correction = None
        self.model = input_model
        self.location = location
        self.reference = apcorr_table.apcorr_table
        self.apcorr_sizeunits = self.reference.radius_units
        self._convert_size_units()

        # Set up some placeholders
        self.apcorr = None
        self.apcorr_size = None
        self.apcorr_correction = None

    def approximate(self):
        """
        Implement the approximation method.

        Has no effect for this class.  The base class requires that it be implemented,
        but the equivalent functionality for this class is performed in `find_apcorr_func`.
        """
        pass

    def _convert_size_units(self):
        """If the SIZE or Radius column is in units of arcseconds, convert to pixels."""
        if self.apcorr_sizeunits.startswith("arcsec"):
            # compute_scale returns scale in degrees
            if self.location is not None:
                scale_degrees = compute_scale(
                    self.model.meta.wcs,
                    self.location,
                    disp_axis=self.model.meta.wcsinfo.dispersion_direction,
                )
                scale_arcsec = scale_degrees * 3600.00
                self.reference.radius /= scale_arcsec
            else:
                raise ValueError(
                    "If the size column for the input APCORR reference file is "
                    "in units with arcseconds, a location "
                    "(RA, DEC, wavelength) must be provided in order to compute "
                    "a pixel scale to convert arcseconds to pixels."
                )

    def apply(self, spec_table):
        """
        Apply interpolated aperture correction to extraction results in-place.

        Parameters
        ----------
        spec_table : `~fits.FITS_rec`
            Table of aperture corrections values from apcorr reference file.
        """
        # check if MIRI data and correct the residual fringe flux and surface brightness columns.
        if self.model.meta.instrument.name == "MIRI":
            flux_cols_to_correct = (
                "flux",
                "flux_error",
                "surf_bright",
                "sb_error",
                "rf_flux",
                "rf_surf_bright",
            )
        else:
            flux_cols_to_correct = ("flux", "flux_error", "surf_bright", "sb_error")
        var_cols_to_correct = (
            "flux_var_poisson",
            "flux_var_rnoise",
            "flux_var_flat",
            "sb_var_poisson",
            "sb_var_rnoise",
            "sb_var_flat",
        )

        for i, row in enumerate(spec_table):
            correction = self.apcorr_correction[i]
            for col in flux_cols_to_correct:
                row[col] *= correction
            for col in var_cols_to_correct:
                row[col] *= correction * correction

    def match_wavelengths(self, wavelength_ifu):
        """
        Interpolate aperture correction and radial size onto input wavelengths.

        Interpolated correction and size values are stored in `self.apcorr`
        and `self.size`.

        Parameters
        ----------
        wavelength_ifu : ndarray
            Input wavelength array.
        """
        # given the ifu wavelength value - redefine the apcor func and radius to this wavelength
        # apcor reference data
        self.wavelength = self.reference.wavelength.flatten()
        self.size = self.reference.radius
        self.apcorr = self.reference.apcorr

        dim = self.apcorr.shape[0]
        size_match = np.zeros((dim, wavelength_ifu.shape[0]))
        apcorr_match = np.zeros((dim, wavelength_ifu.shape[0]))
        self.apcorr_correction = []  # set up here defined in find_apcor_func
        # loop over each radius dependent plane and interpolate to ifu wavelength
        for i in range(dim):
            radi = self.size[i, :]
            frad = interp1d(self.wavelength, radi, bounds_error=False, fill_value="extrapolate")
            radius_match = frad(wavelength_ifu)
            size_match[i, :] = radius_match

            appi = self.apcorr[i, :]
            fap = interp1d(self.wavelength, appi, bounds_error=False, fill_value="extrapolate")
            ap_match = fap(wavelength_ifu)
            apcorr_match[i, :] = ap_match

        self.apcorr = apcorr_match
        self.size = size_match

    def find_apcorr_func(self, iwave, radius_ifu):
        """
        Interpolate aperture correction onto a specific wavelength and radius.

        The correction is appended to `self.apcorr_correction`.

        Parameters
        ----------
        iwave : int
            Index for the wavelength.
        radius_ifu : float
            Radius to use.
        """
        # at ifu wavelength plane (iwave), the extraction radius is radius_ifu
        # pull out the radius values (self.size)  to use in the apcor ref file for this iwave
        # self.size and self.apcorr have already been interpolated in wavelength to match
        # the ifu wavelength range.

        radius_apcor = self.size[:, iwave]
        temparray = self.apcorr[:, iwave]
        fap = interp1d(radius_apcor, temparray, fill_value="extrapolate")
        correction = fap(radius_ifu)
        self.apcorr_correction.append(correction)


class ApCorr(ApCorrBase):
    """Default aperture correction class for use with most spectroscopic modes."""

    size_key = "size"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def approximate(self):
        """
        Generate an approximate function for interpolating apcorr values.

        Returns
        -------
        callable
            The approximation function.
        """
        wavelength = self.reference["wavelength"][: self.reference["nelem_wl"]]
        size = self.reference["size"][: self.reference["nelem_size"]]
        apcorr = self.reference["apcorr"][
            : self.reference["nelem_wl"], : self.reference["nelem_size"]
        ]

        # by default RectBivariateSpline is 3rd order,
        # fails for size_wl=3 as in e.g. the test data
        wl_sortidx = np.argsort(wavelength)
        apcorr = apcorr[wl_sortidx, :]
        wavelength = wavelength[wl_sortidx]

        return RectBivariateSpline(size, wavelength, apcorr.T, ky=1, kx=1)


def select_apcorr(input_model):
    """
    Select appropriate Aperture correction class based on input DataModel.

    Parameters
    ----------
    input_model : `~/jwst.datamodels.JwstDataModel`
        Input data on which the aperture correction is to be applied.

    Returns
    -------
    class
        Aperture correction class.
    """
    if input_model.meta.instrument.name == "MIRI":
        if "MRS" in input_model.meta.exposure.type:
            return ApCorrRadial
        else:
            return ApCorr

    if input_model.meta.instrument.name == "NIRCAM":
        return ApCorr

    if input_model.meta.instrument.name == "NIRISS":
        return ApCorr

    if input_model.meta.instrument.name == "NIRSPEC":
        if input_model.meta.exposure.type.upper() == "NRS_IFU":
            return ApCorrRadial
        else:
            return ApCorrPhase
