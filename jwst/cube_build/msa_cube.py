"""Work horse routines used for building MSA spectra cubes."""

import copy
import logging

import numpy as np
from gwcs import wcstools
from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels import dqflags

from jwst.assign_wcs import pointing
from jwst.assign_wcs.util import compute_footprint_nrs_slit
from jwst.cube_build import cube_lib
from jwst.cube_build.cube_match_sky_driz import cube_wrapper_driz  # c extension

log = logging.getLogger(__name__)

__all__ = ["MSACubeData"]


class MSACubeData:
    """
    Combine MSA slitlet data onto a regular grid.

    Parameters
    ----------
    **pars
        Dictionary of parameters controlling how the cube is built.
    """

    def __init__(
        self,
        **pars,
    ):

        self.source = pars.get("source")
        self.linear_wave = pars.get("linear_wave")
        self.scalexy = pars.get("scalexy")
        self.scalew = pars.get("scalew")
        self.slit_frac = pars.get("slit_frac")
        self.ra_center = pars.get("ra_center")
        self.dec_center = pars.get("dec_center")
        self.cube_pa = pars.get("cube_pa")
        self.nspax_x = pars.get("nspax_x")
        self.nspax_y = pars.get("nspax_y")
        self.debug_spaxel = pars.get("debug_spaxel")

        self.spaxel_x, self.spaxel_y, self.spaxel_z = [
            int(val) for val in self.debug_spaxel.split()
        ]

        self.coord_system = pars.get("coord_system")
        self.wavemin = pars.get("wavemin")
        self.wavemax = pars.get("wavemax")
        self.wavemin_ref = pars.get("wavemin_ref")
        self.wavemax_ref = pars.get("wavemax_ref")
        self.readnoise_weight = pars.get("readnoise_weight")
        self.wavelength_table = pars.get("wavelength_table")
        self.suffix = pars.get("suffix")

        self.cdelt1 = None
        self.cdelt2 = None
        self.cdelt3 = None
        self.rot_angle = None  # rotation angle between Ra-Dec and MSA local instrument plane

    def setup_msacube_wcs(self):
        """
        Set up the WCS of the MSA cube.

        Loops over every slitlet in the source dictionary, filters out bad pixels,
        and collects the spatial footprint corners and wavelength range of each slit.
        The first slit with more than ``_MIN_VALID_LAM_SIZE`` valid wavelength points
        is used to determine the rotation angle between the MSA plane and the sky.

        For ``msaalign`` cubes, the rotation angle is computed from the slit-frame-to-world
        transform of that slit. For ``skyalign`` cubes, the rotation angle is taken from
        ``cube_pa`` if provided, otherwise it defaults to ``None``.

        The final wavelength range is clamped to the reference file bounds unless the
        user has explicitly set ``wavemin`` or ``wavemax``, which take precedence.

        Sets
        ----
        self.cdelt1 : float
            Spatial pixel scale along axis 1 (arcseconds).
        self.cdelt2 : float
            Spatial pixel scale along axis 2 (arcseconds).
        self.cdelt3 : float
            Spectral pixel scale along the wavelength axis (microns), set only
            when ``linear_wave`` is True.
        self.rot_angle : float or None
            Rotation angle in degrees between the MSA plane and the sky.
            ``None`` if no slit had enough valid data or ``cube_pa`` was not set.
        self.cube_wcs : `~jwst.cube_build.cube_lib.CubeWCS`
            Data class containing the full WCS geometry of the output cube.

        Raises
        ------
        ValueError
            If ``lam_med`` is NaN for the slit selected for rotation, meaning the
            slit has no valid wavelength data after DQ filtering.
        ValueError
            If the slit-frame-to-world transform returns NaN sky coordinates,
            indicating an invalid WCS at the median wavelength.
        ValueError
            If ``dra`` and ``ddec`` are both zero, meaning the two slit reference
            points map to the same sky position.

        Notes
        -----
        For NIRSpec the units along/across the slice dimension are meters.

        If no slit exceeds ``_MIN_VALID_LAM_SIZE`` valid wavelength points,
        ``self.rot_angle`` will remain ``None`` and the cube will be built
        without a rotation correction.
        """
        self.cdelt1 = self.scalexy
        self.cdelt2 = self.scalexy
        if self.linear_wave:
            self.cdelt3 = self.scalew

        # lists used to determine min and max ranges
        corner_a = []
        corner_b = []
        lambda_min = []
        lambda_max = []

        # Define the rotation angle
        # If coord_system = msaalign then the angle is between the ra-dec and alpha beta
        # coord system using the first input model. Use first file in first band to set up
        # rotation angle.
        # Compute the rotation angle between local MSA system  and RA-DEC

        num = len(self.source["slit_models"])
        log.info(f"Number of slits:  {num}")

        for i in range(num):
            slit = self.source["slit_models"][i]
            wcs = copy.deepcopy(slit.meta.wcs)  # TODO do we need a deep copy ?
            bbox = wcs.bounding_box

            grid = wcstools.grid_from_bounding_box(bbox)
            ra, dec, lam = np.array(wcs(*grid))  # use to get lambda regions

            dq = slit.dq
            bad1 = np.bitwise_and(dq, dqflags.pixel["DO_NOT_USE"]).astype(bool)
            bad2 = np.bitwise_and(dq, dqflags.pixel["NON_SCIENCE"]).astype(bool)
            good_data = np.where(~bad1 & ~bad2)
            lam = lam[good_data]

            min_valid_lam_size = 1000  # is this the right value ?
            rotation_found = False
            # Use the first slit with enough valid data to determine the rotation angle.
            # The first slit is not always usable — it may have no valid wavelength data.
            use_for_rotation = (lam.size > min_valid_lam_size) and not rotation_found
            if use_for_rotation:
                rotation_found = True
                if self.coord_system == "msaalign":
                    s2w = wcs.get_transform("slit_frame", "world")
                    lam_med = np.nanmedian(lam)

                    if np.isnan(lam_med):
                        non_nan = np.sum(~np.isnan(lam))
                        raise ValueError(
                            f"lam_med is NaN for slit {i} — lam has {lam.size} elements, "
                            f"{non_nan} non-NaN. Slit may have no valid data."
                        )

                    temp_ra1, temp_dec1, _ = s2w(0, 0, lam_med)
                    temp_ra2, temp_dec2, _ = s2w(0, 0.005, lam_med)

                    if np.isnan([temp_ra1, temp_dec1, temp_ra2, temp_dec2]).any():
                        raise ValueError(
                            f"s2w transform returned NaN for slit {i} at lam_med={lam_med:.4f}. "
                            f"Got ra1={temp_ra1}, dec1={temp_dec1}, "
                            f"ra2={temp_ra2}, dec2={temp_dec2}. Check WCS validity."
                        )

                    dra = (temp_ra2 - temp_ra1) * np.cos(temp_dec1 * np.pi / 180.0)
                    ddec = temp_dec2 - temp_dec1

                    if dra == 0 and ddec == 0:
                        raise ValueError(
                            f"dra and ddec are both zero for slit {i} — "
                            "the two slit points map to the same sky position."
                        )

                    rot_angle = 90 + np.arctan2(dra, ddec) * 180.0 / np.pi
                    log.info(f"Rotation angle between MSA plane and sky: {rot_angle:.4f}")
                else:
                    if self.cube_pa is not None:
                        rot_angle = self.cube_pa  # None if not set, which is the correct default
                    else:
                        rot_angle = None
                self.rot_angle = rot_angle

            if len(good_data[0]) > 0:
                lmin = np.nanmin(lam)
                lmax = np.nanmax(lam)
                footprint, lam = compute_footprint_nrs_slit(slit)
                ca1 = float(footprint[0][0])
                cb1 = float(footprint[0][1])
                ca2 = float(footprint[1][0])
                cb2 = float(footprint[1][1])
                ca3 = float(footprint[2][0])
                cb3 = float(footprint[2][1])
                ca4 = float(footprint[3][0])
                cb4 = float(footprint[3][1])

                corner_a.append(ca1)
                corner_b.append(cb1)

                corner_a.append(ca2)
                corner_b.append(cb2)

                corner_a.append(ca3)
                corner_b.append(cb3)

                corner_a.append(ca4)
                corner_b.append(cb4)

                lambda_min.append(lmin)
                lambda_max.append(lmax)

        final_lam_min = min(lambda_min)
        final_lam_max = max(lambda_max)

        # final_lam_min and final_lam_max are determined from the data
        # Check if the user set wavemin. If so use that value. It takes precedence.
        if self.wavemin is not None:
            final_lam_min = self.wavemin
        # If the user has not set a wavemin value, then check that the one determined from the data
        # is not larger than the value given in the reference file.
        else:
            if self.wavemin_ref < final_lam_min:
                final_lam_min = (
                    self.wavemin_ref
                )  # TODO Check with NIRSPEC if we should have this check

        # user set wavemax use this values.
        if self.wavemax is not None:
            final_lam_max = self.wavemax
        # If the user has not set a wavemin value,
        # then check that the one determined from the data
        # is not larger than the value given in the reference file.
        else:
            if final_lam_max > self.wavemax_ref:
                final_lam_max = (
                    self.wavemax_ref
                )  # TODO Check with NIRSPEC if we should have this check

        self.cube_wcs = cube_lib.set_wcs(
            corner_a,
            corner_b,
            final_lam_min,
            final_lam_max,
            self.ra_center,
            self.dec_center,
            self.rot_angle,
            self.nspax_x,
            self.nspax_y,
            self.cdelt1,
            self.cdelt2,
            self.cdelt3,
            self.linear_wave,
            self.wavelength_table,
        )

        cube_lib.print_cube_geometry(self.coord_system, self.cube_wcs)

    def build_msacube(self):
        """
        Construct a 3D Micro-Shutter Array (MSA) spectral data cube.

        Iterates over all slit models for the current target, maps detector pixels
        to the target WCS output frame, and drizzles the data into a 3D spaxel grid.

        Processing Steps:
        1. Loop through all target slits (`self.source['slit_models']`).
        2. Map detector pixels to output coordinates via
        `map_source_pixels_to_output_frame`.
        3. Drizzle pixel fluxes onto the output cube grid using:
        - Accelerated C extension (`cube_wrapper_driz`)
        4. Compute final spaxel fluxes using `find_spaxel_flux`.
        5. Package flux, variance, weight, and DQ maps into an IFUCubeModel
        via `setup_final_model`.

        Returns
        -------
        final_cube : stdatamodels.jwst.datamodels.IFUCubeModel
            The fully combined 3D IFU data cube containing flux, variance,
            data quality, and spatial/spectral WCS headers.

        Notes
        -----
        - Spaxel debugging log outputs can be triggered by setting `self.spaxel_x`,
        `self.spaxel_y`, and `self.spaxel_z` prior to calling this method.
        """
        # self.cube_wcs = cube_wcs
        total_num = self.cube_wcs.naxis1 * self.cube_wcs.naxis2 * self.cube_wcs.naxis3
        self.spaxel_flux = np.zeros(total_num, dtype=np.float64)
        self.spaxel_weight = np.zeros(total_num, dtype=np.float64)
        self.spaxel_var = np.zeros(total_num, dtype=np.float64)
        self.spaxel_iflux = np.zeros(total_num, dtype=np.float64)
        self.spaxel_dq = np.zeros(total_num, dtype=np.uint32)
        self.spaxel_exposure_counter = np.zeros(total_num, dtype=np.uint32)

        # we need bunit in the out header - is there another way to do this ? #TODO
        # bunit_data = self.source["bunit"]
        # bunit_err = self.source["bunit_err"]
        nxyplane = self.cube_wcs.naxis1 * self.cube_wcs.naxis2

        if self.spaxel_z == -1 and self.spaxel_x == -1 and self.spaxel_y == -1:
            debug_spaxel_index = -1

        elif self.spaxel_z < 0 or self.spaxel_x < 0 or self.spaxel_y < 0:
            log.info("Incorrect input for Debug Spaxel values. Counting starts at 0")
            debug_spaxel_index = -1
            log.info(f"{self.spaxel_z} {self.spaxel_x}  {self.spaxel_y}")
        else:
            spaxel_z = self.spaxel_z
            spaxel_x = self.spaxel_x
            spaxel_y = self.spaxel_y
            debug_spaxel_index = spaxel_z * (nxyplane) + spaxel_y * self.naxis1 + spaxel_x
            log.info(
                f"Printing debug information for cube spaxel:  {spaxel_x} {spaxel_y} {spaxel_z}"
            )
            log.info(f"Debug spaxel {debug_spaxel_index}")

        num = self.source["number_slits"]

        # start looping over the slits for the source
        for i in range(num):
            print("Working on slit", i)
            slit = self.source["slit_models"][i]
            results = self.map_source_pixels_to_output_frame(slit)

            (
                x_array,
                y_array,
                flux,
                err,
                dq,
                var_rnoise,
                xi_array,
                eta_array,
                wave,
                dwave,
                corner,
                wave_pixel,
                delta_pixel,
            ) = results

            npt = flux.size
            # Check that the slit has valid pixels
            if npt > 0:
                cdelt3_mean = np.nanmean(self.cube_wcs.cdelt3_normal)
                xi1, eta1, xi2, eta2, xi3, eta3, xi4, eta4 = corner

            linear = 1
            dummy = x_array.copy()  # place holder needed for cube_build - valid for IFU data
            result = cube_wrapper_driz(
                self.readnoise_weight,
                self.cube_wcs.xcoord,
                self.cube_wcs.ycoord,
                self.cube_wcs.zcoord,
                xi_array,
                eta_array,
                wave,
                flux,
                err,
                dq,
                var_rnoise,
                dummy,
                xi1,
                eta1,
                xi2,
                eta2,
                xi3,
                eta3,
                xi4,
                eta4,
                dwave,
                self.cube_wcs.cdelt3_normal,
                self.cube_wcs.cdelt1,
                self.cube_wcs.cdelt2,
                cdelt3_mean,
                linear,
                x_array,
                y_array,
                debug_spaxel_index,
            )

            (
                this_spaxel_flux,
                this_spaxel_weight,
                this_spaxel_var,
                this_spaxel_iflux,
                this_spaxel_dq,
                this_spaxel_exp_counter,
            ) = result
            self.spaxel_flux = self.spaxel_flux + np.asarray(this_spaxel_flux, np.float64)
            self.spaxel_weight = self.spaxel_weight + np.asarray(this_spaxel_weight, np.float64)
            self.spaxel_var = self.spaxel_var + np.asarray(this_spaxel_var, np.float64)
            self.spaxel_iflux = self.spaxel_iflux + np.asarray(this_spaxel_iflux, np.float64)
            self.spaxel_dq.astype(np.uint)
            self.spaxel_dq = np.bitwise_or(self.spaxel_dq, this_spaxel_dq)
            self.spaxel_exposure_counter = self.spaxel_exposure_counter + np.asarray(
                this_spaxel_exp_counter, np.int32
            )
            result = None
            del result
            del (
                this_spaxel_flux,
                this_spaxel_weight,
                this_spaxel_var,
                this_spaxel_iflux,
                this_spaxel_dq,
                this_spaxel_exp_counter,
            )
            # end looping over the slits for the source

        self.spaxel_flux, self.spaxel_var = cube_lib.find_spaxel_flux(
            self.spaxel_iflux, self.spaxel_flux, self.spaxel_weight, self.spaxel_var
        )
        self.spaxel_flux, self.spaxel_var, self.spaxel_dq = cube_lib.set_final_dq_flags(
            self.spaxel_dq, self.spaxel_flux, self.spaxel_var, self.spaxel_weight
        )

        model_ref = self.source["slit_models"][0]
        final_cube = cube_lib.setup_final_model(
            model_ref,
            self.cube_wcs,
            self.spaxel_flux,
            self.spaxel_iflux,
            self.spaxel_var,
            self.spaxel_dq,
            None,
            self.cube_wcs.rot_angle,
            self.coord_system,
            "drizzle",
            "drizzle",
            0,
            0,
            0,
            0,
            None,
        )

        # self.spaxel_exposure_counter,
        # bunit_data,
        # bunit_err,
        # )

        return final_cube

    def map_source_pixels_to_output_frame(self, slit):
        """
        Map detector pixels from a NIRSpec slit model to the target output sky frame.

        Transforms detector (x, y) coordinates through the slit model WCS to derive
        tangent-plane standard coordinates (xi, eta), central wavelengths, wavelength
        bin widths (Delta lambda), and 4-corner sky projections for valid science pixels.

        Parameters
        ----------
        slit : `jwst.datamodels.SlitModel`
            The input NIRSpec slit data model containing detector pixel data,
            flags, and WCS information.

        Returns
        -------
        x_array : ndarray of float
            Flattened array of 1D pixel x-coordinates on the detector for good pixels.
        y_array : ndarray of float
            Flattened array of 1D pixel y-coordinates on the detector for good pixels.
        flux_array : ndarray of float
            Flux values for good pixels.
        error_array : ndarray of float
            Uncertainty/error values for good pixels.
        dq_array : ndarray of uint32
            Data quality flags associated with the good pixels.
        var_rnoise_array : ndarray of float
            Read-noise variance for good pixels.
        xi_array : ndarray of float
            Standard tangent plane xi coordinates corresponding to
            the pixel centers.
        eta_array : ndarray of float
            Standard tangent plane eta coordinates  corresponding to
            the pixel centers.
        lam_array : ndarray of float
            Central wavelengths (lambda) for good pixels.
        lam_delta_array : ndarray of float
            Estimated spectral dispersion step (Delta lambda) spanned by each pixel.
        corner : list of ndarray
            List containing xi_1, eta_1, xi_2, eta_2, xi_3, eta_3, xi_4, eta_4]
            representing the projected standard coordinates of all 4 pixel corners.
        wave_pixel : ndarray of float
            2D grid matching detector dimensions populated with wavelengths for good
            pixels and NaN for bad pixels.
        dwave_pixel : ndarray of float
            2D grid matching detector dimensions populated with Delta lambda for good
            pixels and NaN for bad pixels.

        Notes
        -----
        - Pixel filtering rejects invalid values (`NaN`/Inf) and pixels flagged with
        `DO_NOT_USE` or `NON_SCIENCE` in the Data Quality (DQ) array.
        - If no valid pixels are found, empty/NaN-filled structures are returned with a warning.
        """
        # Loop over each slit - map from the detector to the sky (output frame)
        #   1. for each x,y in a slit map to --> sky and store the following
        #   a. x,y, flux, error, dq at x,y position on input data
        #   b. ra, dec, wavelength
        #   c. corner in tangent plane of pixel (xi, eta)
        #   d. delta lambda (estimation of wavelength bin)  at x,y
        #

        # map pixels to sky, forming the ra,dec, corners
        flux = slit.data
        err = slit.err
        dq = slit.dq
        var_rnoise = slit.var_rnoise
        wcs = slit.meta.wcs
        bbox = wcs.bounding_box

        grid = wcstools.grid_from_bounding_box(bbox)
        x = grid[0, :, :]
        y = grid[1, :, :]

        # set up transforms  -> detector -> slit frame --> world
        d2s = wcs.get_transform("detector", "slit_frame")
        s2w = wcs.get_transform("slit_frame", "world")

        # On the detector the x axis is dispersion, y axis is along the slit
        # For a given pixel find the edges in the along slit direction.
        # The detector to slit transforms result across values always being zero:center of the slit.

        # We need to find transformed pixel corners the sky coordinate
        #    (really the tangent to the sky coordinate).
        # We want to map the pixel edges slit, this will give us the pixel edges in the along slit
        # direction.
        across, along1, _ = d2s(x, y - 0.4999)
        across, along2, _ = d2s(x, y + 0.4999)

        # We want a delta lambda for each pixel. Here we use full 'wcs'. We add a 1/2 pixel in the
        # x dimension (dispersion direction).
        _, _, lam1 = wcs(x - 0.4999, y)
        _, _, lam2 = wcs(x + 0.4999, y)
        d_pixel = np.abs(lam2 - lam1)

        # Map pixel to sky (this is the center of the pixel used in drizzle).
        across_slit, along_slit, lam_pixel = d2s(x, y)
        ra_pixel, dec_pixel, lam = s2w(across_slit, along_slit, lam_pixel)
        xi, eta = cube_lib.radec2std(
            self.cube_wcs.crval1, self.cube_wcs.crval2, ra_pixel, dec_pixel, self.cube_wcs.rot_angle
        )

        # slit width is in fractions of slit - see jwst.assign_wcs.util.compute_footprint_nrs_slit
        # a full slit width = 0.5 in these equations
        # the slit_frac parameter allows the user to shrink the slit width in a similar fashion
        # as pix_frac shrinks the 'drop size' of the pixel in classic drizzle.

        frac_slit = self.slit_frac * 0.5

        # We need to add the slit_width to the across coordinates to find the pixel corners.
        # bottom left corner: (where right/left is across the slit)
        # across is an array of zeros
        ra1, dec1, _ = s2w(across - frac_slit, along1, lam1)
        # bottom right corner:
        ra2, dec2, _ = s2w(across + frac_slit, along1, lam1)
        # upper right corner:
        ra3, dec3, _ = s2w(across + frac_slit, along2, lam2)
        # upper left corner:
        ra4, dec4, _ = s2w(across - frac_slit, along2, lam2)

        # selecting good data. This is probably overkill. But it accounts for pixel edges that fall
        # off the valid region.

        bad1 = np.bitwise_and(dq, datamodels.dqflags.pixel["DO_NOT_USE"]).astype(bool)
        bad2 = np.bitwise_and(dq, datamodels.dqflags.pixel["NON_SCIENCE"]).astype(bool)
        good = np.where(
            ~bad1
            & ~bad2
            & np.isfinite(flux)
            & np.isfinite(lam)
            & np.isfinite(lam1)
            & np.isfinite(lam2)
            & np.isfinite(ra1)
            & np.isfinite(ra2)
            & np.isfinite(ra3)
            & np.isfinite(ra4)
        )
        ngood = len(good[0])
        wave_pixel = np.full_like(lam_pixel, np.nan, dtype=np.float64)
        wave_pixel[good] = lam_pixel[good]

        dwave_pixel = np.full_like(d_pixel, np.nan, dtype=np.float64)
        dwave_pixel[good] = d_pixel[good]

        if ngood == 0:
            log.warning("WARNING: slit contains no valid data")
            return (
                np.array([]),
                np.array([]),
                np.array([]),
                np.array([]),
                np.array([], dtype=np.uint32),
                np.array([]),
                np.array([]),
                np.array([]),
                np.array([]),
                np.array([]),
                [],
                wave_pixel,
                dwave_pixel,
            )

        # Extract good pixels.
        x_array = x[good].flatten()
        y_array = y[good].flatten()
        flux_array = flux[good].flatten()
        error_array = err[good].flatten()
        dq_array = dq[good].flatten()
        var_rnoise_array = var_rnoise[good].flatten()  # used in readnoise weighting
        xi_array = xi[good].flatten()
        eta_array = eta[good].flatten()
        lam_array = lam[good].flatten()
        lam_delta_array = d_pixel[good].flatten()

        ra1 = ra1[good]
        ra2 = ra2[good]
        ra3 = ra3[good]
        ra4 = ra4[good]

        dec1 = dec1[good]
        dec2 = dec2[good]
        dec3 = dec3[good]
        dec4 = dec4[good]

        xi1, eta1 = cube_lib.radec2std(
            self.cube_wcs.crval1, self.cube_wcs.crval2, ra1, dec1, self.cube_wcs.rot_angle
        )
        xi2, eta2 = cube_lib.radec2std(
            self.cube_wcs.crval1, self.cube_wcs.crval2, ra2, dec2, self.cube_wcs.rot_angle
        )
        xi3, eta3 = cube_lib.radec2std(
            self.cube_wcs.crval1, self.cube_wcs.crval2, ra3, dec3, self.cube_wcs.rot_angle
        )
        xi4, eta4 = cube_lib.radec2std(
            self.cube_wcs.crval1, self.cube_wcs.crval2, ra4, dec4, self.cube_wcs.rot_angle
        )
        corner = [xi1, eta1, xi2, eta2, xi3, eta3, xi4, eta4]

        return (
            x_array,
            y_array,
            flux_array,
            error_array,
            dq_array,
            var_rnoise_array,
            xi_array,
            eta_array,
            lam_array,
            lam_delta_array,
            corner,
            wave_pixel,
            dwave_pixel,
        )

    def setup_final_model(
        self,
        model_ref,
        cube_wcs,
        rot_angle,
        spaxel_flux,
        spaxel_iflux,
        spaxel_var,
        spaxel_dq,
        spaxel_exposure_counter,
        bunit,
        bunit_err,
    ):
        """
        Reshape 1D spaxel arrays into a 3D spectral cube and populate the IFUCubeModel.

        Takes 1D flattened arrays resulting from the spatial and spectral mapping stage,
        reshapes them into 3D (wave, dec, ra) data arrays, applies quality masking,
        populates meta WCS header keywords, attaches a transformation WCS object,
        and assigns photometric unit attributes.

        Parameters
        ----------
        cube_wcs : tuple
            A tuple of 12 WCS parameters containing:
            (naxis1, naxis2, naxis3, cdelt1, cdelt2, cdelt3,
             crpix1, crpix2, crpix3, crval1, crval2, crval3).
        rot_angle : float or None
            Rotation angle of the output cube in degrees. If `None`, defaults to `0.0`.
        spaxel_flux : numpy.ndarray
            1D array containing cumulative flux per spaxel.
        spaxel_iflux : numpy.ndarray
            1D array containing weight or count mapping information per spaxel.
        spaxel_var : numpy.ndarray
            1D array containing cumulative variance per spaxel.
        spaxel_dq : numpy.ndarray
            1D array containing data quality flags per spaxel.
        spaxel_exposure_counter : numpy.ndarray
            1D array containing exposure counts per spaxel, mapped to the final DQ array.
        bunit : str
            Physical unit of the main science data extension (e.g., `'MJy/sr'`).
        bunit_err : str
            Physical unit of the error extension.

        Returns
        -------
        cube_model : jwst.datamodels.IFUCubeModel
            Fully initialized 3D IFU Data Model with populated WCS metadata,
            photometry info, wavetable, and AST-based WCS object.

        Notes
        -----
        Output data is an IFUCube model - we may want to define MSACubeModel - but this model will
        work for now.
        """
        index = np.where(spaxel_iflux == 0)
        spaxel_flux[index] = np.nan
        spaxel_var[index] = np.nan

        # Rearrange values from 1d vectors into 3d cubes

        (
            naxis1,
            naxis2,
            naxis3,
            cdelt1,
            cdelt2,
            cdelt3,
            crpix1,
            crpix2,
            crpix3,
            crval1,
            crval2,
            crval3,
        ) = cube_wcs

        flux = spaxel_flux.reshape((naxis3, naxis2, naxis1))
        wmap = spaxel_iflux.reshape((naxis3, naxis2, naxis1))
        var = spaxel_var.reshape((naxis3, naxis2, naxis1))
        dq = spaxel_dq.reshape((naxis3, naxis2, naxis1))
        exp_count = spaxel_exposure_counter.reshape((naxis3, naxis2, naxis1))
        dnu = np.where((dq & dqflags.pixel["DO_NOT_USE"]) != 0)
        flux[dnu] = np.nan
        var[dnu] = np.nan
        err = np.sqrt(var)

        if self.linear_wave:
            pixels = np.arange(naxis3)

            # Calculate wavelengths
            # We add 1 to 'pixels' to convert 0-based Python indexing to 1-based FITS indexing
            wavelength_table = crval3 + (pixels + 1 - crpix3) * cdelt3
            wave = np.asarray(wavelength_table, dtype=np.float32)
            num = len(wave)
            alldata = np.array([(wave[None].T,)], dtype=[("wavelength", "<f4", (num, 1))])
            # always write the wavetable
            cube_model = datamodels.IFUCubeModel(
                data=flux, dq=dq, err=var, weightmap=wmap, wavetable=alldata, expcount=exp_count
            )

        else:
            wave = np.asarray(self.wavelength_table, dtype=np.float32)
            num = len(wave)
            alldata = np.array([(wave[None].T,)], dtype=[("wavelength", "<f4", (num, 1))])

            cube_model = datamodels.IFUCubeModel(
                data=flux, dq=dq, err=err, weightmap=wmap, wavetable=alldata, expcount=exp_count
            )

        # Write the information in header we want to keep track of.
        # TODO - this is not working as expected.
        cube_model.update(model_ref)  # to get the history in the header

        cube_model.meta.wcsinfo.crval1 = crval1
        cube_model.meta.wcsinfo.crval2 = crval2
        cube_model.meta.wcsinfo.crpix1 = crpix1
        cube_model.meta.wcsinfo.crpix2 = crpix2

        cube_model.meta.wcsinfo.cdelt1 = cdelt1 / 3600.0
        cube_model.meta.wcsinfo.cdelt2 = cdelt2 / 3600.0
        # Now that we've got a pixel scale, set photometric area keywords
        cube_model.meta.photometry.pixelarea_arcsecsq = cdelt1 * cdelt2
        cube_model.meta.photometry.pixelarea_steradians = (
            cube_model.meta.photometry.pixelarea_arcsecsq * 2.3504e-11
        )

        cube_model.meta.wcsinfo.crval3 = crval3
        cube_model.meta.wcsinfo.cdelt3 = cdelt3
        cube_model.meta.wcsinfo.ctype3 = "WAVE"
        cube_model.meta.wcsinfo.crpix3 = crpix3

        cube_model.meta.wcsinfo.ctype1 = "RA---TAN"
        cube_model.meta.wcsinfo.ctype2 = "DEC--TAN"
        cube_model.meta.wcsinfo.cunit1 = "deg"
        cube_model.meta.wcsinfo.cunit2 = "deg"

        cube_model.meta.wcsinfo.cunit3 = "um"
        cube_model.meta.wcsinfo.wcsaxes = 3
        cube_model.meta.wcsinfo.pc1_1 = -1
        cube_model.meta.wcsinfo.pc1_2 = 0
        cube_model.meta.wcsinfo.pc1_3 = 0

        cube_model.meta.wcsinfo.pc2_1 = 0
        cube_model.meta.wcsinfo.pc2_2 = 1
        cube_model.meta.wcsinfo.pc2_3 = 0

        cube_model.meta.wcsinfo.pc3_1 = 0
        cube_model.meta.wcsinfo.pc3_2 = 0
        cube_model.meta.wcsinfo.pc3_3 = 1

        if rot_angle is None:
            rot_angle = 0.0

        cube_model.meta.wcsinfo.pc1_1 = -np.cos(rot_angle * np.pi / 180.0)
        cube_model.meta.wcsinfo.pc1_2 = np.sin(rot_angle * np.pi / 180.0)
        cube_model.meta.wcsinfo.pc2_1 = np.sin(rot_angle * np.pi / 180.0)
        cube_model.meta.wcsinfo.pc2_2 = np.cos(rot_angle * np.pi / 180.0)

        cube_model.meta.ifu.flux_extension = "SCI"
        cube_model.meta.ifu.error_extension = "ERR"
        cube_model.meta.ifu.error_type = "ERR"
        cube_model.meta.ifu.dq_extension = "DQ"
        cube_model.meta.ifu.exp_count_extension = "EXPCOUNT"
        cube_model.meta.ifu.weighting = self.weighting

        cube_model.meta.bunit_data = bunit
        cube_model.meta.bunit_err = bunit_err

        cube_model.meta.ifu.slit_frac = self.slit_frac
        # stick in values of 0, otherwise it is NaN and
        # fits file can not be written because these
        # values are defined in ifucube.schema.yaml
        cube_model.meta.ifu.weight_power = 0
        cube_model.meta.ifu.roi_wave = 0
        cube_model.meta.ifu.roi_spatial = 0

        # set WCS information
        wcsobj = pointing.create_fitswcs(cube_model)
        cube_model.meta.wcs = wcsobj
        cube_model.meta.wcs.bounding_box = ((0, naxis1 - 1), (0, naxis2 - 1), (0, naxis3 - 1))
        return cube_model
