"""Work horse routines used for building MSA spectra cubes."""

import copy
import logging
import math

import numpy as np
from gwcs import wcstools
from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels import dqflags

from jwst.assign_wcs import pointing
from jwst.assign_wcs.util import compute_footprint_nrs_slit, wrap_ra
from jwst.cube_build import coord
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
        self.weighting = pars.get("weighting")
        self.wavelength_table = pars.get("wavelength_table")
        self.suffix = pars.get("suffix")

        self.cdelt1 = None
        self.cdelt2 = None
        self.cdelt3 = None
        self.crpix1 = None
        self.crpix2 = None
        self.crpix3 = None
        self.crval1 = None
        self.crval2 = None
        self.crval3 = None
        self.naxis1 = None
        self.naxis2 = None
        self.naxis3 = None
        self.cdelt3_normal = None
        self.rot_angle = None  # rotation angle between Ra-Dec and MSA local instrument plane

        self.a_min = 0
        self.a_max = 0
        self.b_min = 0
        self.b_max = 0
        self.lambda_min = 0
        self.lambda_max = 0
        self.xcoord = None
        self.ycoord = None
        self.zcoord = None

    def find_footprint(self):
        """
        Set up the WCS of the MSA cube.

        Loop over every slitlet contained in the MultiExposureModel and find the WCS
        of the output cube that contains all the data.

        Returns
        -------
        corner_a : list
            List of corner ra values of slit
        corner_b : list
            List of corner dec values of slit
        final_lam_min : float
            Minimum wavelength
        final_lam_max : float
            Maximum wavelength
        rot_angle : float
            Rotation angle

        Notes
        -----
        For NIRSPEC the units along/across slice dimension are meters.

        Module sets the wcs class parameters.
        If the coordinate system is ``skyalign``/``msaalign``, then the min and max of
        RA (degrees), dec (degrees), and lambda (microns) are returned
        for internal calculations.
        """
        self.cdelt1 = self.scalexy
        self.cdelt2 = self.scalexy
        if self.linear_wave:
            self.cdelt3 = self.scalew

        # Define the rotation angle
        # If coord_system = msaalign then the angle is between the ra-dec and alpha beta
        # coord system using the first input model. Use first file in first band to set up
        # rotation angle.
        # Compute the rotation angle between local MSA system  and RA-DEC

        corner_a = []
        corner_b = []

        lambda_min = []
        lambda_max = []

        num = len(self.source["slit_models"])
        log.info(f"Number of slits:  {num}")

        for i in range(num):
            slit = self.source["slit_models"][i]
            wcs = copy.deepcopy(slit.meta.wcs)  # TODO do we need a deep copy
            bbox = wcs.bounding_box

            grid = wcstools.grid_from_bounding_box(bbox)
            ra, dec, lam = np.array(wcs(*grid))  # use to get lambda regions

            dq = slit.dq
            bad1 = np.bitwise_and(dq, dqflags.pixel["DO_NOT_USE"]).astype(bool)
            bad2 = np.bitwise_and(dq, dqflags.pixel["NON_SCIENCE"]).astype(bool)
            good_data = np.where(~bad1 & ~bad2)
            lam = lam[good_data]

            if (
                i == 0
            ):  # only use first slit to figure out the rotation between msa along slit and sky.
                # we want the same rotation applied to all the data.
                if self.coord_system == "msaalign":
                    # set up transforms
                    s2w = wcs.get_transform("slit_frame", "world")
                    lam_med = np.nanmedian(lam)
                    temp_ra1, temp_dec1, lam_temp = s2w(0, 0, lam_med)
                    temp_ra2, temp_dec2, lam_temp = s2w(0, 0.005, lam_med)

                    dra, ddec = (
                        (temp_ra2 - temp_ra1) * np.cos(temp_dec1 * np.pi / 180.0),
                        (temp_dec2 - temp_dec1),
                    )
                    rot_angle = 90 + np.arctan2(dra, ddec) * 180.0 / np.pi
                    log.info(f"Rotation angle between the msa plane and th sky  {rot_angle}")
                else:
                    if self.cube_pa is not None:
                        rot_angle = self.cube_pa
                    else:
                        rot_angle = None

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

                corner_b.append(cb3)
                corner_a.append(ca3)

                corner_b.append(cb4)
                corner_a.append(ca4)

                lambda_min.append(lmin)
                lambda_max.append(lmax)

        final_lam_min = min(lambda_min)
        final_lam_max = max(lambda_max)

        return (corner_a, corner_b, final_lam_min, final_lam_max, rot_angle)

    def set_slit_wcs(self, corner_ra, corner_dec, lambda_min, lambda_max, rot_angle):
        """
        Set up the 3D spatial/spectral World Coordinate System (WCS).

        This method calculates the spatial field of view (RA/Dec mapped to tangent plane
        coordinates xi and eta) and sets the corresponding FITS WCS standard attributes
        (CRVAL, CRPIX, NAXIS) along with the spaxel center arrays for spatial (x, y)
        and spectral (z) axes.

        Parameters
        ----------
        corner_ra : array-like of float
            Right Ascension values at the footprint corners, in degrees.
        corner_dec : array-like of float
            Declination values at the footprint corners, in degrees.
        lambda_min : float
            Minimum wavelength threshold for the spectral axis.
        lambda_max : float
            Maximum wavelength threshold for the spectral axis.
        rot_angle : float
            Position or rotation angle for the sky-to-tangent plane transformation, in degrees.

        Notes
        -----
        Updates the following instance attributes:
        - `crval1`, `crval2`, `crval3` : Reference point coordinates.
        - `crpix1`, `crpix2`, `crpix3` : Reference point pixel indices (1-indexed).
        - `naxis1`, `naxis2`, `naxis3` : Cube pixel dimensions along x, y, and wavelength.
        - `a_min`, `a_max`, `b_min`, `b_max` : Min/Max spatial bounds in standard coordinates.
        - `xcoord`, `ycoord`, `zcoord` : 1D arrays of spaxel center positions.
        - `lambda_min`, `lambda_max` : Adjusted bounds for linear wavelength axis.
        - `cdelt3_normal` : Array of wavelength bin step sizes for flux normalization.

        - Uses `coord.radec2std` to transform RA/Dec into projected standard coordinates xi, eta.
        - Ensures the spatial projection grid is centered at (xi=0, eta=0) by padding `naxis1`
        and `naxis2` to odd dimensions.
        - If `self.ra_center` or `self.dec_center` are pre-defined, they override the mean
        calculated tangent reference points.
        - Supports both linear (`self.linear_wave = True`) and tabular non-linear spectral grids.
        """
        ra_min = min(corner_ra)
        ra_max = max(corner_ra)
        dec_min = min(corner_dec)
        dec_max = max(corner_dec)
        dec_ave = (dec_min + dec_max) / 2.0

        # we can not average ra values because of the convergence
        # of hour angles.
        ravalues = np.zeros(2)
        ravalues[0] = ra_min
        ravalues[1] = ra_max
        ra_ave = np.nanmean(wrap_ra(ravalues))

        self.crval1 = ra_ave
        self.crval2 = dec_ave

        if self.ra_center is not None:
            self.crval1 = self.ra_center
        if self.dec_center is not None:
            self.crval2 = self.dec_center

        # find the 4 corners - tangent plane through crval1, crval2
        xi_corner = []
        eta_corner = []
        num = len(corner_ra)

        for i in range(num):
            xi, eta = coord.radec2std(
                self.crval1, self.crval2, corner_ra[i], corner_dec[i], rot_angle
            )
            xi_corner.append(xi)
            eta_corner.append(eta)

        xi_min = min(xi_corner)
        xi_max = max(xi_corner)
        eta_min = min(eta_corner)
        eta_max = max(eta_corner)

        xilimit = max(np.abs(xi_min), np.abs(xi_max))
        etalimit = max(np.abs(eta_min), np.abs(eta_max))

        # Adjust the cube to be slightly bigger I want the cube to be centered
        # on xi and eta = 0  If no data is mapped to spaxel the value will be NAN.
        # xi, eta is set by  by crval1, crval2.
        # number min

        na = math.ceil(xilimit / self.cdelt1) + 1
        nb = math.ceil(etalimit / self.cdelt2) + 1

        if self.nspax_x is not None:
            na = self.nspax_x
        if self.nspax_y is not None:
            nb = self.nspax_y

        # full range of xi and eta values
        xi_min = 0.0 - (na * self.cdelt1) - (self.cdelt1 / 2.0)
        xi_max = (na * self.cdelt1) + (self.cdelt1 / 2.0)

        eta_min = 0.0 - (nb * self.cdelt2) - (self.cdelt2 / 2.0)
        eta_max = (nb * self.cdelt2) + (self.cdelt2 / 2.0)

        self.a_min = xi_min
        self.a_max = xi_max
        self.b_min = eta_min
        self.b_max = eta_max
        # find the CRPIX1 CRPIX2 - xi and eta centered at 0,0
        self.crpix1 = float(na) + 1.0
        self.crpix2 = float(nb) + 1.0

        self.naxis1 = na * 2 + 1  # + 1 to set the center  = 0
        self.naxis2 = nb * 2 + 1

        # center of spaxels
        self.xcoord = np.zeros(self.naxis1)
        xstart = xi_min + self.cdelt1 / 2.0
        self.xcoord = np.arange(
            start=xstart, stop=xstart + self.naxis1 * self.cdelt1, step=self.cdelt1
        )

        self.ycoord = np.zeros(self.naxis2)
        ystart = eta_min + self.cdelt2 / 2.0
        self.ycoord = np.arange(
            start=ystart, stop=ystart + self.naxis2 * self.cdelt2, step=self.cdelt2
        )

        # depending on the naxis and cdelt values the x,ycoord can have 1 more element than naxis.
        # Clean up arrays dropping extra values at the end.

        # TODO - this could be removed just left for testing.
        # A check that at crpix values the xi,eta values are close to zero
        test1 = self.xcoord[int(self.crpix1 - 1)]
        test2 = self.ycoord[int(self.crpix2 - 1)]

        # print('xi coordinate must be (close to ) zero', test1)  #TODO remove after testing
        # print('eta coordinate must be (close to)  zero', test2)
        if np.abs(test1) > 5e-6 or np.abs(test2) > 5e-6:
            log.warning(" Warning xi, eta grid center is not at zero")

        self.xcoord = self.xcoord[0 : self.naxis1]
        self.ycoord = self.ycoord[0 : self.naxis2]

        # now wavelength axis
        if self.linear_wave:
            range_lambda = lambda_max - lambda_min
            self.naxis3 = int(math.ceil(range_lambda / self.cdelt3))
            # adjust max based on integer value of naxis3
            lambda_max = lambda_min + (self.naxis3) * self.cdelt3
            self.lambda_min = lambda_min
            self.lambda_max = lambda_max
            self.zcoord = np.zeros(self.naxis3)
            # CRPIX3 for FITS is 1 (center of first pixel)
            # CRVAL3 then is lambda_min + self.cdelt3/ 2.0, which is also zcoord[0]
            # Note that these are all values at the center of a spaxel
            self.crval3 = lambda_min + self.cdelt3 / 2.0
            self.crpix3 = 1.0
            zstart = lambda_min + self.cdelt3 / 2.0
            self.zcoord = zstart + np.arange(self.naxis3) * self.cdelt3
        else:
            self.zcoord = np.asarray(self.wavelength_table, dtype=np.float64)
            self.naxis3 = len(self.zcoord)
            self.crval3 = self.zcoord[0]
            self.crpix3 = 1.0
            # Set default cdelt3 step estimate for non-linear tables
            self.cdelt3 = np.mean(np.diff(self.zcoord)) if self.naxis3 > 1 else 0.0

        # Vectorized setup for cdelt3_normal normalization array
        self.cdelt3_normal = np.empty(self.naxis3)
        self.cdelt3_normal[:-1] = np.diff(self.zcoord)
        self.cdelt3_normal[-1] = self.cdelt3_normal[-2]

    def print_geometry(self):
        """Print to log the general properties of the size of the MSA cube."""
        log.info("MSA Cube Geometry:")

        log.info("axis#  Naxis  CRPIX    CRVAL       CDELT(arcsec)  Min & Max (xi, eta arcsec)")
        log.info(
            "Axis 1 %5d  %5.2f %12.8f %12.8f %12.8f %12.8f",
            self.naxis1,
            self.crpix1,
            self.crval1,
            self.cdelt1,
            self.a_min,
            self.a_max,
        )
        log.info(
            "Axis 2 %5d  %5.2f %12.8f %12.8f %12.8f %12.8f",
            self.naxis2,
            self.crpix2,
            self.crval2,
            self.cdelt2,
            self.b_min,
            self.b_max,
        )
        if self.linear_wave:
            log.info("axis#  Naxis  CRPIX    CRVAL      CDELT(microns)  Min & Max (microns)")
            log.info(
                "Axis 3 %5d  %5.2f %12.8f %12.8f %12.8f %12.8f",
                self.naxis3,
                self.crpix3,
                self.crval3,
                self.cdelt3,
                self.lambda_min,
                self.lambda_max,
            )

        if not self.linear_wave:
            log.info("Non-linear wavelength dimension; CDELT3 variable")
            log.info("axis#  Naxis  CRPIX    CRVAL     Min & Max (microns)")
            log.info(
                "Axis 3 %5d  %5.2f %12.8f %12.8f %12.8f",
                self.naxis3,
                self.crpix3,
                self.crval3,
                self.wavelength_table[0],
                self.wavelength_table[self.naxis3 - 1],
            )

        if self.rot_angle is not None:
            log.info("Rotation angle between RA-Dec and Slicer-Plane %12.8f", self.rot_angle)

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
        total_num = self.naxis1 * self.naxis2 * self.naxis3

        self.spaxel_flux = np.zeros(total_num, dtype=np.float64)
        self.spaxel_weight = np.zeros(total_num, dtype=np.float64)
        self.spaxel_var = np.zeros(total_num, dtype=np.float64)
        self.spaxel_iflux = np.zeros(total_num, dtype=np.float64)
        self.spaxel_dq = np.zeros(total_num, dtype=np.uint32)
        self.spaxel_exposure_counter = np.zeros(total_num, dtype=np.uint32)

        # we need bunit in the out header - is there another way to do this ? #TODO
        bunit_data = self.source["bunit"]
        bunit_err = self.source["bunit_err"]
        nxyplane = self.naxis1 * self.naxis2

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
                cdelt3_mean = np.nanmean(self.cdelt3_normal)
                xi1, eta1, xi2, eta2, xi3, eta3, xi4, eta4 = corner

            if self.weighting == "readnoise":
                weight_type = 1
            else:
                weight_type = 0

            linear = 1
            dummy = x_array.copy()  # place holder needed for cube_build - valid for IFU data

            result = cube_wrapper_driz(
                weight_type,
                self.xcoord,
                self.ycoord,
                self.zcoord,
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
                self.cdelt3_normal,
                self.cdelt1,
                self.cdelt2,
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

        self.find_spaxel_flux(
            self.spaxel_iflux, self.spaxel_flux, self.spaxel_weight, self.spaxel_var
        )
        cube_wcs = (
            self.naxis1,
            self.naxis2,
            self.naxis3,
            self.cdelt1,
            self.cdelt2,
            self.cdelt3,
            self.crpix1,
            self.crpix2,
            self.crpix3,
            self.crval1,
            self.crval2,
            self.crval3,
        )

        model_ref = self.source["slit_models"][0]
        final_cube = self.setup_final_model(
            model_ref,
            cube_wcs,
            self.rot_angle,
            self.spaxel_flux,
            self.spaxel_iflux,
            self.spaxel_var,
            self.spaxel_dq,
            self.spaxel_exposure_counter,
            bunit_data,
            bunit_err,
        )

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

        x_array = []
        y_array = []
        flux_array = []
        error_array = []
        dq_array = []
        xi_array = []
        eta_array = []
        lam_array = []
        corner = []
        lam_delta_array = []
        var_rnoise_array = []

        # map pixels to sky, forming the ra,dec, corners
        flux = slit.data
        err = slit.err
        dq = slit.dq
        var_rnoise = slit.var_rnoise
        wcs = copy.deepcopy(slit.meta.wcs)

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
        w_pixel = lam_pixel
        ra_pixel, dec_pixel, lam = s2w(across_slit, along_slit, lam_pixel)
        xi, eta = coord.radec2std(self.crval1, self.crval2, ra_pixel, dec_pixel, self.rot_angle)

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
        valid1 = np.isfinite(flux)
        good = np.where(~np.isnan(lam))
        valid2 = np.isfinite(lam)
        valid3 = np.isfinite(ra1)
        valid4 = np.isfinite(ra2)
        valid31 = np.isfinite(ra3)
        valid41 = np.isfinite(ra4)

        valid5 = np.isfinite(lam1)
        valid6 = np.isfinite(lam2)
        bad1 = np.bitwise_and(dq, datamodels.dqflags.pixel["DO_NOT_USE"]).astype(bool)
        bad2 = np.bitwise_and(dq, datamodels.dqflags.pixel["NON_SCIENCE"]).astype(bool)
        good = np.where(
            ~bad1 & ~bad2 & valid1 & valid2 & valid3 & valid4 & valid5 & valid6 & valid31 & valid41
        )
        ngood = len(good[0])

        wave_pixel = np.full_like(w_pixel, np.nan, dtype=np.float64)
        wave_pixel[good] = w_pixel[good]

        dwave_pixel = np.full_like(d_pixel, np.nan, dtype=np.float64)
        dwave_pixel[good] = d_pixel[good]

        if ngood == 0:
            log.warning("WARNING: slit contains no valid data")
            return (
                x_array,
                y_array,
                np.array(flux_array),
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

        flux = flux[good]
        dq = dq[good]
        x = x[good]
        y = y[good]
        err = err[good]
        var_rnoise = var_rnoise[good]  # used in readnoise weighting

        along_slit = along_slit[good]
        across_slit = across_slit[good]
        xi = xi[good]
        eta = eta[good]
        lam = lam[good]
        ra1 = ra1[good]
        ra2 = ra2[good]
        ra3 = ra3[good]
        ra4 = ra4[good]

        dec1 = dec1[good]
        dec2 = dec2[good]
        dec3 = dec3[good]
        dec4 = dec4[good]

        lam1 = lam1[good]
        lam2 = lam2[good]
        delta_lam = np.abs(lam2 - lam1)

        # store information
        x_array.append(x)
        y_array.append(y)
        flux_array.append(flux)
        error_array.append(err)
        dq_array.append(dq)
        var_rnoise_array.append(var_rnoise)

        xi_array.append(xi)
        eta_array.append(eta)
        lam_array.append(lam)
        lam_delta_array.append(delta_lam)

        xi1, eta1 = coord.radec2std(self.crval1, self.crval2, ra1, dec1, self.rot_angle)
        xi2, eta2 = coord.radec2std(self.crval1, self.crval2, ra2, dec2, self.rot_angle)
        xi3, eta3 = coord.radec2std(self.crval1, self.crval2, ra3, dec3, self.rot_angle)
        xi4, eta4 = coord.radec2std(self.crval1, self.crval2, ra4, dec4, self.rot_angle)
        corner = [xi1, eta1, xi2, eta2, xi3, eta3, xi4, eta4]

        x_array = np.array(x_array).flatten()
        y_array = np.array(y_array).flatten()
        flux_array = np.array(flux_array).flatten()
        error_array = np.array(error_array).flatten()
        dq_array = np.array(dq_array).flatten()
        var_rnoise_array = np.array(var_rnoise_array).flatten()
        xi_array = np.array(xi_array).flatten()
        lam_array = np.array(lam_array).flatten()
        eta_array = np.array(eta_array).flatten()
        lam_delta_array = np.array(lam_delta_array).flatten()

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

    def find_spaxel_flux(self, spaxel_iflux, spaxel_flux, spaxel_weight, spaxel_var):
        """
        Calculate the final normalized flux and variance for each spaxel.

        Normalizes the weighted flux and variance arrays after mapping overlap
        between the detector array and the sky. Spaxels with no contributing points
        (where spaxel_iflux == 0) are ignored to prevent division-by-zero errors.

        Parameters
        ----------
        spaxel_iflux : numpy.ndarray
            Array containing the counts/number of detector elements contributing to each spaxel.
        spaxel_flux : numpy.ndarray
             Array of weighted cumulative flux values. Modified in-place.
        spaxel_weight : numpy.ndarray
             Array of total weights accumulated for each spaxel.Modified in-place.
        spaxel_var : numpy.ndarray
             Array of cumulative variance values. Modified in-place.

        Notes
        -----
            Modifies `spaxel_flux` and `spaxel_var` in-place.
        """
        # after overlapping is performed between the detector mapped array and sky find the
        # final spaxel values for flux and variance.

        # Don't apply any normalization if no points contributed to a spaxel
        # (i.e., don't divide by zero)

        good = spaxel_iflux > 0

        # Normalize the weighted sum of pixel fluxes by the sum of the weights
        spaxel_flux[good] = spaxel_flux[good] / spaxel_weight[good]
        # Normalize the variance by  the weights
        spaxel_var[good] = spaxel_var[good] / (spaxel_weight[good] * spaxel_weight[good])

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
