"""Work horse routines used for building ifu spectra cubes."""

import logging
import math
import warnings

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.stats import circmean
from gwcs import wcstools
from gwcs.utils import to_index
from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels import dqflags

from jwst.assign_wcs import nirspec, pointing
from jwst.assign_wcs.util import compute_footprint_nrs_slit
from jwst.assign_wcs.util import wrap_ra
from jwst.cube_build import coord
from jwst.cube_build.cube_match_sky_driz import cube_wrapper_driz  # c extension
from jwst.model_blender import blendmeta

log = logging.getLogger(__name__)

__all__ = ["MSACubeData", "IncorrectInputError", "IncorrectParameterError"]


class MSACubeData:
    """
    Combine MSA slitlet data onto a regular grid.

    Parameters
    ----------
    input_models : `~jwst.datamodels.MultiExposureModel`


    **pars: dict
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
        self.rot_angle = None  # rotation angle between Ra-Dec and IFU local instrument plane

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

        Notes
        -----
        For NIRSPEC the units along/across slice dimension are meters.

        If the coordinate system is ``skyalign``/``ifualign``, then the min and max of
        RA (degrees), dec (degrees), and lambda (microns) are returned
        for internal calculations.
        """

        print('find_footprint', self.scalexy)
        self.cdelt1 = self.scalexy
        self.cdelt2 = self.scalexy
        if self.linear_wave:
            self.cdelt3 = self.scalew

        # Define the rotation angle
        # If coord_system = msaalign then the angle is between the ra-dec and alpha beta
        # coord system using the first input model. Use first file in first band to set up
        # rotation angle.
        # Compute the rotation angle between local IFU system  and RA-DEC

        corner_a = []
        corner_b = []
    
        lambda_min = []
        lambda_max = []

        num = len(self.source['slitname'])

        for i in range(num):
            slit = self.source['slit_models'][i]
            wcs  = self.source['wcs'][i]
            bbox = wcs.bounding_box

            grid = wcstools.grid_from_bounding_box(bbox)
            ra,dec,lam = np.array(wcs(*grid)) # use to get lambda regions

            dq = slit.dq
            bad1 = np.bitwise_and(dq, dqflags.pixel['DO_NOT_USE']).astype(bool)
            bad2 = np.bitwise_and(dq, dqflags.pixel['NON_SCIENCE']).astype(bool)
            good_data = np.where(~bad1 & ~bad2)
            lam = lam[good_data]

            if i == 0: # only use first slit to figure out the rotation between msa along slit and sky.
                # we want the same rotation applied to all the data. 
                if self.coord_system == 'msaalign':
                    # set up transforms
                    s2w = wcs.get_transform('slit_frame', 'world')
                    lam_med = np.nanmedian(lam)
                    temp_ra1, temp_dec1, lam_temp = s2w(0, 0, lam_med)
                    temp_ra2, temp_dec2, lam_temp = s2w(0, 0.005, lam_med)
                
                    dra, ddec = (temp_ra2 - temp_ra1) * np.cos(temp_dec1 * np.pi / 180.0), (temp_dec2 - temp_dec1)
                    rot_angle = 90 + np.arctan2(dra, ddec) * 180. / np.pi
                    print('Rotation angle between msa plane and sky: ',rot_angle)
                else:
                    if self.cube_pa is not None:
                        rot_angle = self.cube_pa
                    else:
                        rot_angle = None
                

        if len(good_data[0]) > 0:
            lmin  = np.nanmin(lam)
            lmax  = np.nanmax(lam)
            footprint, lam = compute_footprint_nrs_slit(slit)
            ca1 = float(footprint[0][0])
            cb1 = float(footprint[0][1])
            ca2 = float(footprint[1][0])
            cb2 = float(footprint[1][1])
            ca3 = float(footprint[2][0])
            cb3 = float(footprint[2][1])
            ca4 = float(footprint[3][0])
            cb4 = float(footprint[3][1])

            print(ca1,ca2,ca3,ca4)
            print(cb1,cb2,cb3,cb4)
            
            rac = [ca1, ca2, ca3, ca4, ca1]
            decc= [cb1, cb2, cb3, cb4, cb1]

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

        final_a_min = min(corner_a)
        final_a_max = max(corner_a)
        final_b_min = min(corner_b)
        final_b_max = max(corner_b)
        final_lam_min = min(lambda_min)
        final_lam_max = max(lambda_max)
        dec_ave = (final_b_min + final_b_max)/2
        rarange = (final_a_max - final_a_min) * math.cos(math.radians(dec_ave))
        decrange = final_b_max - final_b_min
        print(' ra min max, ra coverage(arc seconds)', final_a_min, final_a_max, rarange*3600.0)
        print(' dec min max, dec coverage(arc seconds)', final_b_min, final_b_max, decrange*3600.0)
        print(' lambda min max (microns)', final_lam_min, final_lam_max,)

        return (corner_a, corner_b, final_lam_min, final_lam_max, rot_angle)



    def set_slit_wcs(self, corner_ra, corner_dec, lambda_min, lambda_max, rot_angle):
                 #Cdelt1, cdelt2, cdelt3, rot_angle, align_msa,
                 #ra_center, dec_center, nspax_x, nspax_y):
        """ 
        Based on the ra,dec and wavelength footprint set up the size
        of the cube in the tangent plane projected coordinate system.
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

        print('in set slit wcs', self.cdelt1, self.cdelt2)
        if self.ra_center is not None:
            self.crval1 = self.ra_center
        if self.dec_center is not None:
            self.crval2 = self.dec_center
        
        # find the 4 corners - tangent plane through crval1, crval2
        xi_corner = []
        eta_corner = []
        num = len(corner_ra)

        for i in range(num):
            xi, eta = coord.radec2std(self.crval1, self.crval2,
                            corner_ra[i], corner_dec[i], rot_angle)
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

        print('xi limit', xilimit[0], etalimit[0])
        na = math.ceil(xilimit[0] / self.cdelt1) + 1
        nb = math.ceil(etalimit[0] / self.cdelt2) + 1

        if self.nspax_x is not None:
            na =self. nspax_x
        if self.nspax_y is not None:
            nb = self.nspax_y
    
        # full range of xi and eta values 
        xi_min = 0.0 - (na * self.cdelt1) - (self.cdelt1 / 2.0)
        xi_max = (na * self.cdelt1) + (self.cdelt1 / 2.0)

        eta_min = 0.0 - (nb * self.cdelt2) - (self.cdelt2 / 2.0)
        eta_max = (nb * self.cdelt2) + (self.cdelt2 / 2.0)

        # find the CRPIX1 CRPIX2 - xi and eta centered at 0,0
        self.crpix1 = float(na) + 1.0
        self.crpix2 = float(nb) + 1.0
    
        self.naxis1 = na * 2 + 1  # + 1 to set the center  = 0 
        self.naxis2 = nb * 2 + 1

        # center of spaxels
        self.xcoord = np.zeros(self.naxis1)
        xstart = xi_min + self.cdelt1 / 2.0
        self.xcoord = np.arange(start=xstart, stop=xstart + self.naxis1 * self.cdelt1,
                                step=self.cdelt1)
    
        self.ycoord = np.zeros(self.naxis2)
        ystart = eta_min + self.cdelt2 / 2.0
        self.ycoord = np.arange(start=ystart, stop=ystart + self.naxis2 * self.cdelt2,
                                step=self.cdelt2)

        # print(' Coordinates to tangent plane:')
        print('Xi ', xi_min, np.max(self.xcoord) + self.cdelt1/2.0)
        print('Eta ', eta_min, np.max(self.ycoord) + self.cdelt2/2.0)
        # depending on the naxis and cdelt values the x,ycoord can have 1 more element than naxis.
        # Clean up arrays dropping extra values at the end.

        # Remove later or keep
        test1 = self.xcoord[int(self.crpix1-1)]
        test2 = self.ycoord[int(self.crpix2-1)]

        print('xi coordinate must be (close to ) zero', test1)
        print('eta coordinate must be (close to)  zero', test2)
        if np.abs(test1) > 1e-6 or np.abs(test2) > 1e-6:
            print(' ERROR xi, eta grid center is not at zero')
        
        self.xcoord = self.xcoord[0:self.naxis1]
        self.ycoord = self.ycoord[0:self.naxis2]

        xv, yv = np.meshgrid(self.xcoord, self.ycoord)
        xcenters = xv.flatten()
        ycenters = yv.flatten()

        # now wavelength axis
        range_lambda = lambda_max - lambda_min
        self.naxis3 = int(math.ceil(range_lambda / self.cdelt3))
                
        # adjust max based on integer value of naxis3
        lambda_max = lambda_min + (self.naxis3) *self.cdelt3

        self.zcoord = np.zeros(self.naxis3)
        # CRPIX3 for FITS is 1 (center of first pixel)
        # CRVAL3 then is lambda_min + self.cdelt3/ 2.0, which is also zcoord[0]
        # Note that these are all values at the center of a spaxel
        self.crval3 = lambda_min + self.cdelt3 / 2.0
        self.crpix3 = 1.0
        zstart = lambda_min + self.cdelt3 / 2.0
        self.zcoord = np.arange(start=zstart, stop=lambda_max, step=self.cdelt3)
        self.zcoord = self.zcoord[0:self.naxis3]
    
        # set up the cdelt3_normal normalizing array used
        self.cdelt3_normal = np.zeros(self.naxis3)
        for j in range(self.naxis3 - 1):
            self.cdelt3_normal[j] = self.zcoord[j + 1] - self.zcoord[j]

        self.cdelt3_normal[self.naxis3 - 1] = self.cdelt3_normal[self.naxis3 - 2]
    

    # _______________________________________________________________________
    def set_geometry(self, corner_a, corner_b, lambda_min, lambda_max):
        """
        Set up the WCS of the cube in the tangent plane.

        Parameters
        ----------
        corner_a : ndarray
            Array of RA corners of the footprint of all input data.
        corner_b : ndarray
            Array of Dec corners of the footprint of all input data.
        lambda_min : float
            Minimum wavelength value of the data.
        lambda_max : float
            Maximum wavelength value of the data.
        """
        ra_min = np.min(corner_a)
        ra_max = np.max(corner_a)

        dec_min = np.min(corner_b)
        dec_max = np.max(corner_b)
        dec_ave = (dec_min + dec_max) / 2.0

        # we can not average ra values because of the convergence
        # of hour angles.
        ravalues = np.zeros(2)
        ravalues[0] = ra_min
        ravalues[1] = ra_max

        # astropy circmean assumes angles are in radians,
        # we have angles in degrees
        ra_ave = circmean(ravalues * u.deg).value % 360

        if self.ra_center is not None:
            self.crval1 = self.ra_center
        else:
            self.crval1 = ra_ave

        if self.dec_center is not None:
            self.crval2 = self.dec_center
        else:
            self.crval2 = dec_ave

        rot_angle = self.rot_angle

        # find the 4 corners
        xi_corner = []
        eta_corner = []
        num = len(corner_a)

        for i in range(num):
            xi, eta = coord.radec2std(self.crval1, self.crval2, corner_a[i], corner_b[i], rot_angle)
            xi_corner.append(xi)
            eta_corner.append(eta)
        xi_min = min(xi_corner)
        xi_max = max(xi_corner)
        eta_min = min(eta_corner)
        eta_max = max(eta_corner)

        # find the CRPIX1 CRPIX2 - xi and eta centered at 0,0
        # to find location of center abs of min values is how many pixels
        # we want a symmetric cube centered on xi,eta = 0
        xilimit = max(np.abs(xi_min), np.abs(xi_max))
        etalimit = max(np.abs(eta_min), np.abs(eta_max))

        na = math.ceil((xilimit / self.cdelt1).item()) + 1
        nb = math.ceil((etalimit / self.cdelt2).item()) + 1

        # if the user set the nspax_x or nspax_y then redefine na, nb
        # it is assumed that both values are ODD numbers
        # We want the central pixel to be the tangent point with na/nb pixels on either
        # side of central pixel.
        if self.nspax_x is not None:
            na = int(self.nspax_x / 2)
        if self.nspax_y is not None:
            nb = int(self.nspax_y / 2)

        xi_min = 0.0 - (na * self.cdelt1) - (self.cdelt1 / 2.0)
        xi_max = (na * self.cdelt1) + (self.cdelt1 / 2.0)

        eta_min = 0.0 - (nb * self.cdelt2) - (self.cdelt2 / 2.0)
        eta_max = (nb * self.cdelt2) + (self.cdelt2 / 2.0)

        self.crpix1 = float(na) + 1.0
        self.crpix2 = float(nb) + 1.0

        self.naxis1 = na * 2 + 1
        self.naxis2 = nb * 2 + 1

        self.a_min = xi_min
        self.a_max = xi_max
        self.b_min = eta_min
        self.b_max = eta_max
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
        self.xcoord = self.xcoord[0 : self.naxis1]
        self.ycoord = self.ycoord[0 : self.naxis2]

        xv, yv = np.meshgrid(self.xcoord, self.ycoord)
        self.xcenters = xv.flatten()
        self.ycenters = yv.flatten()

        # set up the lambda (z) coordinate of the cube
        self.cdelt3_normal = None
        if self.linear_wave:
            self.lambda_min = lambda_min
            self.lambda_max = lambda_max
            range_lambda = self.lambda_max - self.lambda_min
            self.naxis3 = int(math.ceil(range_lambda / self.cdelt3))

            # adjust max based on integer value of naxis3
            self.lambda_max = self.lambda_min + (self.naxis3) * self.cdelt3

            self.zcoord = np.zeros(self.naxis3)
            # CRPIX3 for FITS is 1 (center of first pixel)
            # CRVAL3 then is lambda_min + self.cdelt3/ 2.0, which is also zcoord[0]
            # Note that these are all values at the center of a spaxel
            self.crval3 = self.lambda_min + self.cdelt3 / 2.0
            self.crpix3 = 1.0
            zstart = self.lambda_min + self.cdelt3 / 2.0
            self.zcoord = np.arange(start=zstart, stop=self.lambda_max, step=self.cdelt3)
            self.zcoord = self.zcoord[0 : self.naxis3]
        else:
            self.naxis3 = len(self.wavelength_table)
            self.zcoord = np.asarray(self.wavelength_table)
            self.crval3 = self.wavelength_table[0]
            self.crpix3 = 1.0
        # set up the cdelt3_normal normalizing array used
        cdelt3_normal = np.zeros(self.naxis3)
        for j in range(self.naxis3 - 1):
            cdelt3_normal[j] = self.zcoord[j + 1] - self.zcoord[j]

        cdelt3_normal[self.naxis3 - 1] = cdelt3_normal[self.naxis3 - 2]
        self.cdelt3_normal = cdelt3_normal

    # _______________________________________________________________________
    def print_cube_geometry(self):
        """Print to log the general properties of the size of the IFU cube."""
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


    # ________________________________________________________________________________
    def build_msacube(self):
        """
        Create an MSA cube.

        1. Loop over every band contained in the IFU cube and read in the data
           associated with the band.
        2. :meth:`map_detector_to_outputframe`: Maps the detector data to the cube
           output coordinate system.
        3. For each mapped detector pixel, find the IFU cube spaxel located in the region of
           interest. There are two different routines to do this step,
           both of which use a C extension
           to combine the detector fluxes that fall within a region of influence
           from the spaxel center:

           a. ``src/cube_match_sky*.c``: This routine uses the modified
              Shepard method to determine the weighting function, which weights the detector
              fluxes based on the distance between the detector center and spaxel center.
           b. ``src/cube_match_internal.c`` is only for single exposure, single band cubes, and
              the IFU cube in created in the detector plane. The weighting function is based on
              the overlap of between the detector pixel and spaxel. This method is simplified
              to determine the overlap in the along slice-wavelength plane.

        4. :meth:`find_spaxel_flux`: find the final flux associated with each spaxel
        5. :meth:`setup_final_ifucube_model`

        Returns
        -------
        result : `~stdatamodels.jwst.datamodels.IFUCubeModel`
            An IFU cube of combined IFU image data.
        """

        total_num = self.naxis1 * self.naxis2 * self.naxis3

        self.spaxel_flux = np.zeros(total_num, dtype=np.float64)
        self.spaxel_weight = np.zeros(total_num, dtype=np.float64)
        self.spaxel_var = np.zeros(total_num, dtype=np.float64)
        self.spaxel_iflux = np.zeros(total_num, dtype=np.float64)
        self.spaxel_dq = np.zeros(total_num, dtype=np.uint32)

        spaxel_err = np.zeros(total_num, dtype=np.float64)
        spaxel_exposure_counter = np.zeros(total_num, dtype=np.uint32)
        
        nxyplane = self.naxis1 * self.naxis2

        if self.spaxel_z == -1 and self.spaxel_x == -1 and self.spaxel_y == -1:
            debug_cube_index = -1

        elif self.spaxel_z < 0 or self.spaxel_x < 0 or self.spaxel_y < 0:
            log.info("Incorrect input for Debug Spaxel values. Counting starts at 0")
            debug_cube_index = -1
            log.info(f"{self.spaxel_z} {self.spaxel_x}  {self.spaxel_y}")
        else:
            spaxel_z = self.spaxel_z
            spaxel_x = self.spaxel_x
            spaxel_y = self.spaxel_y
            debug_cube_index = spaxel_z * (nxyplane) + spaxel_y * self.naxis1 + spaxel_x
            log.info(
                f"Printing debug information for cube spaxel:  {spaxel_x} {spaxel_y} {spaxel_z}"
            )

       
            num = self.source['number_slits']
            slit_name = []
            shutter_id = []
            wavelengths = []
            primary_hdu = fits.PrimaryHDU()
            hdul = fits.HDUList([primary_hdu])

            j = 0
            #num = 1 #For testing
            all_results = []
            for i in range(num):
                print(f"\rworking on slit {i} out of {num}", end="", flush=True)                     
                slit = self.source['slit_models'][i]

                results = map_source_pixels_to_output_frame(sname, i,  source, 
                                                    crval1, crval2, rot_angle,
                                                    ra_offset, dec_offset,
                                                    slit_frac,
                                                    check_dq_flags,
                                                    bartol, 
                                                    plot)

            
                x_array, y_array, flux, err, dq, var_rnoise, \
                    xi_array, eta_array, wave, dwave, corner, \
                    wave_pixel, delta_pixel = results


                npt = flux.size
                # Check that the slit has valid pixels 
                if npt > 0:
                    cdelt3_mean = np.nanmean(cdelt3_norm)
                    xi1, eta1, xi2, eta2, xi3, eta3, xi4, eta4 = corner
            
                    bunit_data = source['bunit']
                    bunit_err = source['bunit_err']

            if run_c_driz:
                linear = 1
                instrument = 1
                flag_dq_plane = 0
                start_region = 0
                end_region = 0
                overlap_partial = 0
                overlap_full = 0
                slice_no = 0 
                dummy = x_array.copy()                

                result = cube_wrapper_driz(instrument, flag_dq_plane,
                                           start_region, end_region,
                                           overlap_partial, overlap_full,
                                           xcoord, ycoord, zcoord,
                                           xi_array, eta_array, wave , flux, err, dummy,
                                           xi1, eta1, xi2, eta2, xi3, eta3, xi4, eta4,
                                           dwave, cdelt3_normal,
                                           cdelt1, cdelt2, cdelt3_mean,linear,
                                           x_array, y_array, 
                                           debug_spaxel_index)

                this_spaxel_flux, this_spaxel_weight, this_spaxel_var, this_spaxel_iflux, this_spaxel_dq = result
                spaxel_flux = spaxel_flux + np.asarray(this_spaxel_flux, np.float64)
                spaxel_weight = spaxel_weight + np.asarray(this_spaxel_weight, np.float64)
                spaxel_var = spaxel_var + np.asarray(this_spaxel_var, np.float64)
                spaxel_iflux = spaxel_iflux + np.asarray(this_spaxel_iflux, np.float64)
                spaxel_dq.astype(np.uint)
                spaxel_dq = np.bitwise_or(spaxel_dq, this_spaxel_dq)
                result = None
                del result
                del this_spaxel_flux, this_spaxel_weight, this_spaxel_var, this_spaxel_iflux, this_spaxel_dq
                
            else: 
                match_driz_msa(i, xcoord, ycoord, zcoord, x_array, y_array,
                               wave , flux, err, dq, var_rnoise, 
                               xi1, eta1, xi2, eta2, xi3, eta3, xi4, eta4,
                               dwave, cdelt3_normal, cdelt1, cdelt2, 
                               naxis1, naxis2, naxis3, total_num, npt,
                               spaxel_flux, spaxel_weight, spaxel_var, spaxel_iflux,
                               spaxel_exposure_counter,
                               weight_type, debug_spaxel_index, all_results)

   #  print('data written to overlaps.npy', all_results)
    matrix_data = np.array(all_results)
    # Save as a binary file
    file_overlap = 'overlaps/' + str(sname[0]) + '_' + str(slit_frac) + '_' + 'overlaps.npy'
    print(file_overlap)
    np.save(file_overlap, matrix_data)
    
    find_spaxel_flux(run_c_driz, spaxel_iflux, spaxel_flux, spaxel_weight, spaxel_var)
    cube_wcs = (naxis1, naxis2, naxis3,  cdelt1, cdelt2, cdelt3, crpix1, crpix2, crpix3, crval1, crval2, crval3)


    print('Send rot angle to setup_final_model', rot_angle)
        
    final_cube = setup_final_model(sname, cube_wcs,rot_angle, 
                                   spaxel_flux, spaxel_iflux, spaxel_var, spaxel_dq,
                                   spaxel_exposure_counter, 
                                   zcoord, bunit_data, bunit_err, output_name)

    return final_cube


        number_bands = len(self.list_par1)
        k = 0
        for ib in range(number_bands):
            this_par1 = self.list_par1[ib]
            this_par2 = self.list_par2[ib]
            for input_model in self.master_table.FileMap[self.instrument][this_par1][this_par2]:
                # loop over the files that cover the spectral range the cube is for

                self.input_models_this_cube.append(input_model)
                # set up input_model to be first file used to copy in basic header info
                # to ifucube meta data
                if ib == 0 and k == 0:
                    input_model_ref = input_model

                log.debug(f"Working on Band defined by: {this_par1} {this_par2}")
                input_frame = input_model.meta.wcs.available_frames[0]
                if self.interpolation in ["pointcloud", "drizzle"]:
                    pixelresult = self.map_detector_to_outputframe(this_par1, input_model)

                    (
                        coord1,
                        coord2,
                        corner_coord,
                        wave,
                        dwave,
                        flux,
                        err,
                        slice_no,
                        rois_pixel,
                        roiw_pixel,
                        weight_pixel,
                        softrad_pixel,
                        scalerad_pixel,
                        x_det,
                        y_det,
                    ) = pixelresult

                    # by default flag the dq plane based on the FOV of the detector projected to sky
                    flag_dq_plane = 1
                    if self.skip_dqflagging:
                        flag_dq_plane = 0

                    # check that there is valid data returned
                    # If all the data is flagged as DO_NOT_USE - not common- then log warning
                    build_cube = True
                    if wave is None:
                        log.warning(f"No valid data found on file {input_model.meta.filename}")
                        flag_dq_plane = 0
                        build_cube = False

                    # C extension setup
                    start_region = 0
                    end_region = 0


                    instrument = 1

                    result = None
                    weight_type = 0  # default to emsm instead of msm

                    if self.weighting == "drizzle" and build_cube:
                        cdelt3_mean = np.nanmean(self.cdelt3_normal)
                        xi1, eta1, xi2, eta2, xi3, eta3, xi4, eta4 = corner_coord
                        linear = 0
                        if self.linear_wave:
                            linear = 1
                        if debug_cube_index >= 0:
                            log.info(f"Input filename: {input_model.meta.filename}")
                        result = cube_wrapper_driz(
                            instrument,
                            flag_dq_plane,
                            start_region,
                            end_region,
                            self.overlap_partial,
                            self.overlap_full,
                            self.xcoord,
                            self.ycoord,
                            self.zcoord,
                            coord1,
                            coord2,
                            wave,
                            flux,
                            err,
                            slice_no,
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
                            x_det,
                            y_det,
                            debug_cube_index,
                        )

                        spaxel_flux, spaxel_weight, spaxel_var, spaxel_iflux, spaxel_dq = result
                        self.spaxel_flux = self.spaxel_flux + np.asarray(spaxel_flux, np.float64)
                        self.spaxel_weight = self.spaxel_weight + np.asarray(
                            spaxel_weight, np.float64
                        )
                        self.spaxel_var = self.spaxel_var + np.asarray(spaxel_var, np.float64)
                        self.spaxel_iflux = self.spaxel_iflux + np.asarray(spaxel_iflux, np.float64)
                        spaxel_dq.astype(np.uint)
                        self.spaxel_dq = np.bitwise_or(self.spaxel_dq, spaxel_dq)
                        result = None
                        del result
                        del spaxel_flux, spaxel_weight, spaxel_var, spaxel_iflux, spaxel_dq



        self.find_spaxel_flux()
        self.set_final_dq_flags()

        # shove Flux and iflux in the  final IFU cube
        result = self.setup_final_ifucube_model(input_model_ref)
        return result

    # ________________________________________________________________________________
    def map_detector_to_outputframe(self, this_par1, input_model):
        """
        Loop over a file and map the detector pixels to the output cube.

        The output frame is on the sky (RA-Dec).
        Return the coordinates of all the detector pixel in the output frame.
        In addition, an array of pixel fluxes and weighing parameters are
        determined. The pixel flux and weighing parameters are used later in
        the process to find the final flux of a cube spaxel based on the pixel
        fluxes and pixel weighing parameters that fall within the ROI of
        spaxel center

        Parameters
        ----------
        this_par1 : str
           For MIRI, this is the channel number (1, 2, 3, or 4);
           needed for MIRI to distinguish which channel on the detector we have.
           For NIRSPEC, this is the grating name.
        input_model : `~stdatamodels.jwst.datamodels.IFUImageModel`
           Input IFU image model to combine.

        Returns
        -------
        coord1 : ndarray
           Coordinate for ``axis1`` in output cube for mapped pixel
        coord2 : ndarray
           Coordinate for ``axis2`` in output cube for mapped pixel
        wave : ndarray
           Wavelength associated with ```coord1, coord2``
        flux : ndarray
           Flux associated with ``coord1, coord2``
        err : ndarray
           Error associated with ``coord1, coord2``
        rois_det : float
           Spatial ROI size to use
        roiw_det : ndarray
           Spectral ROI size associated with ``coord1,coord2``
        weight_det : ndarray
            Weighting parameter associated with ``coord1,coord2``
        softrad_det : ndarray
            Softrad parameter associated with ``coord1,coord2``
        """
        # initialize alpha_det and beta_det to None. These are filled in
        # if the instrument is MIRI and the weighting is miripsf

        wave_all = None
        slice_no_all = None  # Slice number
        dwave_all = None
        corner_coord_all = None

        # initialize values to be returned to None
        dwave = None
        corner_coord = None
        coord1 = None
        coord2 = None
        wave = None
        flux = None
        err = None
        slice_no = None
        rois_det = None
        roiw_det = None
        weight_det = None
        softrad_det = None
        scalerad_det = None
        x_det = None
        y_det = None
        offsets = self.offsets

        if self.instrument == "MIRI":
            sky_result = self.map_miri_pixel_to_sky(input_model, this_par1, offsets)
            (x, y, ra, dec, wave_all, slice_no_all, dwave_all, corner_coord_all) = sky_result

        elif self.instrument == "NIRSPEC":
            sky_result = self.map_nirspec_pixel_to_sky(input_model, offsets)
            (x, y, ra, dec, wave_all, slice_no_all, dwave_all, corner_coord_all) = sky_result

        # ______________________________________________________________________________
        # The following is for both MIRI and NIRSPEC

        flux_all = input_model.data[y, x]
        err_all = input_model.err[y, x]
        dq_all = input_model.dq[y, x]
        valid2 = np.isfinite(flux_all)

        x_all = x
        y_all = y
        # Preselect only data within a given wavelength range
        # This range is defined to include all pixels for which the chosen wavelength region
        # of interest would have them fall within one of the cube spectral planes
        # Note that the cube lambda refer to the spaxel midpoints, so we must account for both
        # the spaxel width and the ROI size

        if self.linear_wave:
            min_wave_tolerance = self.zcoord[0] - self.roiw
            max_wave_tolerance = self.zcoord[-1] + self.roiw
        else:
            min_wave_tolerance = self.zcoord[0] - np.max(self.roiw_table)
            max_wave_tolerance = self.zcoord[-1] + np.max(self.roiw_table)
        if self.interpolation == "drizzle":
            dmax = np.nanmax(dwave_all)

            if self.linear_wave:
                min_wave_tolerance = self.zcoord[0] - (self.cdelt3 + dmax)
                max_wave_tolerance = self.zcoord[-1] + (self.cdelt3 + dmax)
            else:
                min_wave_tolerance = self.zcoord[0] - self.cdelt3_normal[0]
                max_wave_tolerance = self.zcoord[-1] + self.cdelt3_normal[-1]
        valid_min = np.where(wave_all >= min_wave_tolerance)
        not_mapped_low = wave_all.size - len(valid_min[0])
        valid_max = np.where(wave_all <= max_wave_tolerance)
        not_mapped_high = wave_all.size - len(valid_max[0])
        if not_mapped_low > 0:
            log.info(
                "# of detector pixels not mapped to output plane: "
                f"{not_mapped_low} with wavelength below {min_wave_tolerance}"
            )

        if not_mapped_high > 0:
            log.info(
                "# of detector pixels not mapped to output plane: "
                f"{not_mapped_high} with wavelength above {max_wave_tolerance}"
            )

        # using the DQFlags from the input_image find pixels that should be excluded
        # from the cube mapping

        valid3 = np.logical_and(wave_all >= min_wave_tolerance, wave_all <= max_wave_tolerance)

        # find the location of good data

        bad1 = np.bitwise_and(dq_all, dqflags.pixel["DO_NOT_USE"]).astype(bool)
        bad2 = np.bitwise_and(dq_all, dqflags.pixel["NON_SCIENCE"]).astype(bool)
        good_data = np.where(~bad1 & ~bad2 & valid2 & valid3)

        num_good = len(good_data[0])
        if num_good == 0:  # This can occur if all the pixels on the detector are marked DO_NOT_USE.
            log.warning(f"No valid pixels found on detector {input_model.meta.filename}")
            return (
                coord1,
                coord2,
                corner_coord,
                wave,
                dwave,
                flux,
                err,
                slice_no,
                rois_det,
                roiw_det,
                weight_det,
                softrad_det,
                scalerad_det,
                x_det,
                y_det,
            )

        # good data holds the location of pixels we want to map to cube
        # define variables as numpy arrays (numba needs this defined)
        flux_all_good = flux_all[good_data]
        good_shape = flux_all_good.shape
        flux = np.zeros(good_shape, dtype=np.float64)
        err = np.zeros(good_shape, dtype=np.float64)
        coord1 = np.zeros(good_shape, dtype=np.float64)
        coord2 = np.zeros(good_shape, dtype=np.float64)
        wave = np.zeros(good_shape, dtype=np.float64)
        slice_no = np.zeros(good_shape)

        flux[:] = flux_all_good
        err[:] = err_all[good_data]
        wave[:] = wave_all[good_data]
        slice_no[:] = slice_no_all[good_data]
        x_det = x_all[good_data]
        y_det = y_all[good_data]

        log.debug(f"After removing pixels min and max wave: {np.min(wave)} {np.max(wave)}")

        # based on the wavelength define the sroi, wroi, weight_power and
        # softrad to use in matching detector to spaxel values
        rois_det = np.zeros(wave.shape)
        roiw_det = np.zeros(wave.shape)
        weight_det = np.zeros(wave.shape)
        softrad_det = np.zeros(wave.shape)
        scalerad_det = np.zeros(wave.shape)

        # ________________________________________________________________________
        # if working with msm or emsm need roi sizes and other parameters defined:
        if self.weighting == "msm" or self.weighting == "emsm":
            if self.linear_wave:
                rois_det[:] = self.rois
                roiw_det[:] = self.roiw
                weight_det[:] = self.weight_power
                softrad_det[:] = self.soft_rad
                scalerad_det[:] = self.scalerad
            else:
                # for each wavelength find the closest point in the self.wavelength_table
                for iw, w in enumerate(wave):
                    self.find_closest_wave(
                        iw,
                        w,
                        self.wavelength_table,
                        self.rois_table,
                        self.roiw_table,
                        self.softrad_table,
                        self.weight_power_table,
                        self.scalerad_table,
                        rois_det,
                        roiw_det,
                        softrad_det,
                        weight_det,
                        scalerad_det,
                    )

        ra_use = ra[good_data]
        dec_use = dec[good_data]
        coord1, coord2 = coord.radec2std(self.crval1, self.crval2, ra_use, dec_use, self.rot_angle)

        if self.interpolation == "drizzle":
            dwave = np.zeros(good_shape)
            dwave[:] = dwave_all[good_data]
            ra1 = corner_coord_all[0]
            dec1 = corner_coord_all[1]
            ra2 = corner_coord_all[2]
            dec2 = corner_coord_all[3]
            ra3 = corner_coord_all[4]
            dec3 = corner_coord_all[5]
            ra4 = corner_coord_all[6]
            dec4 = corner_coord_all[7]
            ra1 = ra1[good_data]
            dec1 = dec1[good_data]
            ra2 = ra2[good_data]
            dec2 = dec2[good_data]
            ra3 = ra3[good_data]
            dec3 = dec3[good_data]
            ra4 = ra4[good_data]
            dec4 = dec4[good_data]

            xi1, eta1 = coord.radec2std(self.crval1, self.crval2, ra1, dec1, self.rot_angle)

            xi2, eta2 = coord.radec2std(self.crval1, self.crval2, ra2, dec2, self.rot_angle)
            xi3, eta3 = coord.radec2std(self.crval1, self.crval2, ra3, dec3, self.rot_angle)
            xi4, eta4 = coord.radec2std(self.crval1, self.crval2, ra4, dec4, self.rot_angle)

            corner_coord = [xi1, eta1, xi2, eta2, xi3, eta3, xi4, eta4]
        return (
            coord1,
            coord2,
            corner_coord,
            wave,
            dwave,
            flux,
            err,
            slice_no,
            rois_det,
            roiw_det,
            weight_det,
            softrad_det,
            scalerad_det,
            x_det,
            y_det,
        )

    # ______________________________________________________________________
    def map_miri_pixel_to_sky(self, input_model, this_par1, offsets):
        """
        Loop over a MIRI model and map the detector pixels to the output cube.

        The output frame is on the sky (RA-Dec).
        Return the coordinates of all the detector pixel in the output frame
        for every valid input pixel from the IFU image model.

        Parameters
        ----------
        input_model : `~stdatamodels.jwst.datamodels.IFUImageModel`
           Input IFU image model to combine
        this_par1 : str
           For MIRI, this is the channel number.
           needed for MIRI to distinguish which channel on the detector we have.
           For NIRSpec, this is the grating name.
        offsets : dict
           Optional dictionary of RA and Dec offsets to apply

        Returns
        -------
        x, y : float
            The pixel values on the detector
        ra, dec : float
            Detector values mapped to sky
        wave : float
            Wavelength corresponding to pixel
        slice_no : int
            Slice number of the pixel
        dwave : float
            Delta wavelength covered by pixel
        corner_coord : tuple
            The corners of the pixel mapped to ``ra,dec``
        """
        wave = None
        slice_no = None  # Slice number
        dwave = None
        corner_coord = None


        # find the slice number of each pixel and fill in slice_det
        ysize, xsize = input_model.data.shape
        if input_model.hasattr("regions"):
            # expected for oversampled data
            slice_det = input_model.regions
        else:
            # find the slice number of each pixel and fill in slice_det
            slice_det = np.zeros((ysize, xsize), dtype=int)
            det2ab_transform = input_model.meta.wcs.get_transform("detector", "alpha_beta")
            mapper = det2ab_transform.label_mapper.mapper
            start_region = self.instrument_info.get_start_slice(this_par1)
            end_region = self.instrument_info.get_end_slice(this_par1)
            regions = list(range(start_region, end_region + 1))
            for i in regions:
                ys, xs = (mapper == i).nonzero()
                xind = to_index(xs)
                yind = to_index(ys)
                xind = np.ndarray.flatten(xind)
                yind = np.ndarray.flatten(yind)
                slice_det[yind, xind] = i

        # define the x,y detector values of channel to be mapped to desired coordinate system
        xstart, xend = self.instrument_info.get_miri_slice_endpts(this_par1)
        xc_start, xc_end = cube_build_wcs_util.miri_slice_limit_coords(
            input_model.meta.wcs, xstart, xend
        )
        y, x = np.mgrid[:ysize, xc_start:xc_end]
        y = np.reshape(y, y.size)
        x = np.reshape(x, x.size)

        # if self.coord_system == 'skyalign' or self.coord_system == 'ifualign':
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", "invalid value", RuntimeWarning)
            ra, dec, wave = input_model.meta.wcs(x, y)

        # offset the central pixel
        if offsets is not None:
            ra, dec = self.offset_coord(ra, dec, raoffset, decoffset)

        valid1 = ~np.isnan(ra)
        ra = ra[valid1]
        dec = dec[valid1]
        wave = wave[valid1]
        x = x[valid1]
        y = y[valid1]

        xind = to_index(x)
        yind = to_index(y)
        xind = np.ndarray.flatten(xind)
        yind = np.ndarray.flatten(yind)
        slice_no = slice_det[yind, xind]
        input_frame = input_model.meta.wcs.available_frames[0]
        if self.interpolation == "drizzle":
            # Delta wavelengths
            _, _, wave1 = input_model.meta.wcs(x, y - 0.4999)
            _, _, wave2 = input_model.meta.wcs(x, y + 0.4999)
            dwave = np.abs(wave1 - wave2)

            # Pixel corners
            pixfrac = 1.0
            alpha1, beta, _ = input_model.meta.wcs.transform(
                input_frame, "alpha_beta", x - 0.4999 * pixfrac, y
            )
            alpha2, _, _ = input_model.meta.wcs.transform(
                input_frame, "alpha_beta", x + 0.4999 * pixfrac, y
            )
            # Find slice width
            allbetaval = np.unique(beta)
            dbeta = np.abs(allbetaval[1] - allbetaval[0])
            ra1, dec1, _ = input_model.meta.wcs.transform(
                "alpha_beta",
                input_model.meta.wcs.output_frame,
                alpha1,
                beta - dbeta * pixfrac / 2.0,
                wave,
            )
            ra2, dec2, _ = input_model.meta.wcs.transform(
                "alpha_beta",
                input_model.meta.wcs.output_frame,
                alpha1,
                beta + dbeta * pixfrac / 2.0,
                wave,
            )
            ra3, dec3, _ = input_model.meta.wcs.transform(
                "alpha_beta",
                input_model.meta.wcs.output_frame,
                alpha2,
                beta + dbeta * pixfrac / 2.0,
                wave,
            )
            ra4, dec4, _ = input_model.meta.wcs.transform(
                "alpha_beta",
                input_model.meta.wcs.output_frame,
                alpha2,
                beta - dbeta * pixfrac / 2.0,
                wave,
            )

            # now offset the pixel corners
            if offsets is not None:
                ra1, dec1 = self.offset_coord(ra1, dec1, raoffset, decoffset)
                ra2, dec2 = self.offset_coord(ra2, dec2, raoffset, decoffset)
                ra3, dec3 = self.offset_coord(ra3, dec3, raoffset, decoffset)
                ra4, dec4 = self.offset_coord(ra4, dec4, raoffset, decoffset)

            corner_coord = [ra1, dec1, ra2, dec2, ra3, dec3, ra4, dec4]

        sky_result = (x, y, ra, dec, wave, slice_no, dwave, corner_coord)
        return sky_result

    # ______________________________________________________________________
    def map_nirspec_pixel_to_sky(self, input_model, offsets):
        """
        Loop over a NIRSpec model and map the detector pixels to the output cube.

        The output frame is on the sky (RA-Dec).
        Return the coordinates of all the detector pixel in the output frame
        for every valid input pixel from the IFU image model.

        Parameters
        ----------
        input_model : `~stdatamodels.jwst.datamodels.IFUImageModel`
            Input IFU image model to combine
        offsets : ndarray
            RA and Dec offsets to apply to each file

        Returns
        -------
        x, y : float
            The pixel values on the detector
        ra, dec : float
            Detector values mapped to sky
        wave : float
            Wavelength corresponding to pixel
        slice_no : int
            Slice number of the pixel
        dwave : float
            Delta wavelength covered by pixel
        corner_coord : tuple
            Rhe corners of the pixel mapped to ``ra,dec```
        """
        # check if we have an ra and dec offset file
        raoffset = 0.0
        decoffset = 0.0
        # pull out ra dec offset if it exists


        # initialize the ra,dec, and wavelength arrays
        # we will loop over slice_nos and fill in values
        # the flag_det will be set when a slice_no pixel is filled in
        # at the end we will use this flag to pull out valid data

        ysize, xsize = input_model.data.shape
        ra_det = np.zeros((ysize, xsize))
        dec_det = np.zeros((ysize, xsize))
        lam_det = np.zeros((ysize, xsize))
        flag_det = np.zeros((ysize, xsize))
        slice_det = np.zeros((ysize, xsize), dtype=int)
        dwave_det = np.zeros((ysize, xsize))
        ra1_det = np.zeros((ysize, xsize))
        ra2_det = np.zeros((ysize, xsize))
        ra3_det = np.zeros((ysize, xsize))
        ra4_det = np.zeros((ysize, xsize))
        dec1_det = np.zeros((ysize, xsize))
        dec2_det = np.zeros((ysize, xsize))
        dec3_det = np.zeros((ysize, xsize))
        dec4_det = np.zeros((ysize, xsize))

        # determine the slice width using slice 1 and 3
        input_frame = input_model.meta.wcs.available_frames[0]
        slice_wcs1 = nirspec.nrs_wcs_set_input(input_model, 0)
        detector2slicer = slice_wcs1.get_transform(input_frame, "slicer")
        mean_x, mean_y = np.mean(slice_wcs1.bounding_box[0]), np.mean(slice_wcs1.bounding_box[1])
        slice_loc1, _, _ = detector2slicer(mean_x, mean_y)

        slice_wcs3 = nirspec.nrs_wcs_set_input(input_model, 2)
        detector2slicer = slice_wcs3.get_transform(input_frame, "slicer")
        mean_x, mean_y = np.mean(slice_wcs3.bounding_box[0]), np.mean(slice_wcs3.bounding_box[1])
        slice_loc3, _, _ = detector2slicer(mean_x, mean_y)

        across_width = abs(slice_loc1 - slice_loc3)

        # for NIRSPEC each file has 30 slices
        # wcs information access separately for each slice
        nslices = 30
        log.info("Mapping each NIRSpec slice to sky for input file: %s", input_model.meta.filename)

        for ii in range(nslices):
            slice_wcs = nirspec.nrs_wcs_set_input(input_model, ii)

            x, y = wcstools.grid_from_bounding_box(slice_wcs.bounding_box)
            ra, dec, lam = slice_wcs(x, y)

            # the slices are curved on detector so a rectangular region returns NaNs
            valid = ~np.isnan(lam)
            x = x[valid]
            y = y[valid]
            ra = ra[valid]
            dec = dec[valid]
            lam = lam[valid]

            if self.interpolation == "drizzle":
                # Delta wavelengths
                _, _, wave1 = slice_wcs(x - 0.4999, y)
                _, _, wave2 = slice_wcs(x + 0.4999, y)
                dwave = np.abs(wave1 - wave2)

                # Pixel corners
                pixfrac = 1.0
                detector2slicer = slice_wcs.get_transform(input_frame, "slicer")
                slicer2world = slice_wcs.get_transform("slicer", slice_wcs.output_frame)
                across1, along1, lam1 = detector2slicer(x, y - 0.49 * pixfrac)
                across2, along2, lam2 = detector2slicer(x, y + 0.49 * pixfrac)

                # Ensure that our ordering wraps around the footprint instead of crossing
                # footprint on a diagonal
                ra1, dec1, _ = slicer2world(across1 - across_width * pixfrac / 2, along1, lam1)
                ra2, dec2, _ = slicer2world(across1 + across_width * pixfrac / 2, along1, lam1)

                ra3, dec3, _ = slicer2world(across2 + across_width * pixfrac / 2, along2, lam2)
                ra4, dec4, _ = slicer2world(across2 - across_width * pixfrac / 2, along2, lam2)

                # near the slice boundaries the corners can become Nan - do not use pixels with
                # Nan corners
                valid1 = np.logical_and(~np.isnan(ra1), ~np.isnan(ra2))
                valid2 = np.logical_and(~np.isnan(ra3), ~np.isnan(ra4))
                final = np.where(np.logical_and(valid1, valid2))

                x = x[final]
                y = y[final]
                ra = ra[final]
                dec = dec[final]
                lam = lam[final]
                ra1 = ra1[final]
                dec1 = dec1[final]
                ra2 = ra2[final]
                dec2 = dec2[final]
                ra3 = ra3[final]
                dec3 = dec3[final]
                ra4 = ra4[final]
                dec4 = dec4[final]
                dwave = dwave[final]

                xind = to_index(x)
                yind = to_index(y)
                xind = np.ndarray.flatten(xind)
                yind = np.ndarray.flatten(yind)
                ra = np.ndarray.flatten(ra)
                dec = np.ndarray.flatten(dec)
                lam = np.ndarray.flatten(lam)
                ra_det[yind, xind] = ra
                dec_det[yind, xind] = dec
                lam_det[yind, xind] = lam
                flag_det[yind, xind] = 1
                slice_det[yind, xind] = ii + 1

                # fill in corner values
                dwave_det[yind, xind] = dwave
                ra1_det[yind, xind] = ra1
                ra2_det[yind, xind] = ra2
                ra3_det[yind, xind] = ra3
                ra4_det[yind, xind] = ra4

                dec1_det[yind, xind] = dec1
                dec2_det[yind, xind] = dec2
                dec3_det[yind, xind] = dec3
                dec4_det[yind, xind] = dec4

            else:  # not drizzling
                xind = to_index(x)
                yind = to_index(y)
                xind = np.ndarray.flatten(xind)
                yind = np.ndarray.flatten(yind)

                ra = np.ndarray.flatten(ra)
                dec = np.ndarray.flatten(dec)
                lam = np.ndarray.flatten(lam)
                ra_det[yind, xind] = ra
                dec_det[yind, xind] = dec
                lam_det[yind, xind] = lam
                flag_det[yind, xind] = 1
                slice_det[yind, xind] = ii + 1

        # after looping over slices  - pull out valid values

        valid_data = np.where(flag_det == 1)
        y, x = valid_data
        wave = lam_det[valid_data]
        slice_no = slice_det[valid_data]
        dwave = dwave_det[valid_data]

        ra = ra_det[valid_data]
        dec = dec_det[valid_data]
        ra1 = ra1_det[valid_data]
        ra2 = ra2_det[valid_data]
        ra3 = ra3_det[valid_data]
        ra4 = ra4_det[valid_data]
        dec1 = dec1_det[valid_data]
        dec2 = dec2_det[valid_data]
        dec3 = dec3_det[valid_data]
        dec4 = dec4_det[valid_data]

        if offsets is not None:
            # central pixel
            ra, dec = self.offset_coord(ra, dec, raoffset, decoffset)

            # pixel corners
            ra1, dec1 = self.offset_coord(ra1, dec1, raoffset, decoffset)
            ra2, dec2 = self.offset_coord(ra2, dec2, raoffset, decoffset)
            ra3, dec3 = self.offset_coord(ra3, dec3, raoffset, decoffset)
            ra4, dec4 = self.offset_coord(ra4, dec4, raoffset, decoffset)

        corner_coord = [ra1, dec1, ra2, dec2, ra3, dec3, ra4, dec4]
        sky_result = (x, y, ra, dec, wave, slice_no, dwave, corner_coord)
        return sky_result

    # ________________________________________________________________________________
    def find_closest_wave(
        self,
        iw,
        w,
        wavelength_table,
        rois_table,
        roiw_table,
        softrad_table,
        weight_power_table,
        scalerad_table,
        rois_det,
        roiw_det,
        softrad_det,
        weight_det,
        scalerad_det,
    ):
        """
        Given a specific wavelength, find the closest value in the ``wavelength_table``.

        Parameters
        ----------
        iw : int
            Index of wavelength array.
        w : float
            Wavelength array of data.
        wavelength_table : ndarray
            Wavelength array read from :ref:`cubepar_reffile`.
        rois_table : ndarray
            ``rois`` array read from :ref:`cubepar_reffile`.
        roiw_table : ndarray
            ``roiw`` array read from :ref:`cubepar_reffile`.
        softrad_table : ndarray
            Softrad array read from :ref:`cubepar_reffile`.
        weight_power_table : ndarray
            Weight power array read from :ref:`cubepar_reffile`.
        scalerad_table : ndarray
            Scalerad array read from :ref:`cubepar_reffile`.
        rois_det : ndarray
            ``rois`` array of detector pixel for the associated wavelength of the pixel.
        roiw_det : ndarray
            ``roiw`` array of detector pixel for the associated wavelength of the pixel.
        softrad_det : ndarray
            Softrad array of detector pixel for the associated wavelength of the pixel.
        weight_det : ndarray
            Weight array of detector pixel for the associated wavelength of the pixel.
        scalerad_det : ndarray
            Scalerad array of detector pixel for the associated wavelength of the pixel.
        """
        ifound = (np.abs(wavelength_table - w)).argmin()
        rois_det[iw] = rois_table[ifound]
        roiw_det[iw] = roiw_table[ifound]
        softrad_det[iw] = softrad_table[ifound]
        weight_det[iw] = weight_power_table[ifound]
        scalerad_det[iw] = scalerad_table[ifound]

    # ________________________________________________________________________________
    def find_spaxel_flux(self):
        """Depending on the interpolation method, find the flux for each spaxel value."""
        # currently these are the same but in the future there could be a difference in
        # how the spaxel flux is determined according to self.interpolation.
        if self.interpolation == "area":
            good = self.spaxel_iflux > 0
            self.spaxel_flux[good] = self.spaxel_flux[good] / self.spaxel_weight[good]
            self.spaxel_var[good] = self.spaxel_var[good] / (
                self.spaxel_weight[good] * self.spaxel_weight[good]
            )
        elif self.interpolation == "pointcloud" or self.interpolation == "drizzle":
            # Don't apply any normalization if no points contributed to a spaxel
            # (i.e., don't divide by zero)
            good = self.spaxel_iflux > 0

            # Normalize the weighted sum of pixel fluxes by the sum of the weights
            self.spaxel_flux[good] = self.spaxel_flux[good] / self.spaxel_weight[good]
            # Normalize the variance by the square of the weights
            self.spaxel_var[good] = self.spaxel_var[good] / (
                self.spaxel_weight[good] * self.spaxel_weight[good]
            )

    # ________________________________________________________________________________
    def set_final_dq_flags(self):
        """
        Set up the final DQ flags.

        These flags include:

        * Good data (0)
        * NON_SCIENCE
        * DO_NOT_USE.
        """
        # An initial set of dq flags was set in overlap_fov_with_spaxel or
        # overlap_slice_with_spaxel. The initial dq dlags are defined in ifu_cube
        # class:
        # self.overlap_partial = 4  # intermediate flag
        # self.overlap_full  = 2    # intermediate flag
        # self.overlap_hole = dqflags.pixel['DO_NOT_USE']
        # self.overlap_no_coverage = dqflags.pixel['NON_SCIENCE'] (also bitwise and with
        # dqflags.pixel['DO_NOT_USE'] )

        # compare the weight plane and spaxel_dq. The initial spaxel_dq flagging
        # has too small a FOV in NIRSpec line mapping case.

        # flatten to match the size of spaxel_weight
        self.spaxel_dq = np.ndarray.flatten(self.spaxel_dq)

        # the fov is an underestimate. Check the spaxel_weight plane
        # if weight map > 0 then set spaxel_dq to overlap_partial
        under_data = self.spaxel_weight > 0
        self.spaxel_dq[under_data] = self.overlap_partial

        # convert all remaining spaxel_dq of 0 to NON_SCIENCE + DO_NOT_USE
        # these pixel should have no overlap with the data
        non_science = self.spaxel_dq == 0
        self.spaxel_dq[non_science] = np.bitwise_or(
            self.overlap_no_coverage, dqflags.pixel["DO_NOT_USE"]
        )

        # refine where good data should be
        ind_full = np.where(np.bitwise_and(self.spaxel_dq, self.overlap_full))
        ind_partial = np.where(np.bitwise_and(self.spaxel_dq, self.overlap_partial))

        self.spaxel_dq[ind_full] = 0
        self.spaxel_dq[ind_partial] = 0

        location_holes = np.where((self.spaxel_dq == 0) & (self.spaxel_weight == 0))
        self.spaxel_dq[location_holes] = self.overlap_hole

        # one last check. Remove pixels flagged as hole but have 1 adjacent spaxel
        # that has no coverage (NON_SCIENCE).  If NON_SCIENCE flag is next to pixel
        # flagged as hole then set the Hole flag to NON_SCIENCE
        spaxel_dq_temp = self.spaxel_dq
        nxy = self.naxis1 * self.naxis2
        index = np.where(self.spaxel_dq == self.overlap_hole)
        for i in range(len(index[0])):
            iwave = int(index[0][i] / nxy)
            rem = index[0][i] - iwave * nxy
            yrem = int(rem / self.naxis1)
            xrem = rem - yrem * self.naxis1

            found = 0
            ij = 0
            # do not allow holes to occur at the edge of IFU cube
            if yrem == 0 or yrem == (self.naxis2 - 1) or xrem == 0 or xrem == (self.naxis1 - 1):
                spaxel_dq_temp[index[0][i]] = np.bitwise_or(
                    self.overlap_no_coverage, dqflags.pixel["DO_NOT_USE"]
                )
                found = 1
            # flag as NON_SCIENCE instead of hole if left, right, top, bottom pixel
            # is NON_SCIENCE
            xcheck = np.zeros(4, dtype=int)
            ycheck = np.zeros(4, dtype=int)
            # left
            xcheck[0] = xrem - 1
            ycheck[0] = yrem
            # right
            xcheck[1] = xrem + 1
            ycheck[1] = yrem
            # bottom
            xcheck[2] = xrem
            ycheck[2] = yrem - 1
            # top
            xcheck[3] = xrem
            ycheck[3] = yrem + 1

            while (ij < 4) and (found == 0):
                if (
                    xcheck[ij] > 0
                    and xcheck[ij] < self.naxis1
                    and ycheck[ij] > 0
                    and ycheck[ij] < self.naxis2
                ):
                    index_check = iwave * nxy + ycheck[ij] * self.naxis1 + xcheck[ij]
                    # If the nearby spaxel_dq contains overlap_no_coverage
                    # then unmark dq flag as hole. A hole has to have nearby
                    # pixels all in FOV.
                    check = (
                        np.bitwise_and(self.spaxel_dq[index_check], self.overlap_no_coverage)
                        == self.overlap_no_coverage
                    )
                    if check:
                        spaxel_dq_temp[index[0][i]] = np.bitwise_or(
                            self.overlap_no_coverage, dqflags.pixel["DO_NOT_USE"]
                        )
                        found = 1
                ij = ij + 1

        self.spaxel_dq = spaxel_dq_temp
        location_holes = np.where(self.spaxel_dq == self.overlap_hole)
        ave_holes = len(location_holes[0]) / self.naxis3

        if ave_holes < 1:
            log.info("Average # of holes/wavelength plane is < 1")
        else:
            log.info("Average # of holes/wavelength plane: %i", ave_holes)
        log.info("Total # of holes for IFU cube is : %i", len(location_holes[0]))

    # ________________________________________________________________________________
    def setup_final_ifucube_model(self, model_ref):
        """
        Set up the final meta WCS info of IFU cube along with other FITS keywords.

        Parameters
        ----------
        model_ref : `~stdatamodels.jwst.datamodels.IFUImageModel`
            The first IFU image model to use to fill in basic header values.

        Returns
        -------
        result : `~stdatamodels.jwst.datamodels.IFUCubeModel`
            IFU cube datamodel with data arrays filled in.
        """
        status = 0
        # loop over the wavelength planes to confirm each plane has some data
        # for initial or final planes that do not have any data - eliminated them
        # from the IFUcube
        # Rearrange values from 1d vectors into 3d cubes

        flux = self.spaxel_flux.reshape((self.naxis3, self.naxis2, self.naxis1))
        wmap = self.spaxel_iflux.reshape((self.naxis3, self.naxis2, self.naxis1))

        var = self.spaxel_var.reshape((self.naxis3, self.naxis2, self.naxis1))
        dq = self.spaxel_dq.reshape((self.naxis3, self.naxis2, self.naxis1))

        # For MIRI MRS, apply a quality cut to help fix spectral tearing at the ends of
        # each band. This is largely taken care of by the WCS regions file, but there
        # will still be 1-2 possibly problematic planes at the end of each band in
        # multi-band cubes. Do this by looking for how many good spaxels there are at
        # each wavelength and finding outliers from the trend.
        if self.instrument == "MIRI":
            nz = flux.shape[0]
            # Create a vector of the number of good spaxels at each wavelength
            ngood = np.zeros(nz)
            for zz in range(0, nz):
                dqvec = dq[zz, :, :].ravel()
                good = np.where(dqvec == 0)
                ngood[zz] = len(good[0])
            # Find where this vector is non-zero, and compute 1% threshold of those good values
            good = np.where(ngood > 0)
            if len(good[0]) > 0:
                pctile = np.percentile(ngood[good], 3)
                # Figure out where the number of good values were less than 75% of threshold,
                # and zero out those arrays.
                lowcov = (np.where((ngood > 0) & (ngood < 0.75 * pctile)))[0]
                nlowcov = len(lowcov)
                log.info("Number of spectral tear planes adjusted: %i", nlowcov)
                for zz in range(0, nlowcov):
                    flux[lowcov[zz], :, :] = 0
                    wmap[lowcov[zz], :, :] = 0
                    var[lowcov[zz], :, :] = 0
                    dq[lowcov[zz], :, :] = (
                        dqflags.pixel["DO_NOT_USE"] + dqflags.pixel["NON_SCIENCE"]
                    )

        # Set np.nan values wherever the DO_NOT_USE flag is set
        dnu = np.where((dq & dqflags.pixel["DO_NOT_USE"]) != 0)
        flux[dnu] = np.nan
        var[dnu] = np.nan

        var = np.sqrt(var)
        if self.linear_wave:
            crval3 = self.crval3
            cdelt3 = self.cdelt3
            crpix3 = self.crpix3
            pixels = np.arange(self.naxis3)

            # Calculate wavelengths
            # We add 1 to 'pixels' to convert 0-based Python indexing to 1-based FITS indexing
            wavelength_table = crval3 + (pixels + 1 - crpix3) * cdelt3
            wave = np.asarray(wavelength_table, dtype=np.float32)
            num = len(wave)
            alldata = np.array([(wave[None].T,)], dtype=[("wavelength", "<f4", (num, 1))])
            # always write the wavetable
            ifucube_model = datamodels.IFUCubeModel(
                data=flux, dq=dq, err=var, weightmap=wmap, wavetable=alldata
            )
        else:
            wave = np.asarray(self.wavelength_table, dtype=np.float32)
            num = len(wave)
            alldata = np.array([(wave[None].T,)], dtype=[("wavelength", "<f4", (num, 1))])

            ifucube_model = datamodels.IFUCubeModel(
                data=flux, dq=dq, err=var, weightmap=wmap, wavetable=alldata
            )

        ifucube_model.update(model_ref)
        if self.output_name is not None:
            ifucube_model.meta.filename = self.output_name

        # Call model_blender if there are multiple inputs
        if len(self.input_models) > 1:
            saved_model_type = ifucube_model.meta.model_type
            self.blend_output_metadata(ifucube_model)
            # Reset to original
            ifucube_model.meta.model_type = saved_model_type
        # ______________________________________________________________________

        ifucube_model.meta.wcsinfo.crval1 = self.crval1
        ifucube_model.meta.wcsinfo.crval2 = self.crval2
        ifucube_model.meta.wcsinfo.crpix1 = self.crpix1
        ifucube_model.meta.wcsinfo.crpix2 = self.crpix2

        ifucube_model.meta.wcsinfo.cdelt1 = self.cdelt1 / 3600.0
        ifucube_model.meta.wcsinfo.cdelt2 = self.cdelt2 / 3600.0
        # Now that we've got a pixel scale, set photometric area keywords
        ifucube_model.meta.photometry.pixelarea_arcsecsq = self.cdelt1 * self.cdelt2
        ifucube_model.meta.photometry.pixelarea_steradians = (
            ifucube_model.meta.photometry.pixelarea_arcsecsq * 2.3504e-11
        )
        if self.linear_wave:
            ifucube_model.meta.wcsinfo.crval3 = self.crval3
            ifucube_model.meta.wcsinfo.cdelt3 = self.cdelt3
            ifucube_model.meta.wcsinfo.ctype3 = "WAVE"
            ifucube_model.meta.wcsinfo.crpix3 = self.crpix3
            ifucube_model.meta.ifu.roi_spatial = float(self.rois)
            ifucube_model.meta.ifu.roi_wave = float(self.roiw)
            # even though we are writing a WAVE-TAB we
            # do not want to set these parameters:
            #   ctype3="WAVE-TAB",
            #   ps3_0 = 'WCS-TABLE'
            #   ps3_1 = 'wavelength'
            # because some viewers (e.g. DS9) report an incorrect wavelength range
        else:
            ifucube_model.meta.wcsinfo.ctype3 = "WAVE-TAB"
            ifucube_model.meta.wcsinfo.ps3_0 = "WCS-TABLE"
            ifucube_model.meta.wcsinfo.ps3_1 = "wavelength"
            ifucube_model.meta.wcsinfo.crval3 = 1.0
            ifucube_model.meta.wcsinfo.crpix3 = 1.0
            ifucube_model.meta.wcsinfo.cdelt3 = None
            ifucube_model.meta.ifu.roi_wave = np.mean(self.roiw_table)
            ifucube_model.wavedim = f"(1,{num:d})"

        ifucube_model.meta.wcsinfo.ctype1 = "RA---TAN"
        ifucube_model.meta.wcsinfo.ctype2 = "DEC--TAN"
        ifucube_model.meta.wcsinfo.cunit1 = "deg"
        ifucube_model.meta.wcsinfo.cunit2 = "deg"

        ifucube_model.meta.wcsinfo.cunit3 = "um"
        ifucube_model.meta.wcsinfo.wcsaxes = 3
        ifucube_model.meta.wcsinfo.pc1_1 = -1
        ifucube_model.meta.wcsinfo.pc1_2 = 0
        ifucube_model.meta.wcsinfo.pc1_3 = 0

        ifucube_model.meta.wcsinfo.pc2_1 = 0
        ifucube_model.meta.wcsinfo.pc2_2 = 1
        ifucube_model.meta.wcsinfo.pc2_3 = 0

        ifucube_model.meta.wcsinfo.pc3_1 = 0
        ifucube_model.meta.wcsinfo.pc3_2 = 0
        ifucube_model.meta.wcsinfo.pc3_3 = 1

        if self.rot_angle is None:
            self.rot_angle = 0.0
        ifucube_model.meta.wcsinfo.pc1_1 = -np.cos(self.rot_angle * np.pi / 180.0)
        ifucube_model.meta.wcsinfo.pc1_2 = np.sin(self.rot_angle * np.pi / 180.0)
        ifucube_model.meta.wcsinfo.pc2_1 = np.sin(self.rot_angle * np.pi / 180.0)
        ifucube_model.meta.wcsinfo.pc2_2 = np.cos(self.rot_angle * np.pi / 180.0)

        ifucube_model.meta.ifu.flux_extension = "SCI"
        ifucube_model.meta.ifu.error_extension = "ERR"
        ifucube_model.meta.ifu.error_type = "ERR"
        ifucube_model.meta.ifu.dq_extension = "DQ"
        ifucube_model.meta.ifu.weighting = str(self.weighting)
        ifucube_model.meta.ifu.weight_power = self.weight_power

        if self.interpolation == "drizzle":
            # stick in values of 0, otherwise it is NaN and
            # fits file can not be written because these
            # values are defined in ifucube.schema.yaml
            ifucube_model.meta.ifu.weight_power = 0
            ifucube_model.meta.ifu.roi_wave = 0
            ifucube_model.meta.ifu.roi_spatial = 0
            ifucube_model.meta.ifu.weighting = str(self.interpolation)


        # set WCS information
        wcsobj = pointing.create_fitswcs(ifucube_model)
        ifucube_model.meta.wcs = wcsobj
        ifucube_model.meta.wcs.bounding_box = (
            (0, self.naxis1 - 1),
            (0, self.naxis2 - 1),
            (0, self.naxis3 - 1),
        )

        ifucube_model.meta.cal_step.cube_build = "COMPLETE"
        # problem with cube_build - contains only 0 data
        if status == 1:
            ifucube_model.meta.cal_step.cube_build = "SKIPPED"

        result = (ifucube_model, status)
        return result

    # ________________________________________________________________________________
    def blend_output_metadata(self, ifu_cube):
        """
        Create new output metadata based on blending all input metadata.

        Parameters
        ----------
        ifu_cube : `~stdatamodels.jwst.datamodels.IFUCubeModel`
            IFU cube data model
        """
        blendmeta.blendmodels(
            ifu_cube,
            self.input_models_this_cube,
            ignore=[
                "meta.filename",
            ],
        )
        # For moving targets, set RA, Dec equal to the average
        mt_avra = getattr(self.input_models_this_cube[0].meta.wcsinfo, "mt_avra", None)
        mt_avdec = getattr(self.input_models_this_cube[0].meta.wcsinfo, "mt_avdec", None)
        if mt_avra is not None:
            ifu_cube.meta.wcsinfo.mt_ra = mt_avra
            ifu_cube.meta.wcsinfo.mt_dec = mt_avdec
            ifu_cube.meta.target.ra = mt_avra
            ifu_cube.meta.target.dec = mt_avdec

    # ________________________________________________________________________________
    def find_ra_dec_offset(self, filename):
        """
        For the given filename find the RA and Dec offset to apply.

        Parameters
        ----------
        filename : str
           Filename that holds the RA and Dec offset to apply

        Returns
        -------
        raoffset : float
            Right ascension offset value read in from file
        decoffset ; float
            Declination offset valuet read in from file
        """
        index = self.offsets["filename"].index(filename)
        raoffset = self.offsets["raoffset"][index]
        decoffset = self.offsets["decoffset"][index]
        return raoffset, decoffset

    # ________________________________________________________________________________
    def offset_coord(self, ra, dec, raoffset, decoffset):
        """
        Given a RA, Dec, RA offset, and Dec offset, use `~astropy.coordinates.SkyCoord` to apply the offsets.

        Parameters
        ----------
        ra : float
            Right ascension coordinate to offset
        dec : float
            Declination coordinate to offset
        raoffset : float
            Right ascension offset to apply
        decoffset : float
            Declination offset to apply

        Returns
        -------
        raw_new : float
            Right ascension coordinate with offset applied
        dec_new : float
            Declination coordinate with offset applied
        """  # noqa: E501
        coord = SkyCoord(ra, dec, unit="deg")
        coord_new = coord.spherical_offsets_by(raoffset, decoffset)

        ra_new = coord_new.ra.value
        dec_new = coord_new.dec.value

        return ra_new, dec_new


class IncorrectInputError(Exception):
    """Interpolation=area when more than 1 file is used to build cube."""

    pass


class IncorrectParameterError(Exception):
    """Cube building parameter is NaN."""

    pass
