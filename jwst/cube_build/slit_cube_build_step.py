import copy
import logging
import time

from jwst.cube_build import msa_cube
from jwst.lib.pipe_utils import match_nans_and_flags
from jwst.stpipe import Step

__all__ = ["SlitCubeBuildStep"]

log = logging.getLogger(__name__)


class SlitCubeBuildStep(Step):
    """Create a 3-D spectral cube from MSA data."""

    class_alias = "slit_cube_build"

    spec = """
         linear_wave = boolean(default=True) # Toggle between linear (True) and nonlinear (False) wavelength dimensions
         scalexy = float(default=0.0) # cube sample size to use for axis 1 and axis2, arc seconds
         scalew = float(default=0.0) # cube sample size to use for axis 3, microns
         weighting = option('volume','readnoise',default = 'volume') # Type of weighting function
         coord_system = option('skyalign','world','msaalign',default='skyalign') # Output Coordinate system.
         ra_center = float(default=None) # RA center of the IFU cube
         dec_center = float(default=None) # Declination center of the IFU cube
         cube_pa = float(default=None) # The position angle of the desired cube in decimal degrees E from N
         nspax_x = integer(default=None) # The odd integer number of spaxels to use in the x dimension of cube tangent plane.
         nspax_y = integer(default=None) # The odd integer number of spaxels to use in the y dimension of cube tangent plane.
         wavemin = float(default=None)  # Minimum wavelength to be used in the IFUCube
         wavemax = float(default=None)  # Maximum wavelength to be used in the IFUCube
         output_use_model = boolean(default=true) # Use filenames in the output models
         suffix = string(default='s3d')
         debug_spaxel = string(default='-1 -1 -1') # Default not used
       """  # noqa: E501

    reference_file_types = ["cubepar"]

    def process(self, input_data):
        """
        Build an spectral  cube from overlapping MSA slitlet stepped data.

        This is the controlling routine for building MSA spectral cubes.
        It loads and sets the various input data and parameters needed by
        the ``slit_cube_build`` step.

        This routine does the following operations:

        1. Extracts the input parameters from the cubepars reference file and
           merges them with any user-provided values.
        2. Creates the output WCS from the input images and defines the mapping
           between all the input arrays and the output array.
        3. Passes the input data to the function to map all their input data
           to the output array.
        4. Updates the output data model with correct meta data.

        Parameters
        ----------
        input_data : str, `~stdatamodels.jwst.datamodels.MultExposureModel`
           Datamodel file name, MultiExposureModel
           The input is expected to several 2D spectral images, to combine together into a
           spectral cube.

        Returns
        -------
        cube_container : `~jwst.datamodels.container.ModelContainer`
           Container (list) of `~stdatamodels.jwst.datamodels.MSACubeModel`.
        """
        log.info("Starting MSA Cube Building Step")

        t0 = time.time()
        # ________________________________________________________________________________
        # For all parameters convert to a standard format
        # Report read in values to screen
        # ________________________________________________________________________________

        if not self.coord_system.islower():
            self.coord_system = self.coord_system.lower()

        if not self.weighting.islower():
            self.weighting = self.weighting.lower()

        if self.scalexy != 0.0:
            log.info(f"Input Scale of axis 1 and 2 {self.scalexy}")
        if self.scalew != 0.0:
            log.info(f"Input wavelength scale {self.scalew}")

        if self.wavemin is not None:
            log.info(f"Setting minimum wavelength of spectral cube to: {self.wavemin}")
        if self.wavemax is not None:
            log.info(f"Setting maximum wavelength of spectral cube to: {self.wavemax}")

        # check that if self.nspax_x or self.nspax_y is provided they must be odd numbers
        if self.nspax_x is not None:
            if self.nspax_x % 2 == 0:
                log.info(f"Input nspax_x must be an odd number {self.nspax_x}")
                self.nspax_x = self.nspax_x + 1
                log.info(f"Updating nspa by 1. New value {self.nspax_x}")

        if self.nspax_y is not None:
            if self.nspax_y % 2 == 0:
                log.info(f"Input nspax_y must be an odd number {self.nspax_y}")
                self.nspax_y = self.nspax_y + 1
                log.info(f"Updating nspax_y by 1. New value {self.nspax_y}")

        # valid coord_system:
        # 1. skyalign (ra dec) (aka world)
        # 2. msaalign (ifu cube aligned with slicer plane/ MRS local coord system)

        if self.coord_system == "world":  # world and skyalign are the same things
            self.coord_system = "skyalign"

        log.info(f"Coordinate system to use: {self.coord_system}")


        self.slit_frac = 1.0 # change this later to be an input parameter
        
        # Read in the input data and make a copy as needed.
        input_model = self.prepare_output(input_data)

        print("Type data models", input_model)
        # grab the grating and filter of the first slit (all the slits have the same grating and filter)
        grating = input_model.exposures[0].meta.instrument.grating
        filter = input_model.exposures[0].meta.instrument.filter

        # Test using the IFU cube pars reference file
        par_filename = self.get_reference_file(input_model, "cubepar")

        # Check for a valid reference file
        if par_filename == "N/A":
            log.error("No default cube parameters reference file found")
            raise ValueError("The cubepar reference file is required.")

        self.spaxelsize = None
        self.spectralstep = None
        self.wavemin_ref = None
        self.wavemax_ref = None
        self.wavetable = None
        log.info("Reading cube parameter file %s", par_filename)
        self.read_cubepars(par_filename, grating, filter)

        # Override defaults if the user has setup cube parameters
        if self.scalew == 0.0:
            self.scalew = self.spectralstep
        if self.scalexy == 0.0:
            self.scalexy = self.spaxelsize
        print("After reading ref file", self.scalew, self.scalexy)
        # Read in the data and organize organize source file
        self.organize_source(input_model)

        pars = {
            "source": self.source,
            "scalexy": self.scalexy,
            "scalew": self.scalew,
            "weighting": self.weighting,
            "coord_system": self.coord_system,
            "linear_wave": self.linear_wave,
            "ra_center": self.ra_center,
            "dec_center": self.dec_center,
            "cube_pa": self.cube_pa,
            "nspax_x": self.nspax_x,
            "nspax_y": self.nspax_y,
            "slit_frac": self.slit_frac,
            "wavemin": self.wavemin,
            "wavemax": self.wavemax,
            "wavetable": self.wavetable,
            "suffix": self.suffix,
            "debug_spaxel": self.debug_spaxel,
        }

        # Make sure all input models have consistent NaN and DO_NOT_USE values
        match_nans_and_flags(input_model)

        status_cube = 0
        msacube = msa_cube.MSACubeData(**pars)

        result = msacube.find_footprint()
        corner_ra, corner_dec, final_lam_min, final_lam_max, rot_angle = result

        # Check if the user set wavemin. If so use that value.
        if self.wavemin is not None:
            final_lam_min = self.wavemin
        # Wave min is min of (final_lam_min, self.wavemin_ref): from data or from reference file
        else:
            if self.wavemin_ref < final_lam_min:
                final_lam_min = self.wavemin_ref

        # user set wavemax use that values
        if self.wavemax is not None:
            final_lam_max = self.wavemax
        # Set wave max. The maximum values can not be larger than the value provide in the reference file
        # If the user wants it set larger then should set the parameter --wavemax
        else:
            if final_lam_max > self.wavemax_ref:
                final_lam_min = self.wavemax_ref

        msacube.set_slit_wcs(corner_ra, corner_dec, final_lam_min, final_lam_max, rot_angle)

        #print("Mapped coordinate system")
        #print("naxis, cdelt, crpix, crval")
        #print("Axis 1 ", self.naxis1, self.cdelt1, self.crpix1, self.crval1)
        #print("Axis 2", self.naxis2, self.cdelt2, self.crpix2, self.crval2)
        #print("Axis 3", self.naxis3, self.cdelt3, self.crpix3, self.crval3)
        #print("crval1 crval2", self.crval1, self.crval2)
        #print("crpix1 crpix2", self.crpix1, self.crpix2)
        #print("cdelt1 cdelt2", self.cdelt1, self.cdelt2)

        msacube.print_geometry()
        
        status = 0
        result, status = msacube.build_msacube()

        # check if cube_build failed
        # if status == 1:
        #    status_cube = 1

        # irrelevant WCS keywords we will remove from final product
        # rm_keys = ["v2_ref", "v3_ref", "ra_ref", "dec_ref", "roll_ref", "v3yangle", "vparity"]

        # for cube in cube_container:
        #    footprint = cube.meta.wcs.footprint(axis_type="spatial")
        #    update_s_region_keyword(cube, footprint)

        # remove certain WCS keywords that are irrelevant after combine data into IFUCubes
        #    for key in rm_keys:
        #        if key in cube.meta.wcsinfo.instance:
        #            del cube.meta.wcsinfo.instance[key]
        # if status_cube == 1:
        #    record_step_status(cube_container, "cube_build", success=False)
        # else:
        #    record_step_status(cube_container, "cube_build", success=True)

        # t1 = time.time()
        # log.debug(f"Time to build all cubes {t1 - t0}")

        # Output is a new model, so close the input if it was opened here
        # if read_in_models is not input_data:
        #    read_in_models.close()

        # return cube_container

    def organize_source(self, model):

        source = {}
        detector = model.meta.instrument.detector
        source["slitname"] = []
        source["slit_models"] = []
        source["wcs"] = []
        source["source_type"] = []
        source["source_id"] = []
        source["bunit"] = []
        source["bunit_err"] = []
        source["pixelarea_ster"] = []
        source["pixelarea_arcsec"] = []
        source["detector"] = []
        source["slitnum"] = []

        source["number_slits"] = len(model.exposures)
        print("this source has ", len(model.exposures), "slit observations")

        for j, slit in enumerate(model.exposures):
            source_id = slit.source_id
            stype = slit.source_type
            slitname = slit.name
            name = slit.source_name
            bunit = slit.meta.bunit_data
            bunit_err = slit.meta.bunit_err
            parea_steradians = slit.meta.photometry.pixelarea_steradians
            parea_arcsec = slit.meta.photometry.pixelarea_arcsecsq
            wcs = copy.deepcopy(slit.meta.wcs)

            source["slitname"].append(slitname)
            source["source_type"].append(stype)
            source["slit_models"].append(slit)
            source["wcs"].append(wcs)
            source["pixelarea_ster"] = parea_steradians
            source["pixelarea_arcsec"] = parea_arcsec
            source["source_id"].append(source_id)
            source["detector"].append(detector)
            source["bunit"] = bunit
            source["bunit_err"] = bunit_err
            source["slitnum"].append(j)

        self.source = source
        return source

    def read_cubepars(self, par_filename, this_grating, this_filter):
        """
        Read in :ref:`cubepar_reffile`.
        Based on filter and grating read in the appropriate columns in the
        :ref:`cubepar_reffile` and return defaults

        Parameters
        ----------
        par_filename : str
        Cube parameter reference filename
        grating : str
        The grating for the MSA input data
        filter : list
        The filter for the MSA input data

        Returns
        -------

        """
        from astropy.table import Table

        # Read the CUBEPAR extension directly as a table
        cubepar_table = Table.read(par_filename, hdu="CUBEPAR")

        # Filter for grating (disperser) G235H and filter F170LP
        # mask = (cubepar_table["disperser"] == "G235H") & (
        #    cubepar_table["filter"] == "F170LP"

        mask = (cubepar_table["disperser"] == this_grating) & (
            cubepar_table["filter"] == this_filter
        )
        target_table = cubepar_table[mask]
        if len(target_table) == 0:
            raise ValueError("No matching configuration found for G235H and F170LP.")

        # Extract the single Row object (or index [0]) to get scalar values
        target_row = target_table[0]
        print(" Row in cubepar for this configuration", target_row)
        self.spaxelsize = target_row["spaxelsize"]
        self.spectralstep = target_row["spectralstep"]
        self.wavemin_ref = target_row["wavemin"]
        self.wavemax_ref = target_row["wavemax"]

        if "PRISM" in this_grating:
            self.wavetable = Table.read(par_filename, "MULTICHAN_PRISM_DRIZZLE")
        elif "M" in this_grating:
            self.wavetable = Table.read(par_filename, "MULTICHAN_MED_DRIZZLE")
        elif "H" in this_grating:
            self.wavetable = Table.read(par_filename, "MULTICHAN_HIGH_DRIZZLE")
        else:
            print("Invalid grating", this_grating)

        print("spaxelsize", self.spaxelsize)
        print("spectralstep", self.spectralstep)
        print("wave min and max", self.wavemin_ref, self.wavemax_ref)
