"""
Module for building 3D spectral cubes from Micro-Shutter Array (MSA) data.

This module defines the `SlitCubeBuildStep`, which processes individual 2D slit spectra
from JWST MSA observations, maps them into a 3D spaxel grid, and constructs an output
spectral cube[cite: 1].
"""

import logging
import time

from stdatamodels.jwst import datamodels

from jwst.assign_wcs.util import update_s_region_keyword
from jwst.cube_build import msa_cube
from jwst.lib.pipe_utils import match_nans_and_flags
from jwst.stpipe import Step, record_step_status

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
         slit_frac = float(default=1.0) # Slit fraction
         ra_center = float(default=None) # RA center of the MSA cube
         dec_center = float(default=None) # Declination center of the MSA cube
         cube_pa = float(default=None) # The position angle of the desired cube in decimal degrees E from N
         nspax_x = integer(default=None) # The odd integer number of spaxels to use in the x dimension of cube tangent plane.
         nspax_y = integer(default=None) # The odd integer number of spaxels to use in the y dimension of cube tangent plane.
         wavemin = float(default=None)  # Minimum wavelength to be used to build the Cube
         wavemax = float(default=None)  # Maximum wavelength to be used to build the Cube
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
        slit_cube :  str, `~stdatamodels.jwst.datamodels.IFUCubeModel`
        """
        log.info("Starting MSA Cube Building Step")

        t0 = time.time()

        # For all parameters convert to a standard format
        # Report read in values to screen

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
        # 2. msaalign (msa cube aligned with slicer plane/ MRS local coord system)

        if self.coord_system == "world":  # world and skyalign are the same things
            self.coord_system = "skyalign"

        log.info(f"Coordinate system to use: {self.coord_system}")

        if self.slit_frac != 1.0:
            if self.slit_frac <= 0 or self.slit_frac > 1.0:
                log.info(f"Slit_frac must be between 0 to 1. It was set to  {self.slit_frac}")
                log.warning("Redefining slit_frac to 1.0")
                self.slit_frac = 1.0

        # Read in the input data and make a copy as needed.
        input_model = self.prepare_output(input_data)

        # 1. Check data model type
        if not isinstance(input_model, datamodels.MultiExposureModel):
            log.warning(
                f"Input dataset is of type '{type(input_model).__name__}'. "
                f"SlitCubeBuildStep requires a MultiExposureModel. Skipping step."
            )
            if input_model is not input_data:
                input_model.close()
            return input_data
        # 2. Check instrument
        instrument = input_model.meta.instrument.name.upper()
        if instrument != "NIRSPEC":
            log.warning(
                f"Input instrument '{instrument}' is not supported. "
                f"SlitCubeBuildStep is designed specifically for NIRSpec. Skipping step."
            )
            if input_model is not input_data:
                input_model.close()
            return input_data

        # grab the grating and filter of the first slit
        # all the slits for a source have the same grating and filter
        grating_name = input_model.exposures[0].meta.instrument.grating
        filter_name = input_model.exposures[0].meta.instrument.filter

        # Using the  reference file designed for IFU cube data
        par_filename = self.get_reference_file(input_model, "cubepar")

        # Check for a valid reference file
        if par_filename == "N/A":
            log.error("No default cube parameters reference file found")
            raise ValueError("The cubepar reference file is required.")

        self.spaxelsize = None
        self.spectralstep = None
        self.wavemin_ref = None
        self.wavemax_ref = None
        self.wavelength_table = None
        log.info("Reading cube parameter file %s", par_filename)
        self.read_cubepars(par_filename, grating_name, filter_name)

        # Override defaults if the user has setup cube parameters
        if self.scalew == 0.0:
            self.scalew = self.spectralstep
        if self.scalexy == 0.0:
            self.scalexy = self.spaxelsize

        # Read in the data from a source cal file and organize slit data
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
            "wavelength_table": self.wavelength_table,
            "suffix": self.suffix,
            "debug_spaxel": self.debug_spaxel,
        }

        # Make sure all input models have consistent NaN and DO_NOT_USE values
        match_nans_and_flags(input_model)

        msacube = msa_cube.MSACubeData(**pars)

        result = msacube.find_footprint()
        corner_ra, corner_dec, final_lam_min, final_lam_max, rot_angle = result

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
        # If the user has not set a wavemin value,  then check that the one determined from the datsa
        # is not larger than the value given in the reference file.
        else:
            if final_lam_max > self.wavemax_ref:
                final_lam_max = (
                    self.wavemax_ref
                )  # TODO Check with NIRSPEC if we should have this check

        msacube.set_slit_wcs(corner_ra, corner_dec, final_lam_min, final_lam_max, rot_angle)

        msacube.print_geometry()

        slit_cube = msacube.build_msacube()

        # irrelevant WCS keywords we will remove from final product
        rm_keys = ["v2_ref", "v3_ref", "ra_ref", "dec_ref", "roll_ref", "v3yangle", "vparity"]

        footprint = slit_cube.meta.wcs.footprint(axis_type="spatial")
        update_s_region_keyword(slit_cube, footprint)

        # remove certain WCS keywords that are irrelevant after combine data into IFUCubes
        for key in rm_keys:
            if key in slit_cube.meta.wcsinfo.instance:
                del slit_cube.meta.wcsinfo.instance[key]

        record_step_status(slit_cube, "slit_cube_build", success=True)

        t1 = time.time()
        log.info(f"Time to build all cubes {t1 - t0}")

        # Output is a new model, so close the input if it was opened here
        if input_model is not input_data:
            input_model.close()

        return slit_cube

    def organize_source(self, model):
        """
        Re-organize all the slits from an exposure into a single source dictionary.

        Parameters
        ----------
        model : `~stdatamodels.jwst.datamodels.MultiExposureModel`
            The input data model containing exposures and slit metadata.

        Returns
        -------
        dict
            Dictionary containing aggregated metadata and list-based attributes
            for each slit exposure associated with the target source.
        """
        exposures = list(model.exposures)
        log.info(f"This source has {len(exposures)} slit observations")

        # Safely extract units from the first slit if exposures exist
        first_slit_meta = exposures[0].meta if exposures else None

        source = {
            "slit_models": exposures,
            "number_slits": len(exposures),
            "bunit": getattr(first_slit_meta, "bunit_data", None),
            "bunit_err": getattr(first_slit_meta, "bunit_err", None),
        }

        self.source = source
        return source

    
    def read_cubepars(self, par_filename, this_grating, this_filter):
        """
        Read default cube parameters from the cubepar reference file based on the grating & filter.

        Parameters
        ----------
        par_filename : str
            Path to the cubepar reference file.
        this_grating : str
            The optical grating name associated with the MSA data (e.g., 'G235H', 'PRISM').
        this_filter : str
            The optical filter name associated with the MSA data (e.g., 'F170LP').

        Raises
        ------
        ValueError
            If no matching row for the grating and filter combination is found in 
            the reference file.
        """
        from astropy.table import Table

        # Read the CUBEPAR extension directly as a table
        cubepar_table = Table.read(par_filename, hdu="CUBEPAR")

        mask = (cubepar_table["disperser"] == this_grating) & (
            cubepar_table["filter"] == this_filter
        )
        target_table = cubepar_table[mask]
        if len(target_table) == 0:
            raise ValueError(
                f"No matching configuration found for {this_grating} and {this_filter}."
            )

        # Extract the single Row object (or index [0]) to get scalar values
        target_row = target_table[0]
        self.spaxelsize = target_row["spaxelsize"]
        self.spectralstep = target_row["spectralstep"]
        self.wavemin_ref = target_row["wavemin"]
        self.wavemax_ref = target_row["wavemax"]

        if not self.linear_wave:
            if "PRISM" in this_grating:
                wavelength_table = Table.read(par_filename, "MULTICHAN_PRISM_DRIZZLE")
            elif "M" in this_grating:
                wavelength_table = Table.read(par_filename, "MULTICHAN_MED_DRIZZLE")
            elif "H" in this_grating:
                wavelength_table = Table.read(par_filename, "MULTICHAN_HIGH_DRIZZLE")
            else:
                raise ValueError(
                    f"Invalid grating '{this_grating}': Unable to read cube parameters."
                )
            self.wavelength_table = wavelength_table["wavelength"]

        log.info("MSA cube parameters read in from the reference file:")
        log.info(f"Spaxelsize {self.spaxelsize:.3f} ")
        if self.linear_wave:
            log.info(f"Spectralstep {self.spectralstep: .6f} ")
            log.info(
                f" Minimum wavelength {self.wavemin_ref: .2f}"
                f" Maximum wavelength {self.wavemax_ref: .2f}"
            )
        else:
            log.info(
                f" Number of wavelength elements in wavelength table {self.wavelength_table.shape} "
            )
