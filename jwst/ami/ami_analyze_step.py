from stdatamodels.jwst import datamodels

from jwst.stpipe import Step
from . import ami_analyze
from . import utils

import numpy as np
import asdf

__all__ = ["AmiAnalyzeStep"]

# Affine parameters from commissioning
MX_COMM = 9.92820e-01
MY_COMM = 9.98540e-01
SX_COMM = 6.18605e-03
SY_COMM = -7.27008e-03
XO_COMM = 0.0
YO_COMM = 0.0


class BandpassError(Exception):
    pass


class AmiAnalyzeStep(Step):
    """Performs analysis of an AMI mode exposure by applying the LG algorithm."""

    class_alias = "ami_analyze"

    spec = """
        oversample = integer(default=3, min=1) # Oversampling factor
        rotation = float(default=0.0) # Rotation initial guess [deg]
        psf_offset = string(default='0.0 0.0') # PSF offset values to use to create the model array
        rotation_search = string(default='-3 3 1') # Rotation search parameters: start, stop, step
        bandpass = string(default=None) # ASDF file containing array to override filter/source
        usebp = boolean(default=True) # If True, exclude pixels marked DO_NOT_USE from fringe fitting
        firstfew = integer(default=None) # If not None, process only the first few integrations
        chooseholes = string(default=None) # If not None, fit only certain fringes e.g. ['B4','B5','B6','C2']
        affine2d = string(default='commissioning') # ASDF file containing user-defined affine parameters OR 'commssioning'
        run_bpfix = boolean(default=True) # Run Fourier bad pixel fix on cropped data
    """  # noqa: E501

    reference_file_types = ["throughput", "nrm"]

    def save_model(self, model, *args, **kwargs):
        """
        Override save_model to change suffix based on list of results.

        Parameters
        ----------
        model : data model
            The model to save
        *args, **kwargs : tuple, dict
            Arguments to pass to the stpipe Step.save_model method

        Returns
        -------
        output_paths : [str[, ...]]
            List of output file paths the model(s) were saved in.
        """
        if "idx" in kwargs and kwargs.get("suffix", None) is None:
            kwargs["suffix"] = ["ami-oi", "amimulti-oi", "amilg"][kwargs.pop("idx")]
        return Step.save_model(self, model, *args, **kwargs)

    def override_bandpass(self):
        """
        Read bandpass from asdf file and use it to override the default.

        Expects an array of [effstims, wave_m]
        (i.e. np.array((effstims,wave_m)).T) stored as 'bandpass' in asdf file,
        where effstims are normalized countrates (unitless) and wave_m are the
        wavelengths across the filter at which to compute the model (meters).

        Returns
        -------
        bandpass : array
            Array of [countrates, wavelengths]
        """
        try:
            with asdf.open(self.bandpass, lazy_load=False) as af:
                bandpass = np.array(af["bandpass"])

            # assume it is an array of the correct shape
            wavemin = np.min(bandpass[:, 1])
            wavemax = np.max(bandpass[:, 1])
            self.log.info("User-defined bandpass provided:")
            self.log.info("\tOVERWRITING ALL NIRISS-SPECIFIC FILTER/BANDPASS VARIABLES")
            self.log.info(f"Using {bandpass.shape[0]} wavelengths for fit.")
            self.log.info(f"Wavelength min: {wavemin:.3e} \t Wavelength max: {wavemax:.3e}")

            # update attribute and return
            self.bandpass = bandpass

        except FileNotFoundError as e:
            message = f"File {self.bandpass} could not be found at the specified location."
            raise BandpassError(message) from e

        except KeyError as e:
            message1 = 'ASDF file does not contain the required "bandpass" key. '
            message2 = "See step documentation for info on creating a custom bandpass ASDF file."
            raise BandpassError(message1 + message2) from e

        except (IndexError, ValueError) as e:
            message1 = f"Could not use bandpass from {self.bandpass}. It may have the wrong shape. "
            message2 = "See documentation for info on creating a custom bandpass ASDF file."
            raise BandpassError(message1 + message2) from e

        else:
            return bandpass

    def override_affine2d(self):
        """
        Read user-input affine transform from ASDF file.

        Makes an Affine2d object (see utils.Affine2D class).
        Input should contain mx,my,sx,sy,xo,yo,rotradccw.

        Returns
        -------
        affine2d : Affine2d object
            User-defined affine transform
        """
        msg_defaulting = "\t **** DEFAULTING TO USE IDENTITY TRANSFORM ****"
        try:
            with asdf.open(self.affine2d, lazy_load=False) as af:
                affine2d = utils.Affine2d(
                    mx=af["mx"],
                    my=af["my"],
                    sx=af["sx"],
                    sy=af["sy"],
                    xo=af["xo"],
                    yo=af["yo"],
                    rotradccw=af["rotradccw"],
                )
                self.log.info(f"Using affine transform from ASDF file {self.affine2d}")
            # now self.affine2d updated from string to object
            self.affine2d = affine2d

        except FileNotFoundError:
            self.log.info(f"File {self.affine2d} could not be found at the specified location.")
            self.log.info(msg_defaulting)
            affine2d = None

        except KeyError:
            message1 = (
                "ASDF file does not contain all of the required keys: "
                "mx, my, sx, sy ,xo, yo, rotradccw. "
            )
            message2 = "See step documentation for info on creating a custom affine2d ASDF file."
            self.log.info(message1 + message2)
            self.log.info(msg_defaulting)
            affine2d = None

        except (IndexError, TypeError, ValueError):
            message1 = f"Could not use affine2d from {self.affine2d}. "
            message2 = "See documentation for info on creating a custom bandpass ASDF file."
            self.log.info(message1 + message2)
            self.log.info(msg_defaulting)
            affine2d = None

        else:
            return affine2d

        self.affine2d = affine2d
        return affine2d

    def process(self, input_data):
        """
        Perform analysis of an AMI mode exposure by applying the LG algorithm.

        Parameters
        ----------
        input_data : str or datamodel
            Input file name or datamodel

        Returns
        -------
        oifitsmodel : AmiOIModel object
            AMI tables of median observables from LG algorithm fringe fitting in OIFITS format
        oifitsmodel_multi : AmiOIModel object
            AMI tables of observables for each integration
            from LG algorithm fringe fitting in OIFITS format
        amilgmodel : AmiLGFitModel object
            AMI cropped data, model, and residual data from LG algorithm fringe fitting
        """
        # Retrieve the parameter values
        oversample = self.oversample
        rotate = self.rotation
        bandpass = self.bandpass
        usebp = self.usebp
        firstfew = self.firstfew
        chooseholes = self.chooseholes
        affine2d = self.affine2d
        run_bpfix = self.run_bpfix

        # pull out parameters that are strings and change to floats
        psf_offset = [float(a) for a in self.psf_offset.split()]
        rotsearch_parameters = [float(a) for a in self.rotation_search.split()]

        # handle command-line None input interpreted as string
        if str(affine2d).lower() == "none":
            affine2d = None

        self.log.info(f"Oversampling factor = {oversample}")
        self.log.info(f"Initial rotation guess = {rotate} deg")
        self.log.info(f"Initial values to use for psf offset = {psf_offset}")

        # Make sure oversample is odd
        if oversample % 2 == 0:
            raise ValueError("Oversample value must be an odd integer.")

        # Open the input data model. Can be 2D or 3D image
        with datamodels.open(input_data) as input_model:
            # Get the name of the filter throughput reference file to use
            throughput_reffile = self.get_reference_file(input_model, "throughput")
            self.log.info(f"Using filter throughput reference file {throughput_reffile}")

            # Check for a valid reference file or user-provided bandpass
            if (throughput_reffile == "N/A") & (bandpass is None):
                raise RuntimeError(
                    "No THROUGHPUT reference file found. ami_analyze cannot continue."
                )

            # If there's a user-defined bandpass or affine, handle it
            if bandpass is not None:
                bandpass = self.override_bandpass()
            if affine2d is not None:
                if affine2d == "commissioning":
                    affine2d = utils.Affine2d(
                        mx=MX_COMM,
                        my=MY_COMM,
                        sx=SX_COMM,
                        sy=SY_COMM,
                        xo=XO_COMM,
                        yo=YO_COMM,
                        name="commissioning",
                    )
                    self.log.info("Using affine parameters from commissioning.")
                else:
                    affine2d = self.override_affine2d()

            # Get the name of the NRM reference file to use
            nrm_reffile = self.get_reference_file(input_model, "nrm")
            self.log.info(f"Using NRM reference file {nrm_reffile}")

            with (
                datamodels.ThroughputModel(throughput_reffile) as throughput_model,
                datamodels.NRMModel(nrm_reffile) as nrm_model,
            ):
                # Apply the LG+ methods to the data
                oifitsmodel, oifitsmodel_multi, amilgmodel = ami_analyze.apply_lg_plus(
                    input_model,
                    throughput_model,
                    nrm_model,
                    oversample,
                    psf_offset,
                    rotsearch_parameters,
                    bandpass,
                    usebp,
                    firstfew,
                    chooseholes,
                    affine2d,
                    run_bpfix,
                )

        amilgmodel.meta.cal_step.ami_analyze = "COMPLETE"
        oifitsmodel.meta.cal_step.ami_analyze = "COMPLETE"
        oifitsmodel_multi.meta.cal_step.ami_analyze = "COMPLETE"

        return oifitsmodel, oifitsmodel_multi, amilgmodel
