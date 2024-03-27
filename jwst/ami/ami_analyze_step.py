from stdatamodels.jwst import datamodels

from ..stpipe import Step
from . import ami_analyze

__all__ = ["AmiAnalyzeStep"]


class AmiAnalyzeStep(Step):
    """Performs analysis of an AMI mode exposure by applying the LG algorithm."""

    class_alias = "ami_analyze"

    spec = """
        oversample = integer(default=3, min=1)  # Oversampling factor
        rotation = float(default=0.0)           # Rotation initial guess [deg]
        psf_offset = string(default='0.0 0.0') # PSF offset values to use to create the model array
        rotation_search = string(default='-3 3 1') # Rotation search parameters: start, stop, step
        bandpass = any(default=None) # Synphot spectrum or array to override filter/source
        usebp = boolean(default=True) # If True, exclude pixels marked DO_NOT_USE from fringe fitting
        firstfew = integer(default=None) # If not None, process only the first few integrations
        chooseholes = string(default=None) # If not None, fit only certain fringes e.g. ['B4','B5','B6','C2']
        affine2d = any(default=None) # None or user-defined Affine2d object
        run_bpfix = boolean(default=True) # Run Fourier bad pixel fix on cropped data
    """

    reference_file_types = ['throughput', 'nrm']

    def save_model(self, model, *args, **kwargs):
        # Override save_model to change suffix based on list of results
        if "idx" in kwargs and kwargs.get("suffix", None) is None:
            kwargs["suffix"] = ["ami-oi", "amimulti-oi", "amilg"][kwargs.pop("idx")]
        return Step.save_model(self, model, *args, **kwargs)

    def process(self, input):
        """
        Performs analysis of an AMI mode exposure by applying the LG algorithm.

        Parameters
        ----------
        input: string
            input file name

        Returns
        -------
        oifitsmodel: AmiOIModel object
            AMI tables of median observables from LG algorithm fringe fitting in OIFITS format
        oifitsmodel_multi: AmiOIModel object
            AMI tables of observables for each integration from LG algorithm fringe fitting in OIFITS format
        amilgmodel: AmiLGFitModel object
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

        self.log.info(f"Oversampling factor = {oversample}")
        self.log.info(f"Initial rotation guess = {rotate} deg")
        self.log.info(f"Initial values to use for psf offset = {psf_offset}")

        # Make sure oversample is odd
        if oversample % 2 == 0:
            raise ValueError("Oversample value must be an odd integer.")

        # Open the input data model. Can be 2D or 3D image
        with datamodels.open(input) as input_model:
            # Get the name of the filter throughput reference file to use
            throughput_reffile = self.get_reference_file(input_model, 'throughput')
            self.log.info(f'Using filter throughput reference file {throughput_reffile}')

            # Check for a valid reference file or user-provided bandpass
            if (throughput_reffile == 'N/A') & (bandpass is None):
                raise RuntimeError("No THROUGHPUT reference file found. "
                                   "ami_analyze cannot continue.")

            # Get the name of the NRM reference file to use
            nrm_reffile = self.get_reference_file(input_model, 'nrm')
            self.log.info(f'Using NRM reference file {nrm_reffile}')

            with (
                datamodels.ThroughputModel(throughput_reffile) as throughput_model,
                datamodels.NRMModel(nrm_reffile) as nrm_model,
            ):
                # Apply the LG+ methods to the data
                oifitsmodel, oifitsmodel_multi, amilgmodel = ami_analyze.apply_LG_plus(
                    input_model,
                    throughput_model,
                    nrm_model,
                    oversample,
                    rotate,
                    psf_offset,
                    rotsearch_parameters,
                    bandpass,
                    usebp,
                    firstfew,
                    chooseholes,
                    affine2d,
                    run_bpfix,
                )

        amilgmodel.meta.cal_step.ami_analyze = 'COMPLETE'
        oifitsmodel.meta.cal_step.ami_analyze = 'COMPLETE'
        oifitsmodel_multi.meta.cal_step.ami_analyze = 'COMPLETE'

        return oifitsmodel, oifitsmodel_multi, amilgmodel
