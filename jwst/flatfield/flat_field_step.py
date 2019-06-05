from ..stpipe import Step
from .. import datamodels
from . import flat_field

# For the following types of data, it is OK -- and in some cases
# required -- for the extract_2d step to have been run.  For all
# other types of data, the extract_2d step must not have been run.
EXTRACT_2D_IS_OK = [
    "NRS_BRIGHTOBJ",
    "NRS_FIXEDSLIT",
    "NRS_LAMP",
    "NRS_MSASPEC",
    ]

# NIRSpec imaging types (see exp_type2transform in assign_wcs/nirspec.py)
NRS_IMAGING_MODES = [
    "NRS_CONFIRM",
    "NRS_FOCUS",
    "NRS_IMAGE",
    "NRS_MIMF",
    "NRS_MSATA",
    "NRS_TACONFIRM",
    "NRS_TACQ",
    "NRS_TASLIT",
    "NRS_WATA",
    ]
# Supported NIRSpec spectrographic types. No flat fielding for NRS_AUTOFLAT
NRS_SPEC_MODES = [
    "NRS_BRIGHTOBJ",
    "NRS_FIXEDSLIT",
    "NRS_IFU",
    "NRS_MSASPEC",
    ]


__all__ = ["FlatFieldStep"]


class FlatFieldStep(Step):
    """Flat-field a science image using a flatfield reference image.
    """

    spec = """
        save_interpolated_flat = boolean(default=False) # Save interpolated NRS flat
    """

    reference_file_types = ["flat", "fflat", "sflat", "dflat"]

    # Define a suffix for optional saved output of the interpolated flat for NRS
    flat_suffix = 'interpolatedflat'

    def process(self, input):

        input_model = datamodels.open(input)
        exposure_type = input_model.meta.exposure.type.upper()

        self.log.debug("Input is {} of exposure type {}".format(
            input_model.__class__.__name__, exposure_type))

        if input_model.meta.instrument.name.upper() == "NIRSPEC":
            if (exposure_type not in NRS_SPEC_MODES and
                exposure_type not in NRS_IMAGING_MODES):
                self.log.warning("Exposure type is %s; flat-fielding will be "
                                 "skipped because it is not currently "
                                 "supported for this mode.", exposure_type)
                return self.skip_step(input_model)

        # Check whether extract_2d has been run.
        if (input_model.meta.cal_step.extract_2d == 'COMPLETE' and
            not exposure_type in EXTRACT_2D_IS_OK):
            self.log.warning("The extract_2d step has been run, but for "
                             "%s data it should not have been run, so ...",
                             exposure_type)
            self.log.warning("flat fielding will be skipped.")
            return self.skip_step(input_model)

        # Get reference file paths
        reference_file_names = {}
        for reftype in self.reference_file_types:
            reffile = self.get_reference_file(input_model, reftype)
            reference_file_names[reftype] = reffile if reffile != 'N/A' else None

        # Define mapping between reftype and datamodel type
        model_type = dict(
            flat=datamodels.FlatModel,
            fflat=datamodels.NirspecFlatModel,
            sflat=datamodels.NirspecFlatModel,
            dflat=datamodels.NirspecFlatModel,
            )
        if exposure_type == "NRS_MSASPEC":
            model_type["fflat"] = datamodels.NirspecQuadFlatModel

        # Open the relevant reference files as datamodels
        reference_file_models = {}
        for reftype, reffile in reference_file_names.items():
            if reffile is not None:
                reference_file_models[reftype] = model_type[reftype](reffile)
                self.log.debug('Using %s reference file: %s', reftype.upper(), reffile)
            else:
                reference_file_models[reftype] = None

        # Do the flat-field correction
        output_model, interpolated_flats = flat_field.do_correction(
            input_model,
            **reference_file_models,
            )

        # Close the input and reference files
        input_model.close()
        try:
            for model in reference_file_models.values():
                model.close()
        except AttributeError:
            pass

        if self.save_interpolated_flat and interpolated_flats is not None:
            self.log.info("Writing interpolated flat field.")
            self.save_model(interpolated_flats, suffix=self.flat_suffix)
            interpolated_flats.close()

        return output_model

    def skip_step(self, input_model):
        """Set the calibration switch to SKIPPED.

        This method makes a copy of input_model, sets the calibration
        switch for the flat_field step to SKIPPED in the copy, closes
        input_model, and returns the copy.
        """

        result = input_model.copy()
        result.meta.cal_step.flat_field = "SKIPPED"
        input_model.close()
        return result
