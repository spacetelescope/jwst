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
    "NRS_AUTOWAVE",
    "NRS_AUTOFLAT",
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
    "NRS_LAMP",
    "NRS_AUTOWAVE",
]


__all__ = ["FlatFieldStep"]


class FlatFieldStep(Step):
    """Flat-field a science image using a flatfield reference image.

    Attributes
    ----------
    correction_pars : {'flat': DataModel}
        After the step has successfully run, the flat field applied is
        stored.

    use_correction_pars : boolean
        Use the flat stored in `correction_pars`
    """

    class_alias = "flat_field"

    spec = """
        save_interpolated_flat = boolean(default=False) # Save interpolated NRS flat
        user_supplied_flat = string(default=None)  # User-supplied flat
        inverse = boolean(default=False)  # Invert the operation
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
                exposure_type not in EXTRACT_2D_IS_OK):
            self.log.warning("The extract_2d step has been run, but for "
                             "%s data it should not have been run, so ...",
                             exposure_type)
            self.log.warning("flat fielding will be skipped.")
            return self.skip_step(input_model)

        # Retrieve reference files only if no user-supplied flat is specified
        if self.user_supplied_flat is not None:
            self.log.info(
                f'User-supplied flat {self.user_supplied_flat} given.'
                ' Ignoring all flat reference files and flat creation.'
            )
            reference_file_models = {
                'user_supplied_flat': datamodels.open(self.user_supplied_flat)
            }

            # Record the user-supplied flat as the FLAT reference type for recording
            # in the result header.
            self._reference_files_used.append(
                ('flat', reference_file_models['user_supplied_flat'].meta.filename)
            )
        elif self.use_correction_pars:
            self.log.info(f'Using flat field from correction pars {self.correction_pars["flat"]}')
            reference_file_models = {
                'user_supplied_flat': datamodels.open(self.correction_pars['flat'])
            }

            # Record the flat as the FLAT reference type for recording
            # in the result header.
            self._reference_files_used.append(
                ('flat', reference_file_models['user_supplied_flat'].meta.filename)
            )
        else:
            reference_file_models = self._get_references(input_model, exposure_type)

        # Do the flat-field correction
        output_model, flat_applied = flat_field.do_correction(
            input_model,
            **reference_file_models,
            inverse=self.inverse
        )

        # Close the input and reference files
        input_model.close()
        try:
            for model in reference_file_models.values():
                model.close()
        except AttributeError:
            pass

        if self.save_interpolated_flat and flat_applied is not None:
            ff_path = self.save_model(flat_applied, suffix=self.flat_suffix, force=True)
            self.log.info(f'Interpolated flat written to "{ff_path}".')

        if not self.correction_pars:
            self.correction_pars = {}
        self.correction_pars['flat'] = flat_applied

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

    def _get_references(self, data, exposure_type):
        """Retrieve required CRDS reference files

        Parameters
        ----------
        data : DataModel
            The data to base the CRDS lookups on.

        exposure_type : str
            The exposure type keyword value

        Returns
        -------
        reference_file_models : {str: DataModel{,...}}
            Dictionary matching reference file types to open models
        """

        # Get reference file paths
        reference_file_names = {}
        for reftype in self.reference_file_types:
            reffile = self.get_reference_file(data, reftype)
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
                self.log.debug('No reference found for type %s', reftype.upper())
                reference_file_models[reftype] = None

        return reference_file_models
