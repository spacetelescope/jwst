import crds
from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelContainer, SourceModelContainer
from jwst.datamodels.utils.wfss_multispec import make_wfss_multiexposure
from jwst.extract_1d import extract
from jwst.extract_1d.ifu import ifu_extract1d
from jwst.extract_1d.soss_extract import soss_extract
from jwst.stpipe import Step

__all__ = ["Extract1dStep"]


class Extract1dStep(Step):
    """Extract a 1D spectrum from 2D data."""

    class_alias = "extract_1d"

    spec = """
    subtract_background = boolean(default=None)  # subtract background?
    apply_apcorr = boolean(default=True)  # apply aperture corrections?

    extraction_type = option("box", "optimal", None, default="box") # Perform box or optimal extraction
    use_source_posn = boolean(default=None)  # use source coords to center extractions?
    position_offset = float(default=0)  # number of pixels to shift source trace in the cross-dispersion direction
    model_nod_pair = boolean(default=True)  # For optimal extraction, model a negative nod if possible
    optimize_psf_location = boolean(default=True)  # For optimal extraction, optimize source location
    smoothing_length = integer(default=None)  # background smoothing size
    bkg_fit = option("poly", "mean", "median", None, default=None)  # background fitting type
    bkg_order = integer(default=None, min=0)  # order of background polynomial fit
    log_increment = integer(default=50)  # increment for multi-integration log messages
    save_profile = boolean(default=False)  # save spatial profile to disk
    save_scene_model = boolean(default=False)  # save flux model to disk
    save_residual_image = boolean(default=False)  # save residual image to disk

    center_xy = float_list(min=2, max=2, default=None)  # IFU extraction x/y center
    ifu_autocen = boolean(default=False) # Auto source centering for IFU point source data.
    bkg_sigma_clip = float(default=3.0)  # background sigma clipping threshold for IFU
    ifu_rfcorr = boolean(default=True) # Apply 1d residual fringe correction (MIRI MRS only)
    ifu_set_srctype = option("POINT", "EXTENDED", None, default=None) # user-supplied source type
    ifu_rscale = float(default=None, min=0.5, max=3) # Radius in terms of PSF FWHM to scale extraction radii
    ifu_covar_scale = float(default=1.0) # Scaling factor to apply to errors to account for IFU cube covariance

    soss_atoca = boolean(default=True)  # use ATOCA algorithm
    soss_threshold = float(default=1e-2)  # TODO: threshold could be removed from inputs. Its use is too specific now.
    soss_n_os = integer(default=2)  # minimum oversampling factor of the underlying wavelength grid used when modeling trace.
    soss_wave_grid_in = input_file(default = None)  # Input wavelength grid used to model the detector
    soss_wave_grid_out = string(default = None)  # Output wavelength grid solution filename
    soss_estimate = input_file(default = None)  # Estimate used to generate the wavelength grid
    soss_rtol = float(default=1.0e-4)  # Relative tolerance needed on a pixel model
    soss_max_grid_size = integer(default=20000)  # Maximum grid size, if wave_grid not specified
    soss_tikfac = float(default=None)  # regularization factor for NIRISS SOSS extraction
    soss_width = float(default=40.)  # aperture width used to extract the 1D spectrum from the de-contaminated trace.
    soss_bad_pix = option("model", "masking", default="masking")  # method used to handle bad pixels
    soss_modelname = output_file(default = None)  # Filename for optional model output of traces and pixel weights
    """  # noqa: E501

    reference_file_types = ["extract1d", "apcorr", "pastasoss", "specprofile", "speckernel", "psf"]

    def _get_extract_reference_files_by_mode(self, model, exp_type):
        """
        Get extraction reference files with special handling by exposure type.

        Parameters
        ----------
        model : DataModel
            The input model.
        exp_type : str
            Exposure type.

        Returns
        -------
        extract_ref : str
            EXTRACT1D reference file or N/A.
        apcorr_ref : str
            APCORR reference file or N/A.
        psf_ref : str
            PSF reference file or N/A.
        """
        if isinstance(model, ModelContainer):
            model = model[0]

        if exp_type not in extract.WFSS_EXPTYPES:
            extract_ref = self.get_reference_file(model, "extract1d")
        else:
            extract_ref = "N/A"
        if extract_ref != "N/A":
            self.log.info(f"Using EXTRACT1D reference file {extract_ref}")

        if self.apply_apcorr:
            apcorr_ref = self.get_reference_file(model, "apcorr")
        else:
            apcorr_ref = "N/A"
        if apcorr_ref != "N/A":
            self.log.info(f"Using APCORR file {apcorr_ref}")

        try:
            psf_ref = self.get_reference_file(model, "psf")
        except crds.core.exceptions.CrdsLookupError:
            psf_ref = "N/A"
        if psf_ref != "N/A":
            self.log.info(f"Using PSF reference file {psf_ref}")

        return extract_ref, apcorr_ref, psf_ref

    def _extract_soss(self, model):
        """
        Extract NIRISS SOSS spectra.

        Parameters
        ----------
        model : DataModel
            Input model.

        Returns
        -------
        DataModel
            The output spectra.
        """
        # Work on a copy of the model
        model = model.copy()

        # Set the filter configuration
        if model.meta.instrument.filter == "CLEAR":
            self.log.info("Exposure is through the GR700XD + CLEAR (science).")
            soss_filter = "CLEAR"
        else:
            self.log.error(
                "The SOSS extraction is implemented for the CLEAR filter only. "
                f"Requested filter is {model.meta.instrument.filter}."
            )
            self.log.error("extract_1d will be skipped.")
            model.meta.cal_step.extract_1d = "SKIPPED"
            return model

        # Set the subarray mode being processed
        if model.meta.subarray.name == "SUBSTRIP256":
            self.log.info("Exposure is in the SUBSTRIP256 subarray.")
            self.log.info("Traces 1 and 2 will be modelled and decontaminated before extraction.")
            subarray = "SUBSTRIP256"
        elif model.meta.subarray.name == "SUBSTRIP96":
            self.log.info("Exposure is in the SUBSTRIP96 subarray.")
            self.log.info(
                "Traces of orders 1 and 2 will be modelled but only order 1 "
                "will be decontaminated before extraction."
            )
            subarray = "SUBSTRIP96"
        else:
            self.log.error(
                "The SOSS extraction is implemented for the SUBSTRIP256 "
                "and SUBSTRIP96 subarrays only. Subarray is currently "
                f"{model.meta.subarray.name}."
            )
            self.log.error("Extract1dStep will be skipped.")
            model.meta.cal_step.extract_1d = "SKIPPED"
            return model

        # Load reference files.
        pastasoss_ref_name = self.get_reference_file(model, "pastasoss")
        specprofile_ref_name = self.get_reference_file(model, "specprofile")
        speckernel_ref_name = self.get_reference_file(model, "speckernel")

        # Build SOSS kwargs dictionary.
        soss_kwargs = {}
        soss_kwargs["threshold"] = self.soss_threshold
        soss_kwargs["n_os"] = self.soss_n_os
        soss_kwargs["tikfac"] = self.soss_tikfac
        soss_kwargs["width"] = self.soss_width
        soss_kwargs["bad_pix"] = self.soss_bad_pix
        soss_kwargs["subtract_background"] = self.subtract_background
        soss_kwargs["rtol"] = self.soss_rtol
        soss_kwargs["max_grid_size"] = self.soss_max_grid_size
        soss_kwargs["wave_grid_in"] = self.soss_wave_grid_in
        soss_kwargs["wave_grid_out"] = self.soss_wave_grid_out
        soss_kwargs["estimate"] = self.soss_estimate
        soss_kwargs["atoca"] = self.soss_atoca
        # Set flag to output the model and the tikhonov tests
        soss_kwargs["model"] = True if self.soss_modelname else False

        # Run the extraction.
        result, ref_outputs, atoca_outputs = soss_extract.run_extract1d(
            model,
            pastasoss_ref_name,
            specprofile_ref_name,
            speckernel_ref_name,
            subarray,
            soss_filter,
            soss_kwargs,
        )

        # Set the step flag to complete
        if result is None:
            return None
        else:
            result.meta.cal_step.extract_1d = "COMPLETE"
            result.meta.target.source_type = None

            model.close()

            if self.soss_modelname:
                soss_modelname = self.make_output_path(
                    basepath=self.soss_modelname, suffix="SossExtractModel"
                )
                ref_outputs.save(soss_modelname)

        if self.soss_modelname:
            soss_modelname = self.make_output_path(
                basepath=self.soss_modelname, suffix="AtocaSpectra"
            )
            atoca_outputs.save(soss_modelname)

        return result

    def _check_mrs_type(self, model):
        """
        Check if the MIRI MRS data is for a single band.

        Parameters
        ----------
        model : DataModel
            Input model.

        Returns
        -------
        band_check : bool
            Flag is data is for a single MIRI MRS band.
        """
        # check channel is only 1 value and band is only 1 value
        validch = ["1", "2", "3", "4"]
        validband = ["short", "medium", "long"]
        band_check = False
        ch_check = False
        b_check = False
        if model.meta.instrument.channel in validch:
            ch_check = True
        if str(model.meta.instrument.band).lower() in validband:
            b_check = True

        if ch_check and b_check:
            band_check = True

        return band_check

    def _extract_ifu(self, model, exp_type, extract_ref, apcorr_ref):
        """
        Extract IFU spectra from a single datamodel.

        Parameters
        ----------
        model : DataModel
            Input model.
        exp_type : str
            Exposure type.
        extract_ref : str
            Path to the EXTRACT1D reference file or N/A.
        apcorr_ref : str
            Path to the APCORR reference file or N/A.

        Returns
        -------
        DataModel
            The output spectra.
        """
        source_type = model.meta.target.source_type
        if self.ifu_set_srctype is not None and exp_type == "MIR_MRS":
            source_type = self.ifu_set_srctype
            self.log.info(f"Overriding source type and setting it to {self.ifu_set_srctype}")

        if exp_type == "MIR_MRS":
            band_cube = self._check_mrs_type(model)
            if self.ifu_rfcorr and not band_cube:
                message = (
                    "Turning off residual fringe correction for MIRI MRS data "
                    "because the input is not a single IFU band"
                )
                self.log.info(message)
                self.ifu_rfcorr = False
        else:
            self.ifu_rfcorr = False
            self.log.info(
                "Turning off residual fringe correction because it only works on MIRI MRS data"
            )

        result = ifu_extract1d(
            model,
            extract_ref,
            source_type,
            self.subtract_background,
            self.bkg_sigma_clip,
            apcorr_ref,
            self.center_xy,
            self.ifu_autocen,
            self.ifu_rfcorr,
            self.ifu_rscale,
            self.ifu_covar_scale,
        )
        return result

    def _save_intermediate(self, intermediate_model, suffix, idx):
        """
        Save an intermediate output file.

        Parameters
        ----------
        intermediate_model : DataModel
            A model to save
        suffix : str
            Suffix to append to the output filename.
        idx : int or None
            Index to append to the output filename.
        """
        if isinstance(intermediate_model, ModelContainer):
            # Save the profile with the slit name + index + suffix 'profile'
            for model in intermediate_model:
                slit = str(model.name).lower()
                if idx is not None:
                    complete_suffix = f"{slit}_{idx}_{suffix}"
                else:
                    complete_suffix = f"{slit}_{suffix}"
                output_path = self.make_output_path(suffix=complete_suffix)
                self.log.info(f"Saving {suffix} {output_path}")
                model.save(output_path)
        else:
            # Only one profile - just use the index and suffix 'profile'
            if idx is not None:
                complete_suffix = f"{idx}_{suffix}"
            else:
                complete_suffix = suffix
            output_path = self.make_output_path(suffix=complete_suffix)
            self.log.info(f"Saving {suffix} {output_path}")
            intermediate_model.save(output_path)
        intermediate_model.close()

    def process(self, input_data):
        """
        Execute the step.

        Parameters
        ----------
        input_data : DataModel
            The input model.

        Returns
        -------
        DataModel
            This will be `input_model` if the step was skipped; otherwise,
            it will be a model containing 1-D extracted spectra.
        """
        # Open the input and figure out what type of model it is
        if isinstance(input_data, ModelContainer):
            input_model = input_data
        else:
            input_model = datamodels.open(input_data)

        if isinstance(
            input_model,
            (
                datamodels.CubeModel,
                datamodels.ImageModel,
                datamodels.SlitModel,
                datamodels.IFUCubeModel,
                ModelContainer,
                SourceModelContainer,
            ),
        ):
            # Acceptable input type, just log it
            self.log.debug(f"Input is a {str(input_model)}.")
        elif isinstance(input_model, datamodels.MultiSlitModel):
            # If input is multislit, with 3D calints, skip the step
            self.log.debug("Input is a MultiSlitModel")
            if len((input_model[0]).shape) == 3:
                self.log.warning("3D input is unsupported; step will be skipped")
                result = input_model.copy()
                result.meta.cal_step.extract_1d = "SKIPPED"
                return result
        else:
            self.log.error(f"Input is a {str(input_model)}, ")
            self.log.error("which was not expected for extract_1d.")
            self.log.error("The extract_1d step will be skipped.")
            result = input_model.copy()
            result.meta.cal_step.extract_1d = "SKIPPED"
            return result

        if not isinstance(input_model, ModelContainer):
            exp_type = input_model.meta.exposure.type

            # Make the input iterable
            input_model = [input_model]
        else:
            exp_type = input_model[0].meta.exposure.type
        self.log.debug(f"Input for EXP_TYPE {exp_type} contains {len(input_model)} items")

        if len(input_model) > 1 and exp_type in extract.WFSS_EXPTYPES:
            # For WFSS level-3, the input is a single entry of a
            # SourceContainer, which contains a list of multiple
            # SlitModels for a single source. Send the whole container
            # into extract1d and put all results in a single product.
            input_model = [input_model]

        if exp_type == "NIS_SOSS":
            # Data is NIRISS SOSS observation, use its own extraction routines
            self.log.info(
                "Input is a NIRISS SOSS observation, the specialized SOSS "
                "extraction (ATOCA) will be used."
            )

            # There is only one input model for this mode
            model = input_model[0]
            result = self._extract_soss(model)

        else:
            result = ModelContainer()
            for i, model in enumerate(input_model):
                # Get the reference file names
                extract_ref, apcorr_ref, psf_ref = self._get_extract_reference_files_by_mode(
                    model, exp_type
                )

                profile = None
                scene_model = None
                residual = None
                if isinstance(model, datamodels.IFUCubeModel):
                    # Call the IFU specific extraction routine
                    extracted = self._extract_ifu(model, exp_type, extract_ref, apcorr_ref)
                else:
                    # Call the general extraction routine
                    extracted, profile, scene_model, residual = extract.run_extract1d(
                        model,
                        extract_ref,
                        apcorr_ref,
                        psf_ref,
                        self.extraction_type,
                        self.smoothing_length,
                        self.bkg_fit,
                        self.bkg_order,
                        self.log_increment,
                        self.subtract_background,
                        self.use_source_posn,
                        self.position_offset,
                        self.model_nod_pair,
                        self.optimize_psf_location,
                        self.save_profile,
                        self.save_scene_model,
                        self.save_residual_image,
                    )

                # Set the step flag to complete in each model
                extracted.meta.cal_step.extract_1d = "COMPLETE"
                result.append(extracted)
                del extracted

                # Save profile if needed
                if len(input_model) > 1:
                    idx = i
                else:
                    idx = None
                if self.save_profile and profile is not None:
                    self._save_intermediate(profile, "profile", idx)

                # Save model if needed
                if self.save_scene_model and scene_model is not None:
                    self._save_intermediate(scene_model, "scene_model", idx)

                # Save residual if needed
                if self.save_residual_image and residual is not None:
                    self._save_intermediate(residual, "residual", idx)

            # If only one result, return the model instead of the container
            if len(result) == 1:
                result = result[0]

        if (
            (self._input_filename is None)
            and (not isinstance(result, ModelContainer))
            and (result is not None)
            and (self.output_file is None)
        ):
            # Fix output file naming for WFSS multislitmodels that were originally generated
            # from SourceModelContainer
            self._input_filename = result.meta.filename

        # For WFSS, reorder the x1d product to save it in the flat format
        if exp_type in extract.WFSS_EXPTYPES:
            result = make_wfss_multiexposure(result)
            result.meta.cal_step.extract_1d = "COMPLETE"

        return result
