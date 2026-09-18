import logging

import numpy as np
from crds.core.exceptions import CrdsLookupError
from stdatamodels.jwst import datamodels

from jwst.combine_1d.combine_1d_step import Combine1dStep
from jwst.cube_build.cube_build_step import CubeBuildStep
from jwst.datamodels import ModelContainer
from jwst.extract_1d.extract_1d_step import Extract1dStep
from jwst.lib.basic_utils import disable_logging
from jwst.lib.reffile_utils import find_row
from jwst.residual_fringe.utils import fit_residual_fringes_1d
from jwst.spectral_leak.spectral_leak_step import SpectralLeakStep
from jwst.stpipe import Step, record_step_status

__all__ = ["PFPCStep"]

log = logging.getLogger(__name__)


class PFPCStep(Step):
    """Apply a point fixed pattern correction (PFPC) to dithered point source spectra."""

    class_alias = "pfpc"

    spec = """
        skip = boolean(default=True) # By default, skip the step.
        suffix = string(default='pfpc') # Set the output suffix.
        output_use_model = boolean(default=True) # Use input filenames in the output models
    """  # noqa: E501

    reference_file_types = ["pfpc"]

    def process(self, input_data):
        """
        Extract spectra and apply a point fixed pattern correction (PFPC).

        Parameters
        ----------
        input_data : str or `~stdatamodels.jwst.datamodels.JwstDataModel` \
                     or `~jwst.datamodels.container.ModelContainer`
            The input filename, datamodel, or container of 2D spectral data.
            If the input is a container, each contained model is separately
            processed.

        Returns
        -------
        output_model : `~jwst.datamodels.container.ModelContainer`
            Extracted spectra with fixed pattern correction applied. If
            no correction can be applied, an empty container is returned.
        """
        # Read in models. In pipeline context, these may be modified
        # in place to set a 'FAILED' status if the step does not succeed,
        # but they are not returned by the step.
        output_model = self.prepare_output(input_data)
        if isinstance(output_model, ModelContainer):
            models = output_model
        else:
            models = [output_model]

        # Get output file name from self or parent if available
        output_file = self.search_attr("output_file")
        if output_file is None:
            # Set up output path name to include the ASN ID if available
            self.add_asn_id_to_output_name(models)

        # Make an empty container for output: if no products can be made,
        # none are returned
        output_container = ModelContainer()

        # Get the reference file from the first model
        try:
            pfpc_file = self.get_reference_file(models[0], "pfpc")
        except CrdsLookupError:
            pfpc_file = "N/A"
        if pfpc_file == "N/A":
            log.warning("No PFPC reference file found.")
            record_step_status(output_model, "pfpc", success=False)
            return output_container

        with datamodels.open(pfpc_file) as pfpc_model:
            pfpc_table = pfpc_model.pfpc_table

        # Check the first input against the PFPC table and skip processing if no match
        ignore_keys = ["channel"]  # channel may not exactly match before extraction
        first_correction = self._find_correction(
            models[0], pfpc_table, ignore_keys=ignore_keys, require_one=False
        )
        if first_correction is None:
            log.warning("No matching correction found for input models.")
            record_step_status(output_model, "pfpc", success=False)
            return output_container

        log.info("Building cubes and extracting spectra from each exposure")
        spec_by_dither = {}
        for model in models:
            log.info(f"Working on exposure {model.meta.filename}")

            # TODO: check input exptype and process accordingly

            cube_param = {
                "coord_system": "ifualign",
                "output_type": "band",
                "output_file": output_file,
            }
            log.debug(f"Calling the cube_build step with parameters {cube_param}")
            with disable_logging(level=logging.INFO):
                cubes = CubeBuildStep.call(model, **cube_param)

            extract_param = {
                "ifu_autocen": True,
                "ifu_rfcorr": False,
            }
            log.debug(
                f"Calling the extract_1d step to extract spectra with parameters {extract_param}"
            )
            with disable_logging(level=logging.INFO):
                spectra = Extract1dStep.call(cubes, **extract_param)

            patt_num = model.meta.dither.position_number
            if patt_num in spec_by_dither:
                spec_by_dither[patt_num].extend(spectra)
            else:
                spec_by_dither[patt_num] = spectra

        log.info("Correcting for spectral leak at each dither position")
        for patt_num in spec_by_dither:
            with disable_logging(level=logging.INFO):
                spec_by_dither[patt_num] = SpectralLeakStep.call(spec_by_dither[patt_num])

        log.info("Applying PFPC correction to each spectrum")
        spec_by_band = {}
        spec_columns = [
            "FLUX",
            "FLUX_ERROR",
            # "FLUX_VAR_POISSON",
            # "FLUX_VAR_RNOISE",
            # "FLUX_VAR_FLAT",
            "SURF_BRIGHT",
            "SB_ERROR",
            # "SB_VAR_POISSON",
            # "SB_VAR_RNOISE",
            # "SB_VAR_FLAT",
        ]
        for patt_num in spec_by_dither:
            for spec in spec_by_dither[patt_num]:
                table_row = self._find_correction(spec, pfpc_table)
                if table_row is None:
                    log.warning(f"No matching correction found for {spec.meta.filename}")
                    continue

                pfpc_wave = table_row["wavelength"]
                pfpc_corr = table_row["correction"]
                valid_pfpc = np.isfinite(pfpc_wave) & np.isfinite(pfpc_corr)

                spec_wave = spec.spec[0].spec_table["WAVELENGTH"]
                valid_spec = np.isfinite(spec_wave)
                if np.any(valid_pfpc) and np.any(valid_spec):
                    interp_correction = np.interp(
                        spec_wave[valid_spec],
                        pfpc_wave[valid_pfpc],
                        pfpc_corr[valid_pfpc],
                        left=np.nan,
                        right=np.nan,
                    )
                    interp_correction[interp_correction == 0] = np.nan
                else:
                    log.warning(f"No valid correction for {spec.meta.filename}")
                    continue

                # Update spectral flux and variances with the correction
                for column in spec_columns:
                    data = spec.spec[0].spec_table[column]
                    corrected_data = np.full_like(data, np.nan)
                    # if "VAR" in column:
                    #    corrected_data[valid_spec] = data[valid_spec] / interp_correction**2
                    # else:
                    #    corrected_data[valid_spec] = data[valid_spec] / interp_correction
                    corrected_data[valid_spec] = data[valid_spec] / interp_correction
                    spec.spec[0].spec_table[column] = corrected_data

                # Re-sort spectra by channel/band for combining across dither positions
                channel = str(spec.meta.instrument.channel)
                band = str(spec.meta.instrument.band).lower()
                channel_band = f"ch{channel}-{band}"
                if channel_band not in spec_by_band:
                    spec_by_band[channel_band] = [spec]
                else:
                    spec_by_band[channel_band].append(spec)

        log.info("Combining spectra across all dithers for each band")
        combine_param = {"sigma_clip": 4.0}
        log.debug(f"Calling the combine_1d step with parameters {combine_param}")
        for band in spec_by_band:
            # TODO: combine1d does not propagate variance or background
            with disable_logging(level=logging.WARNING):
                combined = Combine1dStep.call(spec_by_band[band], **combine_param)

            # Retrieve data from combined spectrum model
            flux = combined.spec[0].spec_table["FLUX"]
            wave = combined.spec[0].spec_table["WAVELENGTH"]
            sb = combined.spec[0].spec_table["SURF_BRIGHT"]

            # Run residual fringing on combined spectrum
            channel = int(combined.meta.instrument.channel)
            log.debug("Calling fit_residual_fringes_1d to defringe")
            with disable_logging(level=logging.WARNING):
                defringe_flux = fit_residual_fringes_1d(flux, wave, channel=channel)
                defringe_sb = fit_residual_fringes_1d(sb, wave, channel=channel)

            # Reassemble combined spectrum into a multispecmodel
            spec = datamodels.MRSSpecModel((flux.size,))
            spec.spec_table["WAVELENGTH"] = wave
            spec.spec_table["FLUX"] = flux
            spec.spec_table["FLUX_ERROR"] = combined.spec[0].spec_table["ERROR"]
            spec.spec_table["SURF_BRIGHT"] = sb
            spec.spec_table["SB_ERROR"] = combined.spec[0].spec_table["SB_ERROR"]
            spec.spec_table["DQ"] = combined.spec[0].spec_table["DQ"]
            spec.spec_table["RF_FLUX"] = defringe_flux
            spec.spec_table["RF_SURF_BRIGHT"] = defringe_sb

            multispec = datamodels.MRSMultiSpecModel()
            multispec.spec.append(spec)
            multispec.update(spec_by_band[band][0])
            multispec.meta.cal_step.pfpc = "COMPLETE"

            output_container.append(multispec)

        # TODO - close intermediate models

        # Close the input models if necessary: returned models are newly created
        if output_model is not input_data:
            output_model.close()

        return output_container

    @staticmethod
    def _find_correction(model, pfpc_table, ignore_keys=None, require_one=True):
        pfpc_ignore = ["wavelength", "correction"]
        if ignore_keys is not None:
            pfpc_ignore.extend(ignore_keys)

        # Set up a default dictionary with keywords to match,
        # from fields present in the PFPC table
        fields_to_match = {}
        for col_name in pfpc_table.columns.names:
            if col_name.lower() not in pfpc_ignore:
                fields_to_match[col_name] = "N/A"

        # Read model metadata into a flat dict
        model_meta = model.to_flat_dict(include_arrays=False)

        # Find the metadata corresponding to the FITS keyword
        # in the table column names
        for field in fields_to_match:
            meta_field = model.find_fits_keyword(field.upper())

            # Update match values if metadata is present
            if len(meta_field) == 1 and meta_field[0] in model_meta:
                fields_to_match[field] = model_meta[meta_field[0]]

        # Find a matching row in the table from all fields.
        # Allow it to raise an error for multiple rows if require_one is True.
        # Will warn and return None if no match is found.
        table_row = find_row(pfpc_table, fields_to_match, require_one=require_one)

        return table_row
