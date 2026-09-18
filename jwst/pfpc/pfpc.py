import logging

import numpy as np
from stdatamodels.jwst import datamodels

from jwst.combine_1d.combine_1d_step import Combine1dStep
from jwst.cube_build.cube_build_step import CubeBuildStep
from jwst.datamodels import ModelContainer
from jwst.extract_1d.extract_1d_step import Extract1dStep
from jwst.lib.basic_utils import disable_logging
from jwst.lib.reffile_utils import find_row
from jwst.residual_fringe.utils import fit_residual_fringes_1d
from jwst.spectral_leak.spectral_leak_step import SpectralLeakStep

log = logging.getLogger(__name__)

__all__ = ["find_correction", "process_exposures", "apply_correction", "combine_dithers"]

SUPPORTED_EXPTYPES = ["MIR_MRS", "NRS_IFU"]


def find_correction(model, pfpc_table, ignore_keys=None, require_one=True):
    """
    Find a matching correction row in the PFPC table.

    Parameters
    ----------
    model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        Input datamodel.
    pfpc_table : `~astropy.io.fits.fitsrec.FITS_rec`
        Record array containing PFPC corrections and columns to match to the data.
    ignore_keys : list of str or None, optional
        List of columns in the input table to ignore for matching purposes.
    require_one : bool, optional
        If True, an error is raised if more than one row matches the
        input model.

    Returns
    -------
    table_row : ndarray
        The matching row, containing PFPC corrections and corresponding
        wavelengths

    Raises
    ------
    MatchRowError
        When more than one row matches and ``require_one`` is True.
    """
    pfpc_ignore = ["wavelength", "correction"]
    if ignore_keys is not None:
        pfpc_ignore.extend(ignore_keys)

    # Set up a default dictionary with keywords to match,
    # from fields present in the PFPC table
    fields_to_match = {}
    for col_name in pfpc_table.columns.names:
        # TODO: add reference file handling
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


def process_exposures(models, output_file=None):
    """
    Extract spectra from each exposure and band.

    Currently only IFU modes are supported.  Calls ``cube_build`` then
    ``extract_1d`` with standardized parameters.  For MIRI MRS,
    also calls ``spectral_leak``.

    Parameters
    ----------
    models : list or `~jwst.datamodels.container.ModelContainer`
        Input datamodels.
    output_file : str or None, optional
        Output file name, from the step or parent parameters. Passed
        to cube_build to compose the output name.

    Returns
    -------
    spectra : list of `~stdatamodels.jwst.datamodels.JwstDataModel`
        Extracted spectra, one for each exposure and band.
    """
    log.info("Extracting spectra from each exposure")
    spec_by_dither = {}
    is_mrs = False
    for model in models:
        log.info(f"Working on exposure {model.meta.filename}")

        exp_type = str(model.meta.exposure.type).upper()
        if exp_type == "MIR_MRS":
            is_mrs = True
        if exp_type not in SUPPORTED_EXPTYPES:
            raise TypeError(f"Exposure type {model.meta.exposure.type} is not supported")

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
        log.debug(f"Calling the extract_1d step to extract spectra with parameters {extract_param}")
        with disable_logging(level=logging.INFO):
            spectra = Extract1dStep.call(cubes, **extract_param)

        if not isinstance(spectra, ModelContainer):
            spectra = [spectra]
        patt_num = model.meta.dither.position_number
        if patt_num in spec_by_dither:
            spec_by_dither[patt_num].extend(spectra)
        else:
            spec_by_dither[patt_num] = spectra

        # TODO - close intermediate models

    # Do spectral leak correction if needed and assemble a list of spectra to return
    output_spectra = []
    for patt_num in spec_by_dither:
        if is_mrs:
            log.info(f"Correcting for spectral leak at dither position {patt_num}")
            with disable_logging(level=logging.INFO):
                spec_by_dither[patt_num] = SpectralLeakStep.call(spec_by_dither[patt_num])

        output_spectra.extend(spec_by_dither[patt_num])

    return output_spectra


def apply_correction(spec_by_exposure, pfpc_table):
    """
    Apply the PFPC correction to each spectrum.

    Parameters
    ----------
    spec_by_exposure : list of `~stdatamodels.jwst.datamodels.JwstDataModel`
        Extracted spectra, one for each exposure and band.
    pfpc_table : `~astropy.io.fits.fitsrec.FITS_rec`
        Record array containing PFPC corrections and columns to match to the data.

    Returns
    -------
    corrected_spec : list of `~stdatamodels.jwst.datamodels.JwstDataModel`
        Corrected spectra.
    """
    log.info("Applying PFPC correction to each spectrum")
    corrected_spec = []
    spec_columns = [
        "FLUX",
        "FLUX_ERROR",
        "SURF_BRIGHT",
        "SB_ERROR",
    ]
    for spec in spec_by_exposure:
        table_row = find_correction(spec, pfpc_table)
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
            interp_correction = None

        if interp_correction is None or np.all(np.isnan(interp_correction)):
            log.warning(f"No valid correction for {spec.meta.filename}")
            continue

        # Update spectral flux and error with the correction
        for column in spec_columns:
            data = spec.spec[0].spec_table[column]
            corrected_data = np.full_like(data, np.nan)
            corrected_data[valid_spec] = data[valid_spec] / interp_correction
            spec.spec[0].spec_table[column] = corrected_data

        # Keep only the corrected spectra
        corrected_spec.append(spec)

    return corrected_spec


def combine_dithers(corrected_spec):
    """
    Combine PFPC corrected spectra across dither positions.

    Parameters
    ----------
    corrected_spec : list of `~stdatamodels.jwst.datamodels.JwstDataModel`
        PFPC corrected spectra from individual exposures and bands.

    Returns
    -------
    combined_spec : list of `~stdatamodels.jwst.datamodels.JwstDataModel`
        Combined spectra.
    """
    # Sort spectra by band
    spec_by_band = {}
    is_mrs = False
    for spec in corrected_spec:
        if spec.meta.exposure.type == "MIR_MRS":
            is_mrs = True
            channel = str(spec.meta.instrument.channel)
            band = str(spec.meta.instrument.band).lower()
            key = f"ch{channel}-{band}"
        else:
            filt = str(spec.meta.instrument.filter).lower()
            grating = str(spec.meta.instrument.grating).lower()
            key = f"{grating}-{filt}"
        if key not in spec_by_band:
            spec_by_band[key] = ModelContainer([spec])
        else:
            spec_by_band[key].append(spec)

    log.info("Combining spectra across all dithers for each band")
    combine_param = {"sigma_clip": 4.0}
    log.debug(f"Calling the combine_1d step with parameters {combine_param}")
    combined_spec = []
    for band in spec_by_band:
        # TODO: combine1d does not propagate variance or background
        with disable_logging(level=logging.WARNING):
            combined = Combine1dStep.call(spec_by_band[band], **combine_param)

        if not isinstance(combined, datamodels.MultiCombinedSpecModel):
            continue

        # Retrieve data from combined spectrum model
        flux = combined.spec[0].spec_table["FLUX"]
        wave = combined.spec[0].spec_table["WAVELENGTH"]
        sb = combined.spec[0].spec_table["SURF_BRIGHT"]

        # Run residual fringing on combined spectrum
        defringe_flux, defringe_sb = None, None
        if is_mrs:
            channel = int(combined.meta.instrument.channel)
            log.debug("Calling fit_residual_fringes_1d to defringe")
            with disable_logging(level=logging.WARNING):
                defringe_flux = fit_residual_fringes_1d(flux, wave, channel=channel)
                defringe_sb = fit_residual_fringes_1d(sb, wave, channel=channel)

        # Reassemble combined spectrum into a multispecmodel
        if is_mrs:
            spec = datamodels.MRSSpecModel((flux.size,))
            spec.spec_table["RF_FLUX"] = defringe_flux
            spec.spec_table["RF_SURF_BRIGHT"] = defringe_sb
        else:
            spec = datamodels.SpecModel((flux.size,))
        spec.spec_table["WAVELENGTH"] = wave
        spec.spec_table["FLUX"] = flux
        spec.spec_table["FLUX_ERROR"] = combined.spec[0].spec_table["ERROR"]
        spec.spec_table["SURF_BRIGHT"] = sb
        spec.spec_table["SB_ERROR"] = combined.spec[0].spec_table["SB_ERROR"]
        spec.spec_table["DQ"] = combined.spec[0].spec_table["DQ"]

        if is_mrs:
            multispec = datamodels.MRSMultiSpecModel()
        else:
            multispec = datamodels.MultiSpecModel()
        multispec.spec.append(spec)
        multispec.update(spec_by_band[band][0])

        combined_spec.append(multispec)

    return combined_spec
