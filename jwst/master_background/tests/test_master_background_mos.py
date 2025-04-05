import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table
from stdatamodels.jwst.datamodels import ImageModel, MultiSlitModel

from jwst.stpipe import query_step_status
from jwst.assign_wcs import AssignWcsStep
from jwst.extract_2d.extract_2d_step import Extract2dStep
from jwst.extract_2d.tests.test_nirspec import create_nirspec_hdul
from jwst.master_background import MasterBackgroundMosStep
from jwst.master_background import nirspec_utils
from jwst.pixel_replace import PixelReplaceStep
from jwst.resample import ResampleSpecStep
from jwst.extract_1d import Extract1dStep


def create_msa_hdul():
    # Four shutters open on MSA, the first three are background slits only
    # with the fourth having a source
    shutter_data = {
        "slitlet_id": [12, 13, 14, 100],
        "msa_metadata_id": [1, 1, 1, 1],
        "shutter_quadrant": [4, 4, 4, 4],
        "shutter_row": [10, 10, 10, 10],
        "shutter_column": [22, 23, 24, 25],
        "source_id": [0, 0, 0, 1],
        "background": ["Y", "Y", "Y", "N"],
        "shutter_state": ["OPEN", "OPEN", "OPEN", "OPEN"],
        "estimated_source_in_shutter_x": [np.nan, np.nan, np.nan, 0.18283921],
        "estimated_source_in_shutter_y": [np.nan, np.nan, np.nan, 0.31907734],
        "dither_point_index": [1, 1, 1, 1],
        "primary_source": ["N", "N", "N", "Y"],
        "fixed_slit": ["NONE", "NONE", "NONE", "NONE"],
    }

    source_data = {
        "program": [95065],
        "source_id": [2],
        "source_name": ["95065_2"],
        "alias": ["2123"],
        "ra": [53.15],
        "dec": [-27.81],
        "preimage_id": ["95065001_000"],
        "stellarity": [1.0],
    }

    shutter_table = Table(shutter_data)
    source_table = Table(source_data)

    hdul = fits.HDUList()
    hdul.append(fits.PrimaryHDU())
    hdul.append(fits.ImageHDU())
    hdul.append(fits.table_to_hdu(shutter_table))
    hdul.append(fits.table_to_hdu(source_table))
    hdul[2].name = "SHUTTER_INFO"
    hdul[3].name = "SOURCE_INFO"

    return hdul


@pytest.fixture
def nirspec_msa_rate(tmp_path):
    hdul = create_nirspec_hdul()
    hdul[0].header["MSAMETFL"] = str(tmp_path / "test_msa_01.fits")
    hdul[0].header["EFFEXPTM"] = 1.0
    hdul[0].header["DURATION"] = 1.0
    filename = str(tmp_path / "test_nrs_msa_rate.fits")
    hdul.writeto(filename, overwrite=True)
    hdul.close()
    return filename


@pytest.fixture
def nirspec_msa_metfl(tmp_path):
    hdul = create_msa_hdul()
    filename = str(tmp_path / "test_msa_01.fits")
    hdul.writeto(filename, overwrite=True)
    hdul.close()
    return filename


@pytest.fixture
def nirspec_msa_extracted2d(nirspec_msa_rate, nirspec_msa_metfl):
    model = ImageModel(nirspec_msa_rate)
    model = AssignWcsStep.call(model)
    model = Extract2dStep.call(model)
    return model


def mk_multispec(model):
    specs_model = MultiSlitModel()
    specs_model.update(model)
    slits = []
    for slit in model.slits:
        if nirspec_utils.is_background_msa_slit(slit):
            slits.append(slit)
    specs_model.slits.extend(slits)
    specs_model = PixelReplaceStep.call(specs_model)
    specs_model = ResampleSpecStep.call(specs_model)
    specs_model = Extract1dStep.call(specs_model)
    return specs_model


def test_master_background_mos(nirspec_msa_extracted2d):
    model = nirspec_msa_extracted2d

    # RuntimeWarning: invalid value encountered in divide (slit.data /= conversion)
    # RuntimeWarning: overflow encountered in multiply (photom var_flat)
    # RuntimeWarning: overflow encountered in square (flat_field var_flat)
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        result = MasterBackgroundMosStep.call(model)

    # Check that the master_background_mos step was run
    assert query_step_status(result, "master_background") == "COMPLETE"

    finite = np.isfinite(result.slits[-1].data) & np.isfinite(model.slits[-1].data)
    sci_orig = model.slits[-1].data[finite]
    sci_bkgsub = result.slits[-1].data[finite]

    # Check that a background was subtracted from the science data
    assert not np.allclose(sci_orig, sci_bkgsub)

    del model
    del result


def test_create_background_from_multispec(nirspec_msa_extracted2d):
    model = nirspec_msa_extracted2d

    # Insert outliers into one of the background spectra
    nypix = len(model.slits[0].data)
    nxpix = len(model.slits[0].data)
    model.slits[0].data[nypix // 2, nxpix // 2 - 1 : nxpix // 2 + 1] = 10

    specs_model = mk_multispec(model)

    # First check that we can make a master background from the inputs

    # Check that with sigma_clip=None, the outlier is retained
    master_background = nirspec_utils.create_background_from_multispec(specs_model, sigma_clip=None)
    assert np.any(master_background.spec[0].spec_table["surf_bright"] > 1)

    # Confirm that using a median_filter will filter out the outlier
    master_background = nirspec_utils.create_background_from_multispec(specs_model, median_kernel=4)
    assert np.allclose(master_background.spec[0].spec_table["surf_bright"], 1)

    # Confirm that using a sigma clipping when combining background spectra
    # removes the outlier
    master_background = nirspec_utils.create_background_from_multispec(specs_model, sigma_clip=3)
    assert np.allclose(master_background.spec[0].spec_table["surf_bright"], 1)

    del model
    del specs_model


def test_map_to_science_slits(nirspec_msa_extracted2d):
    model = nirspec_msa_extracted2d
    specs_model = mk_multispec(model)

    master_background = nirspec_utils.create_background_from_multispec(specs_model)

    # Check that the master background is expanded to the shape of the input slits
    mb_multislit = nirspec_utils.map_to_science_slits(model, master_background)
    assert mb_multislit.slits[0].data.shape == model.slits[0].data.shape

    # background should be all ones, but won't be expanded to populate the whole
    # 2D array, so any non-zero background pixels should have a value of 1
    slit_data = mb_multislit.slits[0].data
    nonzero = slit_data != 0
    assert np.allclose(slit_data[nonzero], 1)

    del model
    del specs_model


def test_apply_master_background(nirspec_msa_extracted2d):
    model = nirspec_msa_extracted2d
    specs_model = mk_multispec(model)

    master_background = nirspec_utils.create_background_from_multispec(specs_model)
    mb_multislit = nirspec_utils.map_to_science_slits(model, master_background)

    result = nirspec_utils.apply_master_background(model, mb_multislit, inverse=False)

    # where the background is applied to the science it should be 0 and elsewhere 1
    sci_data_orig = model.slits[-1].data
    sci_data_bkgsub = result.slits[-1].data
    diff = sci_data_orig - sci_data_bkgsub
    assert np.any(diff != 0)
    assert np.allclose(diff[diff != 0], 1)

    # Check inverse application
    result = nirspec_utils.apply_master_background(model, mb_multislit, inverse=True)

    # where the background is applied to the science it should be 0 and elsewhere 1
    sci_data_orig = model.slits[-1].data
    sci_data_bkgsub = result.slits[-1].data
    diff = sci_data_orig - sci_data_bkgsub

    # Background subtraction was inverted so the differences will be -1
    assert np.any(diff != 0)
    assert np.allclose(diff[diff != 0], -1)

    del model
    del result
    del specs_model
