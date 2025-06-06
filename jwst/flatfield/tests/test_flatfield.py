import pytest
import numpy as np
from astropy.io import fits
from astropy.table import Table
from numpy.testing import assert_allclose
from stdatamodels.jwst import datamodels

from jwst.assign_wcs import AssignWcsStep
from jwst.assign_wcs.tests.test_nirspec import create_nirspec_ifu_file
from jwst.flatfield import FlatFieldStep
from jwst.flatfield.flat_field_step import NRS_IMAGING_MODES, NRS_SPEC_MODES


@pytest.mark.parametrize(
    "instrument, exptype",
    [
        ("NIRCAM", "NRC_IMAGE"),
        ("NIRCAM", "NRC_WFSS"),
        ("MIRI", "MIR_IMAGE"),
        ("MIRI", "MIR_LRS-FIXEDSLIT"),
        ("MIRI", "MIR_LRS-SLITLESS"),
        ("MIRI", "MIR_MRS"),
        ("NIRISS", "NIS_IMAGE"),
        ("NIRISS", "NIS_WFSS"),
        ("NIRISS", "NIS_SOSS"),
        ("NIRISS", "NIS_AMI"),
        ("FGS", "FGS_IMAGE"),
    ]
    + [("NIRSPEC", exptype) for exptype in NRS_IMAGING_MODES],
)
def test_flatfield_step_interface(instrument, exptype):
    """Test that the basic interface works for data requiring a FLAT reffile"""

    shape = (20, 20)

    data = datamodels.ImageModel(shape)
    data.meta.instrument.name = instrument
    data.meta.exposure.type = exptype
    data.meta.subarray.xstart = 1
    data.meta.subarray.ystart = 1
    data.meta.subarray.xsize = shape[1]
    data.meta.subarray.ysize = shape[0]

    flat = datamodels.FlatModel(shape)
    flat.meta.instrument.name = instrument
    flat.meta.subarray.xstart = 1
    flat.meta.subarray.ystart = 1
    flat.meta.subarray.xsize = shape[1]
    flat.meta.subarray.ysize = shape[0]
    flat.data += 1
    flat.data[0, 0] = np.nan
    flat.err = np.random.random(shape) * 0.05

    # override class attribute so only the `flat` type needs to be overridden
    # in the step call.  Otherwise CRDS calls will be made for the other 3
    # types of flat reference file not used in this test.
    previous_reference_file_types = FlatFieldStep.reference_file_types
    try:
        FlatFieldStep.reference_file_types = ["flat"]
        result = FlatFieldStep.call(data, override_flat=flat)
    finally:
        FlatFieldStep.reference_file_types = previous_reference_file_types

    # flat is unity, so data should match up except where invalid
    is_nan = np.isnan(result.data)
    assert np.all(result.data[~is_nan] == data.data[~is_nan])

    # where data is invalid, error and variance are also invalid
    # and DQ is set to DO_NOT_USE
    assert np.all(np.isnan(result.err[is_nan]))
    assert np.all(np.isnan(result.var_rnoise[is_nan]))
    assert np.all(np.isnan(result.var_poisson[is_nan]))
    assert np.all(np.isnan(result.var_flat[is_nan]))
    assert np.all(result.dq[is_nan] | datamodels.dqflags.pixel["DO_NOT_USE"])

    assert result.var_flat.shape == shape
    assert result.meta.cal_step.flat_field == "COMPLETE"


def exptypes():
    """Generate NRS EXPTYPES from the schema enum, removing spec types"""
    model = datamodels.ImageModel()
    alltypes = set(model.meta.exposure._schema["properties"]["type"]["enum"])
    spectypes = set(NRS_SPEC_MODES)
    return sorted([i for i in (alltypes - spectypes)])


@pytest.mark.parametrize("exptype", exptypes())
def test_nirspec_flatfield_step_interface(exptype):
    """Test that the interface works all NIRSpec types"""

    shape = (20, 20)

    data = datamodels.ImageModel(shape)
    data.meta.observation.date = "2019-01-01"
    data.meta.observation.time = "00:00:00"
    data.meta.instrument.name = "NIRSPEC"
    data.meta.instrument.detector = "NRS1"
    data.meta.instrument.filter = "CLEAR"
    data.meta.instrument.grating = "MIRROR"
    data.meta.exposure.type = exptype
    data.meta.subarray.xstart = 1
    data.meta.subarray.ystart = 1
    data.meta.subarray.xsize = shape[1]
    data.meta.subarray.ysize = shape[0]

    FlatFieldStep.call(data)


def create_nirspec_flats(shape, msa=False, flat_data_value=1.0, flat_err_value=0.1):
    flats = []
    for flat_name in ["f", "s", "d"]:
        if flat_name == "f":
            f_data = np.full(shape[0], flat_data_value)
            f_err = np.full(shape[0], flat_err_value)
        else:
            f_data = np.full(shape[0], 1.0)
            f_err = np.full(shape[0], np.nan)

        # make a fast variation table
        flat_table = Table(
            {
                "slit_name": ["ANY"],
                "nelem": [shape[0]],
                "wavelength": [np.arange(1, shape[0] + 1, dtype=float)],
                "data": [f_data],
                "error": [f_err],
            }
        )

        fflat_hdul = fits.HDUList([fits.PrimaryHDU()])
        if flat_name == "s" and not msa:
            fflat_hdul.append(fits.ImageHDU(data=np.full(shape[-2:], 1.0), name="SCI"))
            fflat_hdul.append(fits.ImageHDU(data=np.full(shape[-2:], 0), name="DQ"))
            fflat_hdul.append(fits.ImageHDU(data=np.full(shape[-2:], 0.0), name="ERR"))

        elif flat_name == "s" and msa:
            fflat_hdul.append(fits.ImageHDU(data=np.full((5, shape[1], shape[2]), 1.0), name="SCI"))
            fflat_hdul.append(fits.ImageHDU(data=np.full((5, shape[1], shape[2]), 0), name="DQ"))
            fflat_hdul.append(fits.ImageHDU(data=np.full((5, shape[1], shape[2]), 0.0), name="ERR"))

            fflat_hdul.append(fits.table_to_hdu(Table({"wavelength": np.arange(1, shape[0] + 1)})))
            fflat_hdul[-1].name = "WAVELENGTH"

        elif flat_name == "d":
            fflat_hdul.append(fits.ImageHDU(data=np.full(shape, 1.0), name="SCI"))
            fflat_hdul.append(fits.ImageHDU(data=np.full(shape, 0), name="DQ"))
            fflat_hdul.append(fits.ImageHDU(data=np.full(shape, 0.0), name="ERR"))
            fflat_hdul.append(fits.table_to_hdu(Table({"wavelength": np.arange(1, shape[0] + 1)})))
            fflat_hdul[-1].name = "WAVELENGTH"

        if flat_name == "f" and msa:
            for quadrant in range(4):
                fflat_hdul.append(
                    fits.ImageHDU(data=np.full(shape, 1.0), name="SCI", ver=(quadrant + 1))
                )
                fflat_hdul.append(
                    fits.ImageHDU(data=np.full(shape, 0), name="DQ", ver=(quadrant + 1))
                )
                fflat_hdul.append(
                    fits.ImageHDU(data=np.full(shape, 0.0), name="ERR", ver=(quadrant + 1))
                )
                hdu = fits.table_to_hdu(Table({"wavelength": [np.arange(1, shape[0] + 1)]}))
                hdu.header["EXTNAME"] = "WAVELENGTH"
                hdu.header["EXTVER"] = quadrant + 1
                fflat_hdul.append(hdu)

                hdu = fits.table_to_hdu(flat_table)
                hdu.header["EXTNAME"] = "FAST_VARIATION"
                hdu.header["EXTVER"] = quadrant + 1
                fflat_hdul.append(hdu)

            flat = datamodels.NirspecQuadFlatModel(fflat_hdul)
        else:
            fflat_hdul.append(fits.table_to_hdu(flat_table))
            fflat_hdul[-1].name = "FAST_VARIATION"

            flat = datamodels.NirspecFlatModel(fflat_hdul)

        fflat_hdul.close()

        flat.meta.instrument.name = "NIRSPEC"
        flat.meta.subarray.xstart = 1
        flat.meta.subarray.ystart = 1
        flat.meta.subarray.xsize = shape[1]
        flat.meta.subarray.ysize = shape[0]

        flats.append(flat)

    return flats


def test_nirspec_bots_flat():
    """Test that the interface works for NIRSpec BOTS data"""
    shape = (3, 20, 20)
    w_shape = (10, 20, 20)

    data = datamodels.SlitModel(shape)
    data.meta.instrument.name = "NIRSPEC"
    data.meta.exposure.type = "NRS_BRIGHTOBJ"
    data.meta.subarray.xstart = 1
    data.meta.subarray.ystart = 1
    data.meta.subarray.xsize = shape[1]
    data.meta.subarray.ysize = shape[0]
    data.xstart = 1
    data.ystart = 1
    data.xsize = shape[1]
    data.ysize = shape[0]
    data.data += 1
    data.wavelength = np.ones(shape[-2:])
    data.wavelength[:] = np.linspace(1, 5, shape[-1], dtype=float)

    flats = create_nirspec_flats(w_shape)
    result = FlatFieldStep.call(
        data,
        override_fflat=flats[0],
        override_sflat=flats[1],
        override_dflat=flats[2],
        override_flat="N/A",
    )

    # null flat, so data is the same, other than nan edges
    nn = ~np.isnan(result.data)
    assert_allclose(result.data[nn], data.data[nn])

    # check that NaNs match in every extension they should
    for ext in ["data", "err", "var_rnoise", "var_poisson", "var_flat"]:
        test_data = getattr(result, ext)
        assert np.all(np.isnan(test_data[~nn]))
    assert np.all(result.dq[~nn] | datamodels.dqflags.pixel["DO_NOT_USE"])

    # error is propagated from non-empty fflat error
    # (no other additive contribution, scale from data is 1.0)
    assert result.var_flat.shape == shape
    assert_allclose(result.var_flat[nn], 0.1**2)
    assert result.meta.cal_step.flat_field == "COMPLETE"

    result.close()
    for flat in flats:
        flat.close()


@pytest.mark.parametrize("srctype", ["POINT", "EXTENDED"])
def test_nirspec_fs_flat(srctype):
    """Test that the interface works for NIRSpec FS data."""
    shape = (20, 20)
    w_shape = (10, 20, 20)

    data = datamodels.MultiSlitModel()
    data.meta.instrument.name = "NIRSPEC"
    data.meta.exposure.type = "NRS_FIXEDSLIT"
    data.meta.subarray.xstart = 1
    data.meta.subarray.ystart = 1
    data.meta.subarray.xsize = shape[1]
    data.meta.subarray.ysize = shape[0]

    data.slits.append(datamodels.SlitModel(shape))
    data.slits[0].data = np.full(shape, 1.0)
    data.slits[0].dq = np.full(shape, 0)
    data.slits[0].var_poisson = np.full(shape, 0.01)
    data.slits[0].var_rnoise = np.full(shape, 0.01)
    data.slits[0].err = np.sqrt(data.slits[0].var_poisson + data.slits[0].var_rnoise)
    data.slits[0].data = np.full(shape, 1.0)
    data.slits[0].wavelength = np.ones(shape[-2:])
    data.slits[0].wavelength[:] = np.linspace(1, 5, shape[-1], dtype=float)
    data.slits[0].source_type = srctype
    data.slits[0].name = "S200A1"
    data.slits[0].xstart = 1
    data.slits[0].ystart = 1
    data.slits[0].xsize = shape[1]
    data.slits[0].ysize = shape[0]

    orig_data = data.copy()
    flat_val = 0.9
    flat_err = 0.1
    flats = create_nirspec_flats(w_shape, flat_data_value=flat_val, flat_err_value=flat_err)
    result = FlatFieldStep.call(
        data,
        override_fflat=flats[0],
        override_sflat=flats[1],
        override_dflat=flats[2],
        override_flat="N/A",
    )

    nn = ~np.isnan(result.slits[0].data)
    assert_allclose(result.slits[0].data[nn], orig_data.slits[0].data[nn] / flat_val, rtol=1e-6)

    # Recover original input before flatfield step.
    rt_data = FlatFieldStep.call(
        result,
        inverse=True,
        override_fflat=flats[0],
        override_sflat=flats[1],
        override_dflat=flats[2],
        override_flat="N/A",
    )
    assert_allclose(rt_data.slits[0].data[nn], orig_data.slits[0].data[nn])
    assert_allclose(rt_data.slits[0].var_poisson[nn], orig_data.slits[0].var_poisson[nn])
    assert_allclose(rt_data.slits[0].var_rnoise[nn], orig_data.slits[0].var_rnoise[nn])
    assert_allclose(rt_data.slits[0].var_flat[nn], 0)  # orig_input has no var_flat
    assert_allclose(rt_data.slits[0].err[nn], orig_data.slits[0].err[nn])

    # check that NaNs match in every extension they should
    for ext in ["data", "err", "var_rnoise", "var_poisson", "var_flat"]:
        test_data = getattr(result.slits[0], ext)
        assert np.all(np.isnan(test_data[~nn]))
    assert np.all(result.slits[0].dq[~nn] | datamodels.dqflags.pixel["DO_NOT_USE"])

    # error is propagated from non-empty fflat error
    # (no other additive contribution, scale from data is 1.0)
    assert result.slits[0].var_flat.shape == shape
    assert_allclose(
        result.slits[0].var_flat[nn], ((1 / flat_val) / flat_val * flat_err) ** 2, rtol=1e-6
    )
    assert result.meta.cal_step.flat_field == "COMPLETE"

    result.close()
    for flat in flats:
        flat.close()


def test_nirspec_msa_flat():
    """Test that the interface works for NIRSpec MSA data."""
    shape = (20, 20)
    w_shape = (10, 20, 20)

    data = datamodels.MultiSlitModel()
    data.meta.instrument.name = "NIRSPEC"
    data.meta.exposure.type = "NRS_MSASPEC"
    data.meta.subarray.xstart = 1
    data.meta.subarray.ystart = 1
    data.meta.subarray.xsize = shape[1]
    data.meta.subarray.ysize = shape[0]

    data.slits.append(datamodels.SlitModel(shape))
    data.slits[0].data = np.full(shape, 1.0)
    data.slits[0].dq = np.full(shape, 0)
    data.slits[0].err = np.full(shape, 0.0)
    data.slits[0].var_poisson = np.full(shape, 0.0)
    data.slits[0].var_rnoise = np.full(shape, 0.0)
    data.slits[0].data = np.full(shape, 1.0)
    data.slits[0].wavelength = np.ones(shape[-2:])
    data.slits[0].wavelength[:] = np.linspace(1, 5, shape[-1], dtype=float)
    data.slits[0].source_type = "UNKNOWN"
    data.slits[0].name = "11"
    data.slits[0].quadrant = 1
    data.slits[0].xcen = 10
    data.slits[0].ycen = 10
    data.slits[0].xstart = 1
    data.slits[0].ystart = 1
    data.slits[0].xsize = shape[1]
    data.slits[0].ysize = shape[0]

    flats = create_nirspec_flats(w_shape, msa=True)
    result = FlatFieldStep.call(
        data,
        override_fflat=flats[0],
        override_sflat=flats[1],
        override_dflat=flats[2],
        override_flat="N/A",
    )

    # null flat, so data is the same, other than nan edges
    nn = ~np.isnan(result.slits[0].data)
    assert_allclose(result.slits[0].data[nn], data.slits[0].data[nn])

    # check that NaNs match in every extension they should
    for ext in ["data", "err", "var_rnoise", "var_poisson", "var_flat"]:
        test_data = getattr(result.slits[0], ext)
        assert np.all(np.isnan(test_data[~nn]))
    assert np.all(result.slits[0].dq[~nn] | datamodels.dqflags.pixel["DO_NOT_USE"])

    # error is propagated from non-empty fflat error
    # (no other additive contribution, scale from data is 1.0)
    assert result.slits[0].var_flat.shape == shape
    assert_allclose(result.slits[0].var_flat[nn], 0.1**2)
    assert result.meta.cal_step.flat_field == "COMPLETE"

    result.close()
    for flat in flats:
        flat.close()


def test_nirspec_ifu_flat():
    """
    Test that the interface works for NIRSpec IFU data.

    Larger data and more WCS operations required for testing make
    this test a little longer than the others.
    """
    shape = (2048, 2048)
    w_shape = (10, 2048, 2048)

    # IFU mode requires WCS information, so make a more realistic model
    hdul = create_nirspec_ifu_file(
        grating="PRISM", filter="CLEAR", gwa_xtil=0.35986012, gwa_ytil=0.13448857, gwa_tilt=37.1
    )
    hdul["SCI"].data = np.ones(shape, dtype=float)

    data = datamodels.IFUImageModel(hdul)
    data = AssignWcsStep.call(data)

    flats = create_nirspec_flats(w_shape)
    result = FlatFieldStep.call(
        data,
        override_fflat=flats[0],
        override_sflat=flats[1],
        override_dflat=flats[2],
        override_flat="N/A",
    )

    # null flat, so data is the same, other than nan edges
    nn = ~np.isnan(result.data)
    assert_allclose(result.data[nn], data.data[nn])

    # check that NaNs match in every extension they should
    for ext in ["data", "err", "var_rnoise", "var_poisson", "var_flat"]:
        test_data = getattr(result, ext)
        assert np.all(np.isnan(test_data[~nn]))
    assert np.all(result.dq[~nn] | datamodels.dqflags.pixel["DO_NOT_USE"])

    # error is propagated from non-empty fflat error
    # (no other additive contribution, scale from data is 1.0)
    assert result.var_flat.shape == shape
    assert_allclose(result.var_flat[nn], 0.1**2)
    assert result.meta.cal_step.flat_field == "COMPLETE"

    result.close()
    for flat in flats:
        flat.close()
