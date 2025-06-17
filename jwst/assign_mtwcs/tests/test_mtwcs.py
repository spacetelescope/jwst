import pytest
from astropy.utils.data import get_pkg_data_filename
from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelLibrary
from jwst.assign_mtwcs import AssignMTWcsStep


@pytest.mark.parametrize("errtype", ["ra", "dec", "both", "none"])
def test_mt_multislit(errtype):
    file_path = get_pkg_data_filename("data/test_mt_asn.json", package="jwst.assign_mtwcs.tests")
    with datamodels.open(file_path) as model:
        # Modify first model in library to have invalid data.
        if errtype in ("ra", "both"):
            del model._models[0].slits[0].meta.wcsinfo.mt_ra
        if errtype in ("dec", "both"):
            del model._models[0].slits[0].meta.wcsinfo.mt_dec
        assert model[0].slits[0].meta.wcs.output_frame.name == "world"

    step = AssignMTWcsStep()
    result = step.run(model)
    assert isinstance(result, ModelLibrary)

    if errtype == "none":
        expected_frame = "moving_target"
    else:
        expected_frame = "world"

    with result:
        zero = result.borrow(0)
        one = result.borrow(1)

        assert len(zero.slits) == 1
        assert zero.slits[0].meta.wcs.output_frame.name == expected_frame
        assert len(one.slits) == 1
        assert one.slits[0].meta.wcs.output_frame.name == expected_frame

        result.shelve(zero, 0, modify=False)
        result.shelve(one, 1, modify=False)


@pytest.mark.parametrize("errtype", ["ra", "dec", "both", "none"])
def test_mt_slitmodel(errtype):
    # borrow wcsinfo from multislitmodel
    file_path = get_pkg_data_filename("data/test_mt_asn.json", package="jwst.assign_mtwcs.tests")
    with datamodels.open(file_path) as model:
        wcs = model._models[0].slits[0].meta.wcs
        wcsinfo = model._models[0].slits[0].meta.wcsinfo

    step = AssignMTWcsStep()
    model = datamodels.SlitModel()
    model.meta.wcs = wcs
    model.meta.wcsinfo = wcsinfo
    model.meta.exposure.type = "NRS_FIXEDSLIT"

    if errtype in ("ra", "both"):
        del model.meta.wcsinfo.mt_ra
    if errtype in ("dec", "both"):
        del model.meta.wcsinfo.mt_dec

    library = ModelLibrary([model])
    result = step.run(library)

    if errtype == "none":
        expected_frame = "moving_target"
    else:
        expected_frame = "world"

    with result:
        zero = result.borrow(0)
        assert zero.meta.wcs.output_frame.name == expected_frame
        result.shelve(zero, 0, modify=False)
