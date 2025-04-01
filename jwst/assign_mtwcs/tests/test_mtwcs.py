from astropy.utils.data import get_pkg_data_filename
from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelLibrary
from jwst.assign_mtwcs import AssignMTWcsStep


def test_mt_multislit():
    file_path = get_pkg_data_filename(
        "data/test_mt_asn.json", package="jwst.assign_mtwcs.tests")
    with datamodels.open(file_path) as model:
        assert model[0].slits[0].meta.wcs.output_frame.name == 'world'
    step = AssignMTWcsStep()
    result = step.run(file_path)
    assert isinstance(result, ModelLibrary)
    with result:
        zero = result.borrow(0)
        one = result.borrow(1)

        assert len(zero.slits) == 1
        assert zero.slits[0].meta.wcs.output_frame.name == 'moving_target'
        assert len(one.slits) == 1
        assert one.slits[0].meta.wcs.output_frame.name == 'moving_target'

        result.shelve(zero, 0, modify=False)
        result.shelve(one, 1, modify=False)
