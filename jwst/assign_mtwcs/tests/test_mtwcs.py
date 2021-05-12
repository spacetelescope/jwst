import os

from jwst import datamodels
from jwst.assign_mtwcs import AssignMTWcsStep
from jwst.assign_mtwcs.tests import data

data_path = os.path.split(os.path.abspath(data.__file__))[0]


def test_mt_multislit():
    file_path = os.path.join(data.__path__[0], 'test_mt_asn.json')
    with datamodels.open(file_path) as model:
        assert model[0].slits[0].meta.wcs.output_frame.name == 'world'
        step = AssignMTWcsStep()
        result = step.run(model)
    assert isinstance(result, datamodels.ModelContainer)
    assert len(result[0].slits) == 1
    assert result[0].slits[0].meta.wcs.output_frame.name == 'moving_target'
    assert len(result[1].slits) == 1
    assert result[1].slits[0].meta.wcs.output_frame.name == 'moving_target'
