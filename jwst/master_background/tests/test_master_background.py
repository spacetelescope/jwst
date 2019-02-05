"""
Unit tests for master background subtraction
"""

import pytest

from jwst.master_background import MasterBackgroundStep
from jwst import datamodels


@pytest.mark.parametrize('input_data', [
    datamodels.ImageModel(),
    datamodels.MultiSlitModel(),
    datamodels.IFUImageModel(),
    datamodels.ModelContainer([datamodels.ImageModel()])
    ])
def test_master_background_init(input_data):
    """Verify an empty imagemodel can run through the step"""
    step = MasterBackgroundStep()
    result = step.run(input_data)

    assert type(input_data) is type(result)
