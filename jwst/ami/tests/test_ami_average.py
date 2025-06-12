"""Unit tests for ami_average module and step."""

import pytest
from jwst.ami.ami_average_step import AmiAverageStep
from stdatamodels.jwst.datamodels import ImageModel


def test_ami_average_deprecated(example_model):
    image_list = []
    for i in range(example_model.data.shape[0]):
        image_list.append(ImageModel(example_model.data[i, :, :]))

    with pytest.deprecated_call():
        AmiAverageStep.call(image_list)
