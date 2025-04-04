"""Unit tests for ami_average module and step."""

import pytest
from jwst.ami.ami_average_step import AmiAverageStep


def test_ami_average_deprecated(example_model):
    with pytest.deprecated_call():
        AmiAverageStep.call(example_model)
