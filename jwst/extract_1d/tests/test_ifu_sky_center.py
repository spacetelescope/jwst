import numpy as np
import pytest
from stdatamodels.jwst import datamodels

from jwst.extract_1d import extract_1d_step


def test_center_radec_is_transformed_to_cube_pixels(monkeypatch):
    model = datamodels.IFUCubeModel((2, 3, 4))

    monkeypatch.setattr(extract_1d_step, "locn_from_wcs", lambda model, ra, dec: (12.5, 13.5))

    center = extract_1d_step.Extract1dStep._resolve_ifu_center(
        model, None, (123.4, -54.3)
    )

    assert center == (12.5, 13.5)
    model.close()


def test_center_xy_and_center_radec_are_mutually_exclusive():
    model = datamodels.IFUCubeModel((2, 3, 4))

    with pytest.raises(ValueError, match="cannot both be specified"):
        extract_1d_step.Extract1dStep._resolve_ifu_center(
            model, (10.0, 11.0), (123.4, -54.3)
        )

    model.close()


def test_invalid_sky_transform_is_rejected(monkeypatch):
    model = datamodels.IFUCubeModel((2, 3, 4))
    monkeypatch.setattr(extract_1d_step, "locn_from_wcs", lambda model, ra, dec: (np.nan, 4.0))

    with pytest.raises(ValueError, match="valid IFU cube position"):
        extract_1d_step.Extract1dStep._resolve_ifu_center(
            model, None, (123.4, -54.3)
        )

    model.close()
