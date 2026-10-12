import types
from unittest.mock import patch

import numpy as np
import pytest

from jwst.cube_build.slit_cube_build_step import SlitCubeBuildStep


# Dummy object to replace STDataModels / MultiExposureModel structures
class DummyExposure:
    def __init__(self, name="slit1", source_id=1, source_type="POINT", source_name="src1"):
        self.name = name
        self.source_id = source_id
        self.source_type = source_type
        self.source_name = source_name
        self.meta = types.SimpleNamespace(
            bunit_data="Jy/sr",
            bunit_err="Jy/sr",
            photometry=types.SimpleNamespace(pixelarea_steradians=1e-12, pixelarea_arcsecsq=0.01),
            wcs=types.SimpleNamespace(dummy_wcs=True),
        )


class DummyModel:
    def __init__(self, detector="NRS1", exposures=None):
        self.exposures = exposures or [DummyExposure()]
        self.meta = types.SimpleNamespace(
            instrument=types.SimpleNamespace(
                name="NIRSPEC", detector=detector, grating="G235H", filter="F170LP"
            )
        )

    def close(self):
        pass


def test_step_defaults():
    """Verify default step parameter initializations."""
    step = SlitCubeBuildStep()
    assert step.linear_wave is True
    assert step.scalexy == 0.0
    assert step.scalew == 0.0
    assert step.weighting == "volume"
    assert step.coord_system == "skyalign"
    assert step.suffix == "s3d"


def test_organize_source_aggregates_metadata():
    """Test that organize_source correctly extracts list attributes from a input model."""
    step = SlitCubeBuildStep()
    dummy_model = DummyModel(
        detector="NRS1",
        exposures=[
            DummyExposure(name="slit1", source_id=101),
            DummyExposure(name="slit2", source_id=101),
        ],
    )

    source_dict = step.organize_source(dummy_model)

    assert source_dict["number_slits"] == 2
    assert source_dict["bunit"] == "Jy/sr"
    assert source_dict["bunit_err"] == "Jy/sr"
    assert len(source_dict["slit_models"]) == 2


def test_process_parameter_normalization(monkeypatch):
    """Verify process method normalizes strings and adjusts even spaxel counts to odd values."""
    step = SlitCubeBuildStep()
    step.coord_system = "WORLD"
    step.weighting = "VOLUME"
    step.nspax_x = 10  # Even number should be adjusted to 11
    step.nspax_y = 8  # Even number should be adjusted to 9
    step.slit_frac = 1.5  # Invalid range should trigger fallback to 1.0

    # Short-circuit prepare_output to avoid file IO / JWST datamodel requirement
    monkeypatch.setattr(step, "prepare_output", lambda input_data: input_data)

    dummy_model = DummyModel()

    # Pass an unsupported instrument to exit execution early right after parameter checking
    dummy_model.meta.instrument.name = "MIRI"

    result = step.process(dummy_model)

    # Check updated step parameters
    assert step.coord_system == "skyalign"
    assert step.weighting == "volume"
    assert step.nspax_x == 11
    assert step.nspax_y == 9
    assert step.slit_frac == 1.0
    # Step returns input_data upon early exit
    assert result == dummy_model


def test_read_cubepars_missing_config_raises_error(monkeypatch):
    """Ensure read_cubepars raises a ValueError when grating/filter combo is not found."""
    step = SlitCubeBuildStep()

    class DummyTable:
        def __getitem__(self, item):
            return self

        def __eq__(self, other):
            return np.array([False, False])

        def __and__(self, other):
            return np.array([False, False])

        def __len__(self):
            return 0

    # Mock astropy.table.Table.read to return our empty dummy structure
    import astropy.table

    monkeypatch.setattr(astropy.table.Table, "read", lambda *args, **kwargs: DummyTable())

    with pytest.raises(ValueError, match="No matching configuration found"):
        step.read_cubepars("dummy_path.fits", "G235H", "F170LP")


def test_process_unsupported_instrument_miri(monkeypatch, caplog):
    """Verify process logs a warning and returns input_data if the instrument is MIRI."""
    step = SlitCubeBuildStep()
    dummy_model = DummyModel()
    dummy_model.meta.instrument.name = "MIRI"

    monkeypatch.setattr(step, "prepare_output", lambda input_data: input_data)

    # Mock isinstance to allow DummyModel through the MultiExposureModel check
    with patch("jwst.cube_build.slit_cube_build_step.isinstance", return_value=True):
        with caplog.at_level("WARNING"):
            result = step.process(dummy_model)

    assert result is dummy_model
    assert "Input instrument 'MIRI' is not supported" in caplog.text
    assert "SlitCubeBuildStep is designed specifically for NIRSpec" in caplog.text


def test_process_invalid_datamodel_type(monkeypatch, caplog):
    """Verify process logs a warning and returns input_data if not a MultiExposureModel."""
    step = SlitCubeBuildStep()
    dummy_model = DummyModel()

    monkeypatch.setattr(step, "prepare_output", lambda input_data: input_data)

    with caplog.at_level("WARNING"):
        result = step.process(dummy_model)

    assert result is dummy_model
    assert "SlitCubeBuildStep requires a MultiExposureModel" in caplog.text
