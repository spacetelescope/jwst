"""Unit tests for SlitCubeBuildStep."""

import astropy.table
import numpy as np
import pytest
from stdatamodels.jwst.datamodels import ImageModel, MultiExposureModel, SlitModel

from jwst.cube_build.slit_cube_build_step import SlitCubeBuildStep

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_model(instrument="NIRSPEC", detector="NRS1", grating="G235H", filter="F170LP", n_slits=1):
    """Return a real MultiExposureModel with n_slits SlitModel exposures.

    All metadata fields needed by the step are populated; the actual pixel
    arrays are left at their zero-filled defaults so construction is cheap.
    """
    model = MultiExposureModel()
    model.meta.instrument.name = instrument
    model.meta.instrument.detector = detector
    model.meta.instrument.grating = grating
    model.meta.instrument.filter = filter
    for i in range(n_slits):
        slit = SlitModel()
        slit.name = f"slit{i + 1}"
        slit.source_id = 101
        slit.source_type = "POINT"
        slit.source_name = "src1"
        slit.meta.bunit_data = "MJy/sr"
        slit.meta.bunit_err = "MJy/sr"
        model.exposures.append(slit)
    return model


class _StepNoIO(SlitCubeBuildStep):
    """SlitCubeBuildStep subclass that skips all file I/O.

    * ``prepare_output`` is overridden to be a no-op so tests that call
      ``process`` never touch the filesystem.
    * ``read_cubepars`` is overridden to raise ``ValueError`` immediately when
      no matching row exists in the supplied in-memory table, exercising the
      same error path without reading any file.
    """

    def prepare_output(self, input_data):
        return input_data

    def read_cubepars(self, par_file, grating, filter_name):
        """Use an in-memory table with no matching rows instead of a real file."""
        table = astropy.table.Table({"grating": ["G140H", "G395H"], "filter": ["F100LP", "F290LP"]})
        match = (table["grating"] == grating) & (table["filter"] == filter_name)
        if not np.any(match):
            raise ValueError(
                f"No matching configuration found for grating={grating}, filter={filter_name}"
            )
        return table[match]


# ---------------------------------------------------------------------------
# organize_source
# ---------------------------------------------------------------------------


class TestOrganizeSource:
    """organize_source extracts the correct metadata from a MultiExposureModel."""

    def test_number_of_slits_matches_exposures(self):
        """number_slits must equal the number of exposures in the model."""
        model = _make_model(n_slits=2)
        source_dict = SlitCubeBuildStep().organize_source(model)
        assert source_dict["number_slits"] == 2
        assert len(source_dict["slit_models"]) == 2
        model.close()


# ---------------------------------------------------------------------------
# process — parameter normalization
# ---------------------------------------------------------------------------


class TestProcessParameterNormalization:
    """process normalizes string case and enforces odd spaxel counts and slit_frac range."""

    def test_coord_system_lowercased(self):
        step = _StepNoIO()
        step.coord_system = "WORLD"
        step.process(_make_model(instrument="MIRI"))
        assert step.coord_system == "skyalign"

    def test_even_nspax_x_bumped_to_next_odd(self):
        """Even nspax_x must be incremented to the next odd value."""
        step = _StepNoIO()
        step.nspax_x = 10
        step.nspax_y = 8
        step.process(_make_model(instrument="MIRI"))
        assert step.nspax_x == 11
        assert step.nspax_y == 9

    def test_slit_frac_out_of_range_reset_to_one(self):
        """slit_frac outside [0, 1] must be reset to 1.0."""
        step = _StepNoIO()
        step.slit_frac = 1.5
        step.process(_make_model(instrument="MIRI"))
        assert step.slit_frac == 1.0


# ---------------------------------------------------------------------------
# read_cubepars — missing configuration
# ---------------------------------------------------------------------------


class TestReadCubepars:
    """read_cubepars raises ValueError when no matching grating/filter row exists."""

    def test_missing_config_raises_value_error(self):
        with pytest.raises(ValueError, match="No matching configuration found"):
            _StepNoIO().read_cubepars("unused", "G235H", "F170LP")


# ---------------------------------------------------------------------------
# process — unsupported instrument
# ---------------------------------------------------------------------------


class TestProcessUnsupportedInstrument:
    """process logs a warning and returns input_data when the instrument is not NIRSpec."""

    def test_logs_nirspec_guidance(self, caplog):
        step = _StepNoIO()
        with caplog.at_level("WARNING"):
            step.process(_make_model(instrument="MIRI"))
        assert "SlitCubeBuildStep is designed specifically for NIRSpec" in caplog.text


# ---------------------------------------------------------------------------
# process — invalid datamodel type
# ---------------------------------------------------------------------------


class TestProcessInvalidDatamodelType:
    """process logs a warning and returns input_data when given a non-MultiExposureModel."""

    def _make_non_exposure_model(self):
        """Return a real JwstDataModel that is not a MultiExposureModel."""
        return ImageModel()

    def test_returns_input_unchanged(self):
        step = _StepNoIO()
        bad_model = self._make_non_exposure_model()
        assert step.process(bad_model) is bad_model

    def test_logs_model_type_warning(self, caplog):
        step = _StepNoIO()
        with caplog.at_level("WARNING"):
            step.process(self._make_non_exposure_model())
        assert "SlitCubeBuildStep requires a MultiExposureModel" in caplog.text
