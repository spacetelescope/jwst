import pytest

from jwst.assign_wcs.assign_wcs_step import AssignWcsStep
from jwst.clean_flicker_noise import tso_median_image as tmi
from jwst.clean_flicker_noise.tests import helpers


def test_nrs_bots_whitelight_flux_count_mismatch(monkeypatch):
    input_model = helpers.make_nrs_bots_rateints()
    wcs_model = AssignWcsStep.call(input_model)

    original_run = tmi.WhiteLightStep.run

    def drop_last_flux(self, *args, **kwargs):
        table = original_run(self, *args, **kwargs)
        return table[:-1]

    monkeypatch.setattr(tmi.WhiteLightStep, "run", drop_last_flux)

    with pytest.raises(ValueError, match="does not match number of integrations"):
        tmi.make_median_image(wcs_model, wcs_model)

    input_model.close()
    wcs_model.close()
