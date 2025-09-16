import numpy as np

from jwst.wfss_contam.sens1d import create_1d_sens, get_photom_data


def test_get_photom_data(photom_ref_model):
    filter_name = "GR150C"
    pupil = "F200W"
    order = 1
    refmodel = photom_ref_model.copy()  # we are going to modify this in-place

    # scramble the wavelengths to ensure they are sorted
    n = len(photom_ref_model.phot_table["wavelength"])
    waves = photom_ref_model.phot_table["wavelength"][0].copy()
    waves_inv = waves[::-1]
    refmodel.phot_table["wavelength"] = [waves_inv] * n

    # set the 0th index of relresponse to a different value to ensure it's also sorted
    # according to the wavelength order
    relresp = refmodel.phot_table["relresponse"][0].copy()
    relresp[0] = 99
    refmodel.phot_table["relresponse"] = [relresp] * n

    ref_waves, ref_relresp = get_photom_data(refmodel, filter_name, pupil, order)

    np.testing.assert_array_equal(ref_waves, waves)
    assert ref_relresp.shape == waves.shape

    # get the flux scaling out of the table: F200W is the 3rd row
    scalar_conversion = refmodel.phot_table["photmjsr"][2]
    assert np.isclose(ref_relresp[0], 10 * scalar_conversion)
    assert np.isclose(ref_relresp[-1], 99 * scalar_conversion)


def test_create_1d_sens(photom_ref_model):
    filter_name = "GR150C"
    pupil = "F200W"
    order = 1
    ref_waves, ref_relresp = get_photom_data(photom_ref_model, filter_name, pupil, order)

    # create a set of wavelengths to interpolate the response onto
    # make it wider than the reference wavelengths to test zero-filling and no_cal mask
    data_waves = np.linspace(1.5, 2.5, 50)

    sens, no_cal = create_1d_sens(data_waves, ref_waves, ref_relresp)

    assert sens.shape == data_waves.shape
    assert no_cal.shape == data_waves.shape

    assert np.all(sens[:10] == 0)  # check zero-filling at the start
    assert np.all(no_cal[:10])  # check no_cal is True at the start
    assert np.all(sens[-10:] == 0)  # check zero-filling at the end
    assert np.all(no_cal[-10:])  # check no_cal is True at the end

    scalar_conversion = photom_ref_model.phot_table["photmjsr"][2]
    assert np.allclose(sens[10:-10], scalar_conversion)
    assert np.all(~no_cal[10:-10])  # check no_cal is False
