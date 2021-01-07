"""
Unit tests for pathloss correction
"""

from jwst.datamodels import MultiSlitModel, PathlossModel
from jwst.pathloss.pathloss import (calculate_pathloss_vector,
                                    get_aperture_from_model,
                                    get_center,
                                    interpolate_onto_grid,
                                    is_pointsource)
from jwst.pathloss.pathloss import do_correction
import numpy as np
import pytest


def test_get_center_ifu():
    """get_center assumes IFU targets are centered @ (0.0, 0.0)"""

    x_pos,y_pos = get_center("NRS_IFU", None)

    assert(x_pos==y_pos==0.0)


def test_get_center_attr_err():
    """if no center provided for modes that are not IFU,
    center is assigned to (0.0, 0.0)"""

    datmod = MultiSlitModel()
    x_pos, y_pos = get_center("NRS_MSASPEC", datmod)

    assert(x_pos==y_pos==0.0)


def test_get_center_exp_type():
    """if exp_type is not in NRS, center is returned (0.0,0.0)"""
    datmod = MultiSlitModel()
    x_pos, y_pos = get_center("NRC_IMAGE", datmod)

    assert(x_pos==y_pos==0.0)


def test_get_center_exptype():
    """ If exptype is "NRS_MSASPEC" | "NRS_FIXEDSLIT" | "NRS_BRIGHTOBJ" and
    source_xpos and source_ypos exist in datamod.slits, make sure it's returned"""

    datmod = MultiSlitModel()
    datmod.slits.append({'source_xpos':1, 'source_ypos':2})

    for exptype in ["NRS_MSASPEC", "NRS_FIXEDSLIT", "NRS_BRIGHTOBJ"]:
        x_pos, y_pos = get_center(exptype, datmod.slits[0])

        assert(x_pos==1)
        assert(y_pos==2)


# Begin get_aperture_from_model tests
def test_get_app_from_model_null():
    """If exp_type isn't the NRS or NIS specific mode,
    routine returns None"""

    datmod = MultiSlitModel()
    datmod.meta.exposure.type = 'NRC_IMAGE'

    result = get_aperture_from_model(datmod, None)

    assert(result is None)


def test_get_aper_from_model_fixedslit():
    """For a given exposures aperture, make sure the correct
    aperture reference data is returned for fixedslit mode"""

    datmod = PathlossModel()
    datmod.apertures.append({'name':'S200A1'})
    datmod.meta.exposure.type = 'NRS_FIXEDSLIT'

    result = get_aperture_from_model(datmod, 'S200A1')

    assert(result == datmod.apertures[0])


def test_get_aper_from_model_msa():
    """For a given exposures aperture, make sure the correct
    aperture reference data is returned for MSA mode"""

    datmod = PathlossModel()
    datmod.apertures.append({'shutters':5})
    datmod.meta.exposure.type = 'NRS_MSASPEC'

    result = get_aperture_from_model(datmod, 5)

    assert(result == datmod.apertures[0])


# Begin calculate_pathloss_vector tests.
def test_calculate_pathloss_vector_pointsource_data():
    """Calculate pathloss vector for 3D pathloss data"""

    datmod = PathlossModel()

    ref_data = {'pointsource_data':np.ones((10,10,10), dtype=np.float32),
                'pointsource_wcs': {'crval2': -0.5, 'crpix2': 1.0, 'cdelt2': 0.05,
                                    'cdelt3': 1, 'crval1': -0.5, 'crpix1': 1.0,
                                    'crpix3': 1.0, 'crval3': 1, 'cdelt1': 0.05}}

    datmod.apertures.append(ref_data)

    wavelength, pathloss, is_inside_slitlet = calculate_pathloss_vector(datmod.apertures[0].pointsource_data,
                                                                        datmod.apertures[0].pointsource_wcs,
                                                                        0.0, 0.0)

    # Wavelength array is calculated with this: crval3 +(float(i+1) - crpix3)*cdelt3
    # Where i is the iteration of np.arange(wavesize) which is the 1st dimension of the pointsource
    # data array.
    wavelength_comparison = np.array([1 + (float(i+1) - 1.0)*1 for i in np.arange(10)])
    assert(np.allclose(wavelength,wavelength_comparison))

    # pathloss vector gets assigned at beginning of calculate_pathloss_vector and in this
    # case, doesnt change (np.zeros(wavesize, dtype=np.float32))
    pathloss_comparison = np.zeros(10, dtype=np.float32)
    assert(np.all(pathloss==pathloss_comparison))

    # With the current wcs values, the logic should be returning False
    assert(is_inside_slitlet==False)

def test_calculate_pathloss_vector_uniform_data():
    """Calculate the pathloss vector for uniform data arrays."""

    datmod = PathlossModel()

    ref_data = {'uniform_data':np.ones((10,), dtype=np.float32),
                'uniform_wcs': {'crpix1': 1.0, 'cdelt1': 1, 'crval1': 1}}

    datmod.apertures.append(ref_data)

    wavelength, pathloss, _ = calculate_pathloss_vector(datmod.apertures[0].uniform_data,
                                                        datmod.apertures[0].uniform_wcs,
                                                        0.0, 0.0)

    # Wavelength array is calculated with this: crval1 +(float(i+1) - crpix1)*cdelt1
    # Where i is the iteration of np.arange(wavesize) which is the shape of the uniform
    # data array.
    comparison = np.array([1 +(float(i+1) - 1)*1 for i in np.arange(10)])
    assert(np.all(wavelength==comparison))

    # The same array is returned in this case
    assert(np.all(datmod.apertures[0].uniform_data == pathloss))


def test_calculate_pathloss_vector_interpolation():
    """Calculate the pathloss vector for when interpolation is necessary."""

    datmod = PathlossModel()

    ref_data = {'pointsource_data':np.ones((10,10,10), dtype=np.float32),
                'pointsource_wcs': {'crval2': -0.5, 'crpix2': 1.0, 'cdelt2': 0.5,
                                    'cdelt3': 1.0, 'crval1': -0.5, 'crpix1': 1.0,
                                    'crpix3': 1.0, 'crval3': 1.0, 'cdelt1': 0.5}}

    datmod.apertures.append(ref_data)

    wavelength, pathloss, is_inside_slitlet = calculate_pathloss_vector(datmod.apertures[0].pointsource_data,
                                                                        datmod.apertures[0].pointsource_wcs,
                                                                        0.0, 0.0)

    # Wavelength array is calculated with this: crval3 +(float(i+1) - crpix3)*cdelt3
    # Where i is the iteration of np.arange(wavesize) which is the 1st dimension of the pointsource
    # data array.
    wavelength_comparison = np.array([1 + (float(i+1) - 1.0)*1 for i in np.arange(10)])
    assert(np.all(wavelength==wavelength_comparison))

    # In this instance we interpolate to get the array for pathloss VS wavelength.
    # With the current values inside of the of the pointsource_wcs starting at line 143 of pathloss.py
    # dx1 = 1 - int(1) = 0.0
    # dx2 = 1 - dx1 = 1.0
    # dy1 = 1 - int(1) = 0.0
    # dy2 = 1 - dx1 = 1.0
    # This means a11 == a12 == a21 == 0 and a22 == 1
    # This simplifies that pathloss vector to:
    # pathloss_vector = (a22*pathloss_ref[:,i,j]) = (1*pathloss_ref[:1,j])
    # Thus pathloss == the input array to the function.
    pathloss_comparison = datmod.apertures[0].pointsource_data
    assert(np.all(pathloss==pathloss_comparison))

    # With the current wcs values, the logic should be returning True
    assert(is_inside_slitlet==True)


def test_is_pointsource():
    """Check to see if object it point source"""

    point_source = None
    result = is_pointsource(point_source)
    assert(result == False)

    point_source = 'point'
    result = is_pointsource(point_source)
    assert(result == True)

    point_source = 'not a point'
    result = is_pointsource(point_source)
    assert(result == False)

def test_do_correction_msa_slit_size_eq_0():
    """If slits have size 0, quit calibration."""

    datmod = MultiSlitModel()
    datmod.slits.append({'data':np.array([])})
    pathlossmod = PathlossModel()
    datmod.meta.exposure.type = 'NRS_MSASPEC'

    result, _ = do_correction(datmod, pathlossmod)
    assert(result.meta.cal_step.pathloss == 'COMPLETE')


def test_do_correction_fixed_slit_exception():
    """If no matching aperture name found, exit."""

    datmod = MultiSlitModel()
    # Give input_model aperture name
    datmod.slits.append({'data':np.array([]), 'name':'S200A1'})
    # Do assign pathloss model aperture with similar name.
    pathlossmod = PathlossModel()
    datmod.meta.exposure.type = 'NRS_FIXEDSLIT'

    result, _ = do_correction(datmod, pathlossmod)
    assert(result.meta.cal_step.pathloss == 'COMPLETE')


def test_do_correction_nis_soss_tso():
    """If observation is tso, skip correction"""

    datmod = MultiSlitModel()
    pathlossmod = PathlossModel()
    datmod.meta.exposure.type = 'NIS_SOSS'
    datmod.meta.visit.tsovisit = True

    result, _ = do_correction(datmod, pathlossmod)
    assert(result.meta.cal_step.pathloss == 'SKIPPED')


def test_do_correction_nis_soss_pupil_position_is_none():
    """If pupil_position is None, skip correction"""

    datmod = MultiSlitModel()
    pathlossmod = PathlossModel()
    datmod.meta.exposure.type = 'NIS_SOSS'
    datmod.meta.visit.tsovisit = False
    datmod.meta.instrument.pupil_position = None

    result, _ = do_correction(datmod, pathlossmod)
    assert(result.meta.cal_step.pathloss == 'SKIPPED')

def test_do_correction_nis_soss_aperture_is_none():
    """If no matching aperture is found, skip correction."""

    datmod = MultiSlitModel()
    # Is FULL an option for NIRISS?
    # The test doesn't care but something to remember.
    datmod.slits.append({'data':np.array([]), 'name':'FULL'})
    # Don't assign pathloss model aperture with similar name
    pathlossmod = PathlossModel()
    datmod.meta.exposure.type = 'NIS_SOSS'
    datmod.meta.visit.tsovisit = False
    datmod.meta.instrument.pupil_position = 1

    result, _ = do_correction(datmod, pathlossmod)
    assert(result.meta.cal_step.pathloss == 'SKIPPED')


@pytest.mark.skip(reason="Fraction calculation in interpolate_onto_grid needs refactoring.")
def test_interpolate_onto_grid():
    # Mock wavelength vector, grid and pathloss vector.
    wavelength_grid = np.arange(1, 101).reshape(10,10) * 1.1
    wavelength_vector = np.arange(1, 11, dtype='float64')
    pathloss_vector = np.arange(1, 11, dtype='float64')

    # Call interpolate onto grid
    result = interpolate_onto_grid(wavelength_grid,
                                   wavelength_vector,
                                   pathloss_vector)

    # Before interpolation is done in interpolate_onto_grid, the vectors are padded
    # so interpolation that happens outside of the grid are NaN.
    extended_pathloss_vector = np.zeros(len(pathloss_vector) + 2)
    extended_pathloss_vector[1:-1] = pathloss_vector
    extended_pathloss_vector[0] = np.nan
    extended_pathloss_vector[-1] = np.nan
    extended_wavelength_vector = np.zeros(len(wavelength_vector) + 2)
    extended_wavelength_vector[1:-1] = wavelength_vector
    extended_wavelength_vector[0] = wavelength_vector[0] - 0.1
    extended_wavelength_vector[-1] = wavelength_vector[-1] + 0.1

    # Call numpy interpolation to get truth.
    result_comparison = np.interp(wavelength_grid, extended_wavelength_vector, extended_pathloss_vector)

    assert np.testing.assert_array_equal(result, result_comparison)
