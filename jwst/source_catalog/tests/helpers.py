import numpy as np
import stdatamodels.jwst.datamodels as dm
from photutils.datasets import make_gwcs


def make_nircam_model():
    """
    Create a NIRCam ImageModel with synthetic data for testing.

    Returns
    -------
    model : dm.ImageModel
        NIRCam image model with data, error, and weight arrays.
    """
    rng = np.random.default_rng(seed=123)
    data = rng.normal(0, 0.5, size=(101, 101))
    data[20:80, 10:20] = 1.4
    data[20:30, 20:45] = 1.4
    data[20:80, 55:65] = 7.2
    data[70:80, 65:87] = 7.2
    data[45:55, 65:87] = 7.2
    data[20:30, 65:87] = 7.2
    data[55:75, 82:92] = 7.2
    data[25:45, 82:92] = 7.2

    wht = np.ones(data.shape)
    wht[0:10, :] = 0.0
    err = np.abs(data) / 10.0
    model = dm.ImageModel(data, wht=wht, err=err)
    model.meta.bunit_data = "MJy/sr"
    model.meta.bunit_err = "MJy/sr"
    model.meta.photometry.pixelarea_steradians = 1.0e-13
    model.meta.wcs = make_gwcs(data.shape)
    model.meta.wcsinfo = {
        "ctype1": "RA---TAN",
        "ctype2": "DEC--TAN",
        "dec_ref": 11.99875540218638,
        "ra_ref": 22.02351763251896,
        "roll_ref": 0.005076934167039675,
        "v2_ref": 86.039011,
        "v3_ref": -493.385704,
        "v3yangle": -0.07385127,
        "vparity": -1,
        "wcsaxes": 2,
        "crpix1": 50,
        "crpix2": 50,
    }
    model.meta.instrument = {
        "channel": "LONG",
        "detector": "NRCALONG",
        "filter": "F444W",
        "lamp_mode": "NONE",
        "module": "A",
        "name": "NIRCAM",
        "pupil": "CLEAR",
    }
    model.meta.exposure.type = "NRC_IMAGE"
    model.meta.observation.date = "2021-01-01"
    model.meta.observation.time = "00:00:00"

    return model


def make_nircam_model_without_apcorr():
    """
    Create a NIRCam ImageModel without aperture correction.

    Returns
    -------
    model : dm.ImageModel
        NIRCam image model with data, error, and weight arrays.
    """
    rng = np.random.default_rng(seed=123)
    data = rng.normal(0, 0.5, size=(101, 101))
    data[20:80, 10:20] = 1.4
    data[20:30, 20:45] = 1.4
    data[20:80, 55:65] = 7.2
    data[70:80, 65:87] = 7.2
    data[45:55, 65:87] = 7.2
    data[20:30, 65:87] = 7.2
    data[55:75, 82:92] = 7.2
    data[25:45, 82:92] = 7.2

    wht = np.ones(data.shape)
    wht[0:10, :] = 0.0
    err = np.abs(data) / 10.0
    model = dm.ImageModel(data, wht=wht, err=err)
    model.meta.bunit_data = "MJy/sr"
    model.meta.bunit_err = "MJy/sr"
    model.meta.photometry.pixelarea_steradians = 1.0e-13
    model.meta.wcs = make_gwcs(data.shape)
    model.meta.wcsinfo = {
        "ctype1": "RA---TAN",
        "ctype2": "DEC--TAN",
        "dec_ref": 11.99875540218638,
        "ra_ref": 22.02351763251896,
        "roll_ref": 0.005076934167039675,
        "v2_ref": 86.039011,
        "v3_ref": -493.385704,
        "v3yangle": -0.07385127,
        "vparity": -1,
        "wcsaxes": 2,
        "crpix1": 50,
        "crpix2": 50,
    }
    model.meta.instrument = {
        "channel": "LONG",
        "detector": "NRCALONG",
        "filter": "F2550WR",
        "lamp_mode": "NONE",
        "module": "A",
        "name": "NIRCAM",
        "pupil": "CLEAR",
    }
    model.meta.exposure.type = "NRC_IMAGE"
    model.meta.observation.date = "2021-01-01"
    model.meta.observation.time = "00:00:00"

    return model


def mock_apcorr():
    """
    Mock a NIRCam image APCORR reference file.

    The table contains F444W/CLEAR entries only, borrowed from
    jwst_nircam_apcorr_0005.fits.

    Returns
    -------
    apcorr_model : `~stdatamodels.jwst.datamodels.NrcImgApcorrModel`
        The APCORR reference file datamodel.
    """
    table = np.array(
        [
            ("F444W", "CLEAR", 0.22483717, 1.0, 4.4558034, 8.0, 13.0),
            ("F444W", "CLEAR", 0.3, 1.2335707, 3.3402927, 8.0, 13.0),
            ("F444W", "CLEAR", 0.4, 1.4828625, 2.5056577, 8.0, 13.0),
            ("F444W", "CLEAR", 0.5, 1.7693356, 2.0051568, 8.0, 13.0),
            ("F444W", "CLEAR", 0.55829775, 2.0, 1.7964456, 8.0, 13.0),
            ("F444W", "CLEAR", 0.6, 2.2247522, 1.6723331, 8.0, 13.0),
            ("F444W", "CLEAR", 0.69345284, 3.0, 1.4497877, 8.0, 13.0),
            ("F444W", "CLEAR", 0.7, 3.073966, 1.4365366, 8.0, 13.0),
            ("F444W", "CLEAR", 0.75860673, 4.0, 1.3297256, 8.0, 13.0),
            ("F444W", "CLEAR", 0.8, 4.6166277, 1.2638301, 8.0, 13.0),
            ("F444W", "CLEAR", 0.8189596, 5.0, 1.2365663, 8.0, 13.0),
            ("F444W", "CLEAR", 0.8460141, 6.0, 1.2030405, 8.0, 13.0),
            ("F444W", "CLEAR", 0.8556371, 7.0, 1.1968731, 8.0, 13.0),
            ("F444W", "CLEAR", 0.8650267, 8.0, 1.1922662, 8.0, 13.0),
            ("F444W", "CLEAR", 0.8779693, 9.0, 1.183854, 8.0, 13.0),
            ("F444W", "CLEAR", 0.89006805, 10.0, 1.1778655, 8.0, 13.0),
            ("F444W", "CLEAR", 0.8982966, 11.0, 1.1784167, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9037039, 12.0, 1.1840535, 8.0, 13.0),
            ("F444W", "CLEAR", 0.90804386, 13.0, 1.1924243, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9120697, 14.0, 1.2025543, 8.0, 13.0),
            ("F444W", "CLEAR", 0.91617095, 15.0, 1.2139562, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9203339, 16.0, 1.2267189, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9240836, 17.0, 1.2416534, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9275178, 18.0, 1.258755, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9308524, 19.0, 1.2778363, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9339979, 20.0, 1.2992089, 8.0, 13.0),
            ("F444W", "CLEAR", 0.93682176, 21.0, 1.3233073, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9393199, 22.0, 1.3504053, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9415356, 23.0, 1.3807377, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9435172, 24.0, 1.4145732, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9453089, 25.0, 1.4522386, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9469987, 26.0, 1.4940228, 8.0, 13.0),
            ("F444W", "CLEAR", 0.94865125, 27.0, 1.5403177, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9502525, 28.0, 1.5917816, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9518031, 29.0, 1.6491717, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9532612, 30.0, 1.7135346, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9545909, 31.0, 1.7861507, 8.0, 13.0),
            ("F444W", "CLEAR", 0.95583403, 32.0, 1.8683583, 8.0, 13.0),
            ("F444W", "CLEAR", 0.95701283, 33.0, 1.961902, 8.0, 13.0),
            ("F444W", "CLEAR", 0.95813316, 34.0, 2.0690663, 8.0, 13.0),
            ("F444W", "CLEAR", 0.95917517, 35.0, 2.192933, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9601479, 36.0, 2.3374317, 8.0, 13.0),
            ("F444W", "CLEAR", 0.96109635, 37.0, 2.5076241, 8.0, 13.0),
            ("F444W", "CLEAR", 0.96206695, 38.0, 2.7104056, 8.0, 13.0),
            ("F444W", "CLEAR", 0.96305627, 39.0, 2.9558666, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9640153, 40.0, 3.2592366, 8.0, 13.0),
            ("F444W", "CLEAR", 0.96491206, 41.0, 3.6436963, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9657537, 42.0, 4.1460023, 8.0, 13.0),
            ("F444W", "CLEAR", 0.96656126, 43.0, 4.8288126, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9673487, 44.0, 5.809133, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9681209, 45.0, 7.333614, 8.0, 13.0),
            ("F444W", "CLEAR", 0.968873, 46.0, 10.026812, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9695971, 47.0, 16.062393, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9702898, 48.0, 41.79031, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9709452, 49.0, -65.53901, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9715617, 50.0, -18.081371, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9721526, 51.0, -10.395007, 8.0, 13.0),
            ("F444W", "CLEAR", 0.972731, 52.0, -7.2501287, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9732986, 53.0, -5.5404916, 8.0, 13.0),
            ("F444W", "CLEAR", 0.97384495, 54.0, -4.4664207, 8.0, 13.0),
            ("F444W", "CLEAR", 0.97436297, 55.0, -3.729307, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9748628, 56.0, -3.1924465, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9753585, 57.0, -2.7842891, 8.0, 13.0),
            ("F444W", "CLEAR", 0.975852, 58.0, -2.463658, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9763338, 59.0, -2.205188, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9767934, 60.0, -1.9924473, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9772291, 61.0, -1.8143551, 8.0, 13.0),
            ("F444W", "CLEAR", 0.97764635, 62.0, -1.6631613, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9780517, 63.0, -1.5332658, 8.0, 13.0),
            ("F444W", "CLEAR", 0.97844714, 64.0, -1.4205109, 8.0, 13.0),
            ("F444W", "CLEAR", 0.978831, 65.0, -1.3217468, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9792027, 66.0, -1.2345515, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9795644, 67.0, -1.1570348, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9799164, 68.0, -1.0876929, 8.0, 13.0),
            ("F444W", "CLEAR", 0.9802554, 69.0, -1.0253146, 8.0, 13.0),
            ("F444W", "CLEAR", 0.98058164, 70.0, -0.968919, 8.0, 13.0),
        ],
        dtype=(
            [
                ("filter", "S12"),
                ("pupil", "S15"),
                ("eefraction", "<f4"),
                ("radius", "<f4"),
                ("apcorr", "<f4"),
                ("skyin", "<f4"),
                ("skyout", "<f4"),
            ]
        ),
    )
    apcorr_model = dm.NrcImgApcorrModel()
    apcorr_model.apcorr_table = table
    return apcorr_model
