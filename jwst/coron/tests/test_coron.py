import numpy as np
import numpy.testing as npt

from jwst import datamodels

from jwst.coron import imageregistration
from jwst.coron import klip


def test_fourier_imshift():
    """ Test of fourier_imshift() in imageregistration.py """

    image = np.zeros((5, 5), dtype=np.float32)
    image[1:4, 1:4] += 1.0
    image[2, 2] += 2.0
    shift = [1.2, 0.6]
    result = imageregistration.fourier_imshift(image, shift)

    truth = np.array(
        [
            [-0.11534368, 0.09685412, -0.31489129, 0.07308346, -0.203199],
            [0.17121927, -0.14163561, 0.5043889, 0.06434016, 0.39715385],
            [0.11266705, -0.06854077, 0.75829127, 2.03638809, 1.36342888],
            [0.04314152, -0.00521786, 0.65394254, 2.4800922, 1.46183018],
            [0.22979281, -0.18803098, 0.7125186, 0.25274376, 0.62498253],
        ]
    )
    npt.assert_allclose(result, truth, atol=1e-6)


def test_shift_subtract():
    """ Test of shift_subtract() in imageregistration.py """

    target = np.arange((15), dtype=np.float32).reshape((3, 5))
    reference = target + 0.1
    reference[1, 0] -= 0.2
    reference[2, 0] += 2.3
    mask = target * 0 + 1
    mask[1, 1] = 0
    mask[1, 2] = 0
    params = (0.2, -0.3, 1)

    result = imageregistration.shift_subtract(params, reference, target, mask)

    truth = np.array(
        [
            0.65830479,
            1.38243832,
            0.46556956,
            1.24234379,
            0.17542511,
            -4.85323936,
            -0.0,
            -0.0,
            -3.06869101,
            -3.70808286,
            -0.03657461,
            2.61462573,
            2.73619353,
            2.81983023,
            2.56214928,
        ]
    )

    npt.assert_allclose(result, truth, atol=1e-6)


def test_align_fourierLSQ():
    """ Test of align_fourierLSQ() in imageregistration.py """

    target = np.arange((15), dtype=np.float32).reshape((3, 5))
    reference = target + 0.1
    reference[1, 0] -= 0.2
    reference[2, 0] += 2.3
    mask = target * 0 + 1
    mask[1, 1] = 0
    mask[1, 2] = 0

    shifts = imageregistration.align_fourierLSQ(reference, target, mask)
    truth = np.array([-0.0899215, -0.01831958, 0.96733475])

    npt.assert_allclose(shifts, truth, atol=1e-6)


def test_align_array():
    """ Test of align_array() in imageregistration.py """

    temp = np.arange((15), dtype=np.float32).reshape((3, 5))
    targ = np.zeros((3, 3, 5))
    targ[:] = temp
    targ[0, 1, 1] += 0.3
    targ[0, 2, 1] += 0.7
    targ[0, 0, 3] -= 0.5
    targ[0, 1, 2] -= 1.3
    targ[1, 0, 1] += 0.7
    targ[1, 2, 1] += 0.2
    targ[1, 2, 3] += 0.8
    targ[1, 2, 2] -= 1.8
    targ[2, 1, 1] += 0.9
    targ[2, 2, 2] -= 0.5
    targ[2, 1, 0] += 0.8
    targ[2, 1, 2] += 0.8

    ref = temp.copy()
    ref[0, 0] += 0.5
    ref[1, 2] -= 0.4
    ref[2, 2] -= 1.4
    ref[0, 2] += 1.3
    ref[1, 4] -= 0.6

    mask = temp * 0 + 1
    mask[1, 1] = 0
    mask[1, 2] = 0
    aligned, shifts = imageregistration.align_array(ref, targ, mask)

    truth_aligned = np.array(
        [
            [
                [0.12432412, 1.08545608, 2.13530046, 2.58060399, 4.11698591],
                [4.81329092, 6.07738984, 5.50928942, 7.7584455, 8.82976601],
                [10.10986446, 11.78015814, 12.07965768, 13.0802872, 14.11918026],
            ],
            [
                [0.27298558, 1.55238983, 2.16767105, 2.83462221, 4.31946497],
                [5.06350107, 5.58666974, 6.89508812, 7.57902246, 9.00320211],
                [10.17131384, 11.2266214, 10.30098165, 13.46529354, 14.46117244],
            ],
            [
                [0.15066264, 1.18760843, 2.14701273, 3.22159633, 4.17955834],
                [5.4048206, 6.54242247, 7.4322509, 7.63745739, 8.60847585],
                [10.19237215, 11.23021861, 11.70579103, 13.20361126, 14.15614128],
            ],
        ]
    )

    truth_shifts = np.array(
        [
            [-0.00773077, -0.01635325, 1.02677575],
            [-0.08423515, -0.01487977, 1.0143966],
            [0.00768275, -0.03127504, 1.02081983],
        ]
    )

    npt.assert_allclose(aligned, truth_aligned, atol=1e-6)
    npt.assert_allclose(shifts, truth_shifts, atol=1e-6)


def test_align_models():
    """ Test of align_models() in imageregistration.py """

    temp = np.arange((15), dtype=np.float32).reshape((3, 5))
    targ = np.zeros((3, 3, 5))
    targ[:] = temp
    targ[0, 1, 1] += 0.3
    targ[0, 2, 1] += 0.7
    targ[0, 0, 3] -= 0.5
    targ[0, 1, 2] -= 1.3
    targ[1, 0, 1] += 0.7
    targ[1, 2, 1] += 0.2
    targ[1, 2, 3] += 0.8
    targ[1, 2, 2] -= 1.8
    targ[2, 1, 1] += 0.9
    targ[2, 2, 2] -= 0.5
    targ[2, 1, 0] += 0.8
    targ[2, 1, 2] += 0.8

    ref = targ.copy()
    ref[1, 2, 3] -= 5.0
    ref[2, 0, 3] -= 1.6

    mask = ref[0, :, :]
    mask = ref[0, :, :] * 0 + 1
    mask[1, 1] = 0
    mask[1, 2] = 0

    targ_mod = datamodels.CubeModel(data=targ)
    mask_mod = datamodels.ImageModel(data=mask)
    ref_mod = datamodels.CubeModel(data=ref)

    am_results = imageregistration.align_models(ref_mod, targ_mod, mask_mod)
    results_sub = am_results.data[:3, :2, 2, :3]

    truth_results_sub = np.array(
        [
            [[10.0, 11.7, 12.0], [10.036278, 11.138131, 10.180669]],
            [[10.053974, 11.1953335, 11.993213], [10.36224, 10.805556, 10.274276]],
            [[9.988604, 11.33026, 11.968155], [10.024722, 10.971058, 10.108071]],
        ]
    )

    npt.assert_allclose(results_sub, truth_results_sub, atol=1e-6)


def test_KLT():
    """ Test of KarhunenLoeveTransform() in klip.py """

    temp = np.arange((15), dtype=np.float32).reshape((3, 5))
    refs = np.zeros((3, 3, 5))
    refs[:] = temp
    refs[0, 1, 1] += 0.3
    refs[0, 2, 1] += 0.7
    refs[0, 0, 3] -= 0.5
    refs[0, 1, 2] -= 1.3
    refs[1, 0, 1] += 0.7
    refs[1, 2, 1] += 0.2
    refs[1, 2, 3] += 0.8
    refs[1, 2, 2] -= 1.8
    refs[2, 1, 1] += 0.9
    refs[2, 2, 2] -= 0.5
    refs[2, 1, 0] += 0.8
    refs[2, 1, 2] += 0.8

    rshape = refs.shape
    nrefs = rshape[0]
    refs = refs.reshape(nrefs, rshape[1] * rshape[2])

    # Make each ref image have zero mean
    for k in range(nrefs):
        refs[k] -= np.mean(refs[k], dtype=np.float64)

    klvect, eigval, eigvect = klip.KarhunenLoeveTransform(refs, normalize=True)

    truth_klvect = np.array(
        [
            [
                0.42288619,
                0.34884115,
                0.30246786,
                0.25254087,
                0.18204953,
                0.10593695,
                0.03757055,
                0.01225229,
                -0.05878713,
                -0.1189963,
                -0.17920546,
                -0.25776279,
                -0.2541062,
                -0.37564538,
                -0.42004212,
            ],
            [
                -0.16042986,
                0.07486249,
                -0.12636643,
                0.11188599,
                -0.09230299,
                0.04286651,
                -0.05806697,
                0.65210377,
                -0.02417611,
                -0.00714439,
                0.00988733,
                -0.22042976,
                -0.59112697,
                0.3104232,
                0.0780142,
            ],
            [
                -0.08789825,
                -0.3672306,
                -0.08482407,
                -0.02768044,
                -0.08174988,
                0.33098098,
                0.35055338,
                0.47863215,
                -0.07560151,
                -0.07406442,
                -0.07252733,
                -0.22908779,
                0.39578645,
                -0.3889097,
                -0.06637896,
            ],
        ]
    )

    npt.assert_allclose(klvect, truth_klvect, atol=1e-6)

    truth_eigval = np.array([59.09269614, 0.22691493, 0.16324608])
    npt.assert_allclose(eigval, truth_eigval, atol=1e-6)

    truth_eigvect = np.array(
        [
            [-0.59148838, -0.78859008, -0.16812845],
            [-0.56851174, 0.55574161, -0.60658525],
            [-0.57178309, 0.2632051, 0.77703743],
        ]
    )
    npt.assert_allclose(eigvect, truth_eigvect, atol=1e-6)


def test_klip():
    """ Test of klip() in klip.py """

    target_model_data = np.array(
        [
            [
                [0.84456396, 0.75297576, 0.79752606, 1.2884853],
                [0.98718196, 0.83042985, 0.76917756, 0.9760479],
                [1.4539411, 1.1395805, 0.9635729, 0.92798233],
            ],
            [
                [0.84456396, 0.75297576, 0.79752606, 1.2884853],
                [0.98718196, 0.83042985, 0.76917756, 0.9760479],
                [1.4539411, 1.1395805, 0.9635729, 0.92798233],
            ],
        ],
        dtype=np.float32,
    )

    target_model = datamodels.CubeModel(data=target_model_data)

    refs_model_data = np.array(
        [
            [
                [
                    [0.8174741, 0.74938107, 0.73527235, 1.3193785],
                    [1.0032778, 0.8247719, 0.78944355, 0.99227476],
                    [1.4609907, 1.1605016, 0.9564753, 0.9186427],
                ],
                [
                    [0.86789674, 0.7998908, 0.8557136, 1.2926395],
                    [0.97756547, 0.7788742, 0.75892323, 0.9819151],
                    [1.495664, 1.1455023, 1.002115, 0.92159164],
                ],
                [
                    [0.84426856, 0.7719569, 0.8088021, 1.2781427],
                    [0.98734635, 0.8125992, 0.77424014, 0.9934157],
                    [1.419994, 1.1546139, 0.961317, 0.95088667],
                ],
            ],
            [
                [
                    [0.8174741, 0.74938107, 0.73527235, 1.3193785],
                    [1.0032778, 0.8247719, 0.78944355, 0.99227476],
                    [1.4609907, 1.1605016, 0.9564753, 0.9186427],
                ],
                [
                    [0.86789674, 0.7998908, 0.8557136, 1.2926395],
                    [0.97756547, 0.7788742, 0.75892323, 0.9819151],
                    [1.495664, 1.1455023, 1.002115, 0.92159164],
                ],
                [
                    [0.84426856, 0.7719569, 0.8088021, 1.2781427],
                    [0.98734635, 0.8125992, 0.77424014, 0.9934157],
                    [1.419994, 1.1546139, 0.961317, 0.95088667],
                ],
            ],
        ],
        dtype=np.float32,
    )

    refs_model = datamodels.QuadModel(data=refs_model_data)

    # Call the KLIP routine
    truncate = 50
    psf_sub, psf_fit = klip.klip(target_model, refs_model, truncate)

    truth_psf_sub_data = np.array(
        [
            [
                [0.00677787, -0.01639403, 0.00753349, -0.00190761],
                [-0.00016282, 0.02299659, -0.00647708, -0.00992276],
                [0.00798155, -0.00867007, -0.00287747, 0.00112234],
            ],
            [
                [0.00677787, -0.01639403, 0.00753349, -0.00190761],
                [-0.00016282, 0.02299659, -0.00647708, -0.00992276],
                [0.00798155, -0.00867007, -0.00287747, 0.00112234],
            ],
        ]
    )

    # psf_fit is currently not used in the code, co not compared here
    npt.assert_allclose(psf_sub.data, truth_psf_sub_data, atol=1e-6)
