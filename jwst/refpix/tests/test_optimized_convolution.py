import numpy as np
from stdatamodels.jwst.datamodels import RampModel, SIRSKernelModel
from jwst.refpix.optimized_convolution import (
    make_kernels,
    get_conv_kernel_coeffs,
    apply_conv_kernel,
)


# create the ConvKernelModel
ckm = {
    "nrcb1": {
        "gamma": np.array(
            [
                [
                    0.8737859 + 0.0j,
                    0.72877103 - 0.01848215j,
                    0.7474708 + 0.00441926j,
                    0.7596158 - 0.01682704j,
                    0.7710808 - 0.00618939j,
                ],
                [
                    0.37835783 + 0.0j,
                    0.27234325 - 0.03058944j,
                    0.38302818 + 0.03056235j,
                    0.36819065 - 0.02578794j,
                    0.3908449 + 0.02115744j,
                ],
                [
                    0.36443716 + 0.0j,
                    0.335223 + 0.02436169j,
                    0.32699308 - 0.02325623j,
                    0.3830375 - 0.01340938j,
                    0.39612782 + 0.00736016j,
                ],
                [
                    0.00335188 + 0.0j,
                    0.01759672 - 0.01073076j,
                    0.04302938 + 0.00353758j,
                    0.08149841 - 0.00643084j,
                    0.07274915 - 0.002046j,
                ],
            ]
        ),
        "zeta": np.array(
            [
                [
                    0.14007446 + 0.0000000e00j,
                    0.2371146 + 1.6455967e-02j,
                    0.22727127 - 5.9413449e-03j,
                    0.2090475 + 7.0676603e-03j,
                    0.20298977 + 2.2992526e-05j,
                ],
                [
                    0.6206608 + 0.0j,
                    0.680701 + 0.02468053j,
                    0.57776874 - 0.03374288j,
                    0.5873975 + 0.01647749j,
                    0.5693782 - 0.02531039j,
                ],
                [
                    0.6543285 + 0.0j,
                    0.6167225 - 0.02665404j,
                    0.6405862 + 0.01494319j,
                    0.57719606 + 0.00970044j,
                    0.57160926 - 0.01088286j,
                ],
                [
                    1.0137521 + 0.0j,
                    0.9492664 + 0.0071805j,
                    0.92866725 - 0.00784425j,
                    0.8868761 - 0.00237024j,
                    0.89918566 - 0.00323711j,
                ],
            ]
        ),
    }
}
conv_kernel_model = SIRSKernelModel(ckm)


def mk_data_mdl(data, instrument, detector):
    # create input_model
    input_model = RampModel(data=data)
    input_model.meta.instrument.name = instrument
    input_model.meta.instrument.detector = detector
    input_model.meta.subarray.name = "FULL"
    return input_model


def test_get_conv_kernel_coeffs():
    detector = "NRCB1"
    gamma, zeta = get_conv_kernel_coeffs(conv_kernel_model, detector)
    assert gamma is not None
    assert zeta is not None


def test_mk_kernels():
    detector = "nothing"
    gaussmooth = 1
    halfwidth = 30
    kernels = make_kernels(conv_kernel_model, detector, gaussmooth, halfwidth)
    assert kernels is None


def test_apply_conv_kernel():
    data = np.zeros((3, 3, 2048, 2048)) + 1.999
    instrument, detector = "NIRCAM", "NRCB1"
    input_model = mk_data_mdl(data, instrument, detector)
    gaussmooth = 1
    halfwidth = 30
    kernels = make_kernels(conv_kernel_model, detector, gaussmooth, halfwidth)
    sigreject = 4
    result = apply_conv_kernel(input_model.data[1, 1, ...], kernels, sigreject=sigreject)
    compare = np.ones((1, 1, 2048, 2048))
    assert compare.all() == result.all()
