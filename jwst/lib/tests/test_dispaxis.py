"""
Test for dispaxis
"""

from jwst.lib import dispaxis


def test_dispaxis_1():
    value = dispaxis.get_dispersion_direction("FGS_IMAGE", "junk", "junk",
                                              "junk")
    assert value is None


def test_dispaxis_2():
    value = dispaxis.get_dispersion_direction("MIR_4QPM", "junk", "junk",
                                              "junk")
    assert value is None


def test_dispaxis_3():
    value = dispaxis.get_dispersion_direction("MIR_LRS-SLITLESS", "junk",
                                              "junk", "junk")
    assert value == 2


def test_dispaxis_4():
    value = dispaxis.get_dispersion_direction("MIR_LRS-FIXEDSLIT", "junk",
                                              "junk", "junk")
    assert value == 2


def test_dispaxis_5():
    value = dispaxis.get_dispersion_direction("NIS_SOSS", "junk", "junk",
                                              "junk")
    assert value == 1


def test_dispaxis_6():
    value = dispaxis.get_dispersion_direction("NRS_FIXEDSLIT", "junk", "junk",
                                              "junk")
    assert value == 1


def test_dispaxis_7():
    value = dispaxis.get_dispersion_direction("NRS_IFU", "junk", "junk",
                                              "junk")
    assert value == 1


def test_dispaxis_8():
    value = dispaxis.get_dispersion_direction("NRS_MSASPEC", "junk", "junk",
                                              "junk")
    assert value == 1


def test_dispaxis_9():
    value = dispaxis.get_dispersion_direction("NIS_WFSS", "junk", "GR150R",
                                              "junk")
    assert value == 2


def test_dispaxis_10():
    value = dispaxis.get_dispersion_direction("NIS_WFSS", "junk", "GR150C",
                                              "junk")
    assert value == 1


def test_dispaxis_11():
    value = dispaxis.get_dispersion_direction("NRC_GRISM", "junk", "junk",
                                              "GRISMR")
    assert value == 1


def test_dispaxis_12():
    value = dispaxis.get_dispersion_direction("NRC_GRISM", "junk", "junk",
                                              "GRISMC")
    assert value == 2


def test_dispaxis_13():
    value = dispaxis.get_dispersion_direction("NRC_TSGRISM", "junk", "junk",
                                              "GRISMR")
    assert value == 1


def test_dispaxis_14():
    value = dispaxis.get_dispersion_direction("NRC_WFSS", "junk", "junk",
                                              "GRISMR")
    assert value == 1


def test_dispaxis_15():
    value = dispaxis.get_dispersion_direction("NRC_WFSS", "junk", "junk",
                                              "GRISMC")
    assert value == 2


def test_dispaxis_16():
    value = dispaxis.get_dispersion_direction("nrc_grism", "junk", "junk",
                                              "missing")
    assert value is None
