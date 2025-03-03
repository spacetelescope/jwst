"""Test the srctype step on exposures with various settings."""

from stdatamodels.jwst import datamodels

from jwst.srctype import srctype
import pytest


def test_background_target_set():

    # An exposure flagged as background target
    step_input = datamodels.ImageModel((10, 10))

    step_input.meta.exposure.type = "NRS_IFU"
    step_input.meta.observation.bkgdtarg = True
    step_input.meta.dither.primary_type = "2-POINT-NOD"
    step_input.meta.target.source_type = "POINT"

    output = srctype.set_source_type(step_input)

    # Result should be EXTENDED regardless of other step_input settings
    assert output.meta.target.source_type == "EXTENDED"


def test_background_target_unset():

    # An exposure without BKGDTARG set at all
    step_input = datamodels.ImageModel((10, 10))

    step_input.meta.exposure.type = "NRS_IFU"
    step_input.meta.dither.primary_type = "2-POINT-NOD"
    step_input.meta.target.source_type = "EXTENDED"

    output = srctype.set_source_type(step_input)

    # If BKGDTARG is missing, next test should be based on the
    # value of PATTTYPE, which in this case should return POINT.
    assert output.meta.target.source_type == "EXTENDED"


def test_nrsifu_nodded():

    # An exposure using a NIRSpec IFU NOD dither pattern
    step_input = datamodels.ImageModel((10, 10))

    step_input.meta.exposure.type = "NRS_IFU"
    step_input.meta.observation.bkgdtarg = False
    step_input.meta.dither.primary_type = "2-POINT-NOD"
    step_input.meta.target.source_type = "EXTENDED"

    output = srctype.set_source_type(step_input)

    # Result should be POINT regardless of step_input setting
    assert output.meta.target.source_type == "EXTENDED"


def test_mirmrs_nodded():

    # An exposure using a MIRI MRS NOD dither pattern
    step_input = datamodels.ImageModel((10, 10))

    step_input.meta.exposure.type = "MIR_MRS"
    step_input.meta.observation.bkgdtarg = False
    step_input.meta.dither.primary_type = "2-POINT"
    step_input.meta.dither.optimized_for = "POINT-SOURCE"
    step_input.meta.target.source_type = "EXTENDED"

    output = srctype.set_source_type(step_input)

    # Result should be POINT regardless of step_input setting
    assert output.meta.target.source_type == "EXTENDED"


def test_user_input():

    # An exposure with the value set upstream by the user
    step_input = datamodels.ImageModel((10, 10))

    step_input.meta.exposure.type = "NRS_IFU"
    step_input.meta.observation.bkgdtarg = False
    step_input.meta.dither.primary_type = "4-POINT"
    step_input.meta.target.source_type = "POINT"

    output = srctype.set_source_type(step_input)

    # Result should be POINT regardless of other step_input settings
    assert output.meta.target.source_type == "POINT"


def test_mrs_unknown():

    # An exposure with upstream step_input UNKNOWN
    step_input = datamodels.ImageModel((10, 10))

    step_input.meta.exposure.type = "MIR_MRS"
    step_input.meta.observation.bkgdtarg = False
    step_input.meta.dither.primary_type = "CYCLING"
    step_input.meta.target.source_type = "UNKNOWN"

    output = srctype.set_source_type(step_input)

    # Result should be EXTENDED regardless of other step_input settings
    assert output.meta.target.source_type == "EXTENDED"


def test_ifu_unknown():

    # An exposure with upstream step_input UNKNOWN
    step_input = datamodels.ImageModel((10, 10))

    step_input.meta.exposure.type = "NRS_IFU"
    step_input.meta.observation.bkgdtarg = False
    step_input.meta.dither.primary_type = "4-POINT"
    step_input.meta.target.source_type = "UNKNOWN"

    output = srctype.set_source_type(step_input)

    # Result should be EXTENDED regardless of other step_input settings
    assert output.meta.target.source_type == "EXTENDED"


def test_exptype():

    # An exposure with an unrecognized EXP_TYPE
    step_input = datamodels.ImageModel((10, 10))

    step_input.meta.exposure.type = "NRS_LAMP"
    step_input.meta.observation.bkgdtarg = False
    step_input.meta.dither.primary_type = "NONE"

    output = srctype.set_source_type(step_input)

    # Result should be UNKNOWN regardless of other step_input settings
    assert output.meta.target.source_type == "UNKNOWN"


def test_no_sourcetype():

    # An exposure without the SRCTYPE keyword present at all
    step_input = datamodels.ImageModel((10, 10))

    step_input.meta.exposure.type = "MIR_LRS-FIXEDSLIT"
    step_input.meta.observation.bkgdtarg = False
    step_input.meta.dither.primary_type = "2-POINT"

    output = srctype.set_source_type(step_input)

    # Result should be POINT regardless of other step_input settings
    assert output.meta.target.source_type == "POINT"


def test_nrs_msaspec():
    """Test for when exposure type is NRS_MSASPEC."""

    step_input = datamodels.MultiSlitModel()
    step_input.meta.exposure.type = "NRS_MSASPEC"

    slits = [{"source_id": 1, "stellarity": 0.9},
             {"source_id": 2, "stellarity": -1},
             {"source_id": 3, "stellarity": 0.5}]

    for slit in slits:
        step_input.slits.append(slit)

    result = srctype.set_source_type(step_input)

    assert result.slits[0].source_type == "POINT"
    assert result.slits[1].source_type == "POINT"
    assert result.slits[2].source_type == "EXTENDED"


def test_nrs_fixedslit():
    """Test for when exposure type is NRS_FIXEDSLIT."""

    step_input = datamodels.MultiSlitModel()
    step_input.meta.exposure.type = "NRS_FIXEDSLIT"
    step_input.meta.instrument.fixed_slit = "S200A1"
    step_input.meta.target.source_type = "POINT"

    slits = [{"name": "S200A2"},
             {"name": "S200A1"},
             {"name": "S1600A1"}]

    for slit in slits:
        step_input.slits.append(slit)

    result = srctype.set_source_type(step_input)

    assert result.slits[0].source_type == "EXTENDED"
    assert result.slits[1].source_type == "POINT"
    assert result.slits[2].source_type == "EXTENDED"

    # Check that source_xpos and source_ypos are reset to zero for primary slit
    # if the source is EXTENDED
    step_input.slits[1].source_xpos = 0.5
    step_input.slits[1].source_ypos = -0.5
    result = srctype.set_source_type(step_input, "EXTENDED")
    assert result.slits[1].source_xpos == 0.0
    assert result.slits[1].source_ypos == 0.0


def test_valid_user_spec():
    """Overwrite source_type with valid user-specified value."""

    step_input = datamodels.ImageModel((10, 10))

    step_input.meta.exposure.type = "NRS_IFU"
    step_input.meta.target.source_type = "POINT"

    # User sets to valid value "EXTENDED"
    output = srctype.set_source_type(step_input, "EXTENDED")

    # Result should be match user override value
    assert output.meta.target.source_type == "EXTENDED"


@pytest.mark.parametrize("exptype", ["NRS_BRIGHTOBJ", "NRC_TSGRISM",
                                     "NIS_SOSS", "MIR_LRS-SLITLESS"])
def test_tso_types(exptype):
    """Test for when visit is tso."""

    step_input = datamodels.ImageModel()
    step_input.meta.observation.bkgdtarg = False
    step_input.meta.dither.primary_type = "NONE"
    step_input.meta.visit.tsovisit = True
    step_input.meta.exposure.type = exptype

    result = srctype.set_source_type(step_input)

    assert result.meta.target.source_type == "POINT"


def test_exptype_is_none():
    """Test for when exposure type is None."""

    with pytest.raises(RuntimeError):
        step_input = datamodels.ImageModel((10, 10))
        step_input.meta.exposure.type = None
        srctype.set_source_type(step_input)
