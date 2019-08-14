"""
Test the srctype step on exposures with various settings
"""
from jwst import datamodels
from jwst.srctype import srctype
import pytest

def test_background_target_set():

    # An exposure flagged as background target
    input = datamodels.ImageModel((10,10))

    input.meta.exposure.type = 'NRS_IFU'
    input.meta.observation.bkgdtarg = True
    input.meta.dither.primary_type = '2-POINT-NOD'
    input.meta.target.source_type = 'POINT'

    output = srctype.set_source_type(input)

    # Result should be EXTENDED regardless of other input settings
    assert output.meta.target.source_type == 'EXTENDED'

def test_background_target_unset():

    # An exposure without BKGDTARG set at all
    input = datamodels.ImageModel((10,10))

    input.meta.exposure.type = 'NRS_IFU'
    input.meta.dither.primary_type = '2-POINT-NOD'
    input.meta.target.source_type = 'EXTENDED'

    output = srctype.set_source_type(input)

    # If BKGDTARG is missing, next test should be based on the
    # value of PATTTYPE, which in this case should return POINT.
    assert output.meta.target.source_type == 'POINT'

def test_nrsifu_nodded():

    # An exposure using a NIRSpec IFU NOD dither pattern
    input = datamodels.ImageModel((10,10))

    input.meta.exposure.type = 'NRS_IFU'
    input.meta.observation.bkgdtarg = False
    input.meta.dither.primary_type = '2-POINT-NOD'
    input.meta.target.source_type = 'EXTENDED'

    output = srctype.set_source_type(input)

    # Result should be POINT regardless of input setting
    assert output.meta.target.source_type == 'POINT'

def test_mirmrs_nodded():

    # An exposure using a MIRI MRS NOD dither pattern
    input = datamodels.ImageModel((10,10))

    input.meta.exposure.type = 'MIR_MRS'
    input.meta.observation.bkgdtarg = False
    input.meta.dither.primary_type = 'POINT-SOURCE'
    input.meta.target.source_type = 'EXTENDED'

    output = srctype.set_source_type(input)

    # Result should be POINT regardless of input setting
    assert output.meta.target.source_type == 'POINT'

def test_user_input():

    # An exposure with the value set upstream by the user
    input = datamodels.ImageModel((10,10))

    input.meta.exposure.type = 'NRS_IFU'
    input.meta.observation.bkgdtarg = False
    input.meta.dither.primary_type = '4-POINT'
    input.meta.target.source_type = 'POINT'

    output = srctype.set_source_type(input)

    # Result should be POINT regardless of other input settings
    assert output.meta.target.source_type == 'POINT'

def test_unknown():

    # An exposure with upstream input UNKNOWN
    input = datamodels.ImageModel((10,10))

    input.meta.exposure.type = 'MIR_MRS'
    input.meta.observation.bkgdtarg = False
    input.meta.dither.primary_type = 'CYCLING'
    input.meta.target.source_type = 'UNKNOWN'

    output = srctype.set_source_type(input)

    # Result should be EXTENDED regardless of other input settings
    assert output.meta.target.source_type == 'EXTENDED'

def test_exptype():

    # An exposure with an unrecognized EXP_TYPE
    input = datamodels.ImageModel((10,10))

    input.meta.exposure.type = 'NRS_LAMP'
    input.meta.observation.bkgdtarg = False
    input.meta.dither.primary_type = 'NONE'

    output = srctype.set_source_type(input)

    # Result should be EXTENDED regardless of other input settings
    assert output.meta.target.source_type == 'EXTENDED'

def test_no_sourcetype():

    # An exposure without the SRCTYPE keyword present at all
    input = datamodels.ImageModel((10,10))

    input.meta.exposure.type = 'MIR_LRS-FIXEDSLIT'
    input.meta.observation.bkgdtarg = False
    input.meta.dither.primary_type = '2-POINT'

    output = srctype.set_source_type(input)

    # Result should be POINT regardless of other input settings
    assert output.meta.target.source_type == 'POINT'

def test_tso_types():

    # Exposure without the SRCTYPE keyword present at all,
    # but it's a TSO mode
    input = datamodels.ImageModel((10,10))
    input.meta.observation.bkgdtarg = False
    input.meta.dither.primary_type = 'NONE'

    # All results should be POINT
    input.meta.exposure.type = 'NRS_BRIGHTOBJ'
    output = srctype.set_source_type(input)
    assert output.meta.target.source_type == 'POINT'

    del input.meta.target.source_type
    input.meta.exposure.type = 'NRC_TSGRISM'
    output = srctype.set_source_type(input)
    assert output.meta.target.source_type == 'POINT'

    del input.meta.target.source_type
    input.meta.exposure.type = 'NIS_SOSS'
    input.meta.visit.tsovisit = True
    output = srctype.set_source_type(input)
    assert output.meta.target.source_type == 'POINT'

    del input.meta.target.source_type
    input.meta.exposure.type = 'MIR_LRS-SLITLESS'
    input.meta.visit.tsovisit = True
    output = srctype.set_source_type(input)
    assert output.meta.target.source_type == 'POINT'

def test_nrs_msaspec():
    """Test for when exposure type is NRS_MSASPEC
    """
    input = datamodels.MultiSlitModel()
    input.meta.exposure.type = "NRS_MSASPEC"

    slits = [{'source_id':1, 'stellarity':0.9},
             {'source_id':2, 'stellarity':-1},
             {'source_id':3, 'stellarity':0.5}]

    for slit in slits:
        input.slits.append(slit)

    result = srctype.set_source_type(input)

    assert(result.slits[0].source_type == 'POINT')
    assert(result.slits[1].source_type == 'UNKNOWN')
    assert(result.slits[2].source_type == 'EXTENDED')

def test_is_tso():
    """ Test for when visit is tso.
    """
    input = datamodels.ImageModel((10,10))
    input.meta.exposure.type = "FGS_DARK"
    input.meta.visit.tsovisit = True
    result = srctype.set_source_type(input)

    assert(result.meta.target.source_type == "POINT")

@pytest.mark.xfail
def test_exptype_is_none():
    """ Test for when exposure type is None.
    """
    input = datamodels.ImageModel((10,10))
    input.meta.exposure.type = None
    srctype.set_source_type(input)