import pytest
from stdatamodels.jwst.datamodels import CubeModel

from jwst.assign_wcs.nircam import _get_tsgrism_source_center


def test_tsgrism_siaf_reference_position_is_converted_to_zero_indexing():
    model = CubeModel((1, 2, 2))
    model.meta.wcsinfo.siaf_xref_sci = 1024.5
    model.meta.wcsinfo.siaf_yref_sci = 35.0

    assert _get_tsgrism_source_center(model) == (1023.5, 34.0)
    model.close()


def test_tsgrism_source_center_requires_siaf_reference_position():
    model = CubeModel((1, 2, 2))
    model.meta.wcsinfo.siaf_xref_sci = None
    model.meta.wcsinfo.siaf_yref_sci = 35.0

    with pytest.raises(ValueError, match="XREF_SCI is missing"):
        _get_tsgrism_source_center(model)
    model.close()
