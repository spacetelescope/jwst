import numpy as np
import pytest
from stdatamodels.jwst.datamodels import ImageModel

from jwst.resample.resample import input_jwst_model_to_dict


def test_ivm_rejects_mismatched_var_rnoise():
    model = ImageModel((4, 5))
    model.var_rnoise = np.ones((3, 5), dtype=np.float32)

    with pytest.raises(ValueError, match="IVM weighting requires VAR_RNOISE"):
        input_jwst_model_to_dict(
            model, weight_type="ivm", enable_var=False, compute_err=None
        )

    model.close()
