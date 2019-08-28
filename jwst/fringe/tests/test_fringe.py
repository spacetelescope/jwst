"""

Unit tests for fringe correction

"""

import pytest
import numpy as np
import numpy.random as rn

from jwst.datamodels import IFUImageModel
from jwst.datamodels import FringeModel
from jwst.fringe import fringe

FRINGE_CONSTANT = 2. # correction will be input data divided by this factor

def test_data_correction( setup_inputs ):
    ''' Test both good and NaN pixels. '''

    shape = (4,5)
    input_model, fringe_model = setup_inputs( shape )

    # Make 1 bad pixel
    input_model.data[0,0] = np.nan
    input_model.err[0,0] = np.nan

    #  Do the correction()
    output_model = fringe.do_correction( input_model, fringe_model )

    # Check that correction was done on pixels with valid values for both
    #    SCI and ERR arrays
    good_pix = np.where( np.isfinite(input_model.data) )

    assert (output_model.data[good_pix] ==
        (input_model.data * FRINGE_CONSTANT)[good_pix]).all()
    assert (output_model.err[good_pix] ==
        (input_model.err * FRINGE_CONSTANT)[good_pix]).all()

    # Check that correction was not done on pixel with NaN values for both SCI
    #     and ERR arrays (i.e. these pixels have not been corrected)
    assert( np.isnan(output_model.data[0,0] ))
    assert( np.isnan(output_model.err[0,0] ))


@pytest.fixture
def setup_inputs():
    ''' Create input and fringe models.'''
    def _setup( shape=(2,2)):

        input_data = (np.ones( shape[0]*shape[1])).reshape(shape) * 6.
        input_err = rn.random_sample( shape )
        input_model = IFUImageModel( data=input_data, err=input_err )

        fringe_data = (np.ones( shape[0]*shape[1])).reshape(shape)/FRINGE_CONSTANT
        fringe_model = FringeModel( data=fringe_data )

        return input_model, fringe_model

    return _setup
