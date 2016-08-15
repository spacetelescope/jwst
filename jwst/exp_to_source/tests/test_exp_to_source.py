from jwst.datamodels import MultiSlitModel

from ..exp_to_source import exp_to_source

INPUT_FILE = '/Users/eisenham/Documents/ssbdev/testdata/jwst_data/NIRSpec/MSA/NRSMOS-MODEL-23_NRS1_SloperPipeline_assign_wcs_extract_2d.fits'


def test_exp_to_source():
    inputs = [MultiSlitModel(INPUT_FILE),
              MultiSlitModel(INPUT_FILE)]
    outputs = exp_to_source(inputs)
    assert len(outputs) == 15
    assert len(outputs[outputs.keys()[0]].slits) == 2
    assert (inputs[0].slits[0].data == outputs['S1001'].slits[0].data).all()
    assert (inputs[1].slits[14].data == outputs['S1015'].slits[1].data).all()
