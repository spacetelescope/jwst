"""

Unit tests for EMI correction

"""


from jwst.emicorr import EmiCorrStep, emicorr
from stdatamodels.jwst.datamodels import Level1bModel, EmiModel


def test_get_subarcase():
    # set up a real subarray
    subarray, readpatt = 'MASK1550', 'FASTR1'
    subname_real = emicorr.get_subarcase(subarray, readpatt)

    # set up a fake configuration
    subarray, fake_readpatt = 'FULL', 'kkkkk'
    subname_fake = emicorr.get_subarcase(subarray, readpatt)

    # test if we get the right configuration
    compare_real = ['MASK1550', 82, 23968]
    compare_fake = ['FULL_kkkkk', None, None]

    assert compare_real == subname_real
    assert compare_fake == subname_fake


def test_mk_reffile_waveform():
    im = Level1bModel()
    emicorr_ref_filename = ''
    mdl = emicorr.mk_reffile_waveform(input_model, emicorr_ref_filename,
                                    save_mdl=False)
