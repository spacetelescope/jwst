import os
import shutil
from astropy.io import fits
from jwst.extract_2d.extract_2d_step import Extract2dStep

BIGDATA = os.environ['TEST_BIGDATA']

def test_extract2d_nrs_msa():
    """

    Regression test of extract_2d step performed on NIRSpec fixed slit data.

    """
    input_files = [BIGDATA+'/nirspec/test_extract_2d/msa/jw00038001001_01101_00001_NRS1_assign_wcs.fits',
                   BIGDATA+'/nirspec/test_extract_2d/msa/msa_configuration.fits'
                   ]

    try:
        os.remove("extract2d2_output.fits")
        for f in input_files:
            os.remove(os.path.basename(f))
    except:
        pass

    try:
        # This is needed because we can't run the step with the
        # MSA configuration file not in the current dir.
        # When this is fixed, the copy can be removed.
        for f in input_file:
            shutil.copyfile(f, os.path.join(".", os.path.basename(f)))
    except:
        raise OSError("Could not copy inputs files")

    Extract2dStep.call('jw00038001001_01101_00001_NRS1_assign_wcs.fits',
                       config_file='extract_2d.cfg',
                       output_file='extract2d2_output.fits'
    )
    h = fits.open('extract2d2_output.fits')
    href = fits.open(BIGDATA+'/nirspec/test_extract_2d/msa/jw00038001001_01101_00001_NRS1_assign_wcs_extract_2d.fits')
    '''
    newh = fits.HDUList([h['primary'],h['sci',1],h['err',1],h['dq',1],
                         h['sci',2],h['err',2],h['dq',2],
                         h['sci',3],h['err',3],h['dq',3],
                         h['sci',4],h['err',4],h['dq',4],
                         h['sci',5],h['err',5],h['dq',5]])
    newhref = fits.HDUList([href['primary'],href['sci',1],href['err',1],href['dq',1],
                            href['sci',2],href['err',2],href['dq',2],
                            href['sci',3],href['err',3],href['dq',3],
                            href['sci',4],href['err',4],href['dq',4],
                            href['sci',5],href['err',5],href['dq',5]])
    
    result = fits.diff.FITSDiff(newh, newhref,
                                ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                                rtol = 0.00001
    )
    '''
    result = fits.diff.FITSDiff(h, href,
                                ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                                rtol = 0.00001)

    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

