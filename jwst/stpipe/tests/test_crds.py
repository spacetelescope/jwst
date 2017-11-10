from __future__ import absolute_import, print_function

from os.path import join, dirname, basename
import shutil
import tempfile

from nose.tools import raises

from astropy.io import fits
import six

from jwst.stpipe import Step
import crds

TMP_DIR = None
TMP_FITS = None


def setup():
    global TMP_DIR, TMP_FITS
    TMP_DIR = tempfile.mkdtemp()
    TMP_FITS = join(TMP_DIR, 'tmp.fits')


def teardown():
    shutil.rmtree(TMP_DIR)

class CrdsStep(Step):
    reference_file_types = ['flat']

    def process(self, input_file):
        from ... import datamodels

        with datamodels.open(input_file) as dm:
            self.ref_filename = self.get_reference_file(dm, 'flat')
        return datamodels.DataModel()


def test_crds_step():
    """Nominal working CrdsStep flat fetch for valid dataset."""
    _run_flat_fetch_on_dataset('data/crds.fits')


@raises(Exception)
def test_crds_step_bad():
    """Run CrdsStep on a dataset with invalid FILTER."""
    _run_flat_fetch_on_dataset('data/crds_bad.fits')


def _run_flat_fetch_on_dataset(dataset_path):
    from ... import datamodels
    step = CrdsStep()
    with datamodels.ImageModel(join(dirname(__file__), dataset_path)) as input_file:
        step.run(input_file)
    print(step.ref_filename)
    assert basename(step.ref_filename) == "jwst_nircam_flat_0105.fits"

def test_crds_step_override():
    """Run CRDS step with override parameter bypassing CRDS lookup."""
    from ... import datamodels

    step = CrdsStep(override_flat=join(dirname(__file__), 'data/flat.fits'))
    with datamodels.ImageModel(join(dirname(__file__), 'data/crds.fits')) as input_file:
        result = step.run(input_file)
    assert step.ref_filename.endswith('data/flat.fits')
    assert result.meta.ref_file.flat.name.endswith('flat.fits')

    result.to_fits(TMP_FITS)

    with fits.open(TMP_FITS) as hdulist:
        assert hdulist[0].header['R_FLAT'].endswith('flat.fits')

@raises(crds.CrdsError)
def test_crds_failed_getreferences_parameter():
    """Run crds.getreferences() with invalid FILTER."""
    header = {
        '_extra_fits.PRIMARY.IRAF-TLM': '2013-12-12T15:56:30',
        'meta.date': '2014-07-22T15:53:19.893683',
        'meta.filename': 'crds.fits',
        'meta.instrument.detector': 'yyyyNRCA1yyyy', # whack this parameter
        'meta.instrument.filter': 'F140M',
        'meta.instrument.name': 'NIRCAM',
        'meta.instrument.pupil': 'CLEAR',
        'meta.observation.date': '2012-04-22',
        'meta.origin': 'NOAO-IRAF FITS Image Kernel July 2003',
        'meta.subarray.name': 'FULL',
        'meta.subarray.xsize': 2048,
        'meta.subarray.xstart': 1,
        'meta.subarray.ysize': 2048,
        'meta.subarray.ystart': 1,
        'meta.telescope': 'JWST'
        }
    crds.getreferences(header, reftypes=["flat"])

@raises(crds.CrdsError)
def test_crds_failed_getreferences_reftype():
    """Run crds.getreferences() with an invalid reftypes list."""
    header = {
        '_extra_fits.PRIMARY.IRAF-TLM': '2013-12-12T15:56:30',
        'meta.date': '2014-07-22T15:53:19.893683',
        'meta.filename': 'crds.fits',
        'meta.instrument.detector': 'NRCA1',
        'meta.instrument.filter': 'F140M',
        'meta.instrument.name': 'NIRCAM',
        'meta.instrument.pupil': 'CLEAR',
        'meta.observation.date': '2012-04-22',
        'meta.origin': 'NOAO-IRAF FITS Image Kernel July 2003',
        'meta.subarray.name': 'FULL',
        'meta.subarray.xsize': 2048,
        'meta.subarray.xstart': 1,
        'meta.subarray.ysize': 2048,
        'meta.subarray.ystart': 1,
        'meta.telescope': 'JWST'
        }
    crds.getreferences(header, reftypes=["foo"])

# def test_crds_failed_getreferences_bad_context():
#     import crds
#     header = {
#         u'_extra_fits.PRIMARY.IRAF-TLM': '2013-12-12T15:56:30',
#         u'meta.date': '2014-07-22T15:53:19.893683',
#         u'meta.filename': 'crds.fits',
#         u'meta.instrument.detector': 'NRCA1',
#         u'meta.instrument.filter': 'F140M',
#         u'meta.instrument.name': 'NIRCAM',
#         u'meta.instrument.pupil': 'CLEAR',
#         u'meta.observation.date': '2012-04-22',
#         u'meta.origin': 'NOAO-IRAF FITS Image Kernel July 2003',
#         u'meta.subarray.name': 'FULL',
#         u'meta.subarray.xsize': 2048,
#         u'meta.subarray.xstart': 1,
#         u'meta.subarray.ysize': 2048,
#         u'meta.subarray.ystart': 1,
#         u'meta.telescope': 'JWST'
#         }
#     assert_raises(crds.getreferences, header, reftypes=["foo"], context="jwst_9942.pmap")


