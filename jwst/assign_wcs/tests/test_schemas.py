import inspect
import sys

from astropy.modeling import models
from astropy import units as u
import pytest

from jwst.datamodels import DistortionModel, ReferenceFileModel
from jwst.datamodels import wcs_ref_models
from jwst.datamodels.wcs_ref_models import _SimpleModel


def find_all_wcs_ref_models_classes():
    clsmembers = inspect.getmembers(sys.modules[wcs_ref_models.__name__], inspect.isclass)
    classes = [cls for name,cls in clsmembers if issubclass(cls, ReferenceFileModel)]
    classes.remove(_SimpleModel)
    return classes


@pytest.fixture
def distortion_model():
    """Create a distortion model that should pass all validation"""
    m = models.Shift(1) & models.Shift(2)
    dist = DistortionModel(model=m, input_units=u.pixel, output_units=u.arcsec)

    dist.meta.reftype = "distortion"
    dist.meta.instrument.name = "NIRCAM"
    dist.meta.instrument.detector = "NRCA1"
    dist.meta.instrument.p_pupil = "F162M|F164N|CLEAR|"
    dist.meta.instrument.pupil = "F162M"
    dist.meta.exposure.p_exptype = "NRC_IMAGE|NRC_TSIMAGE|NRC_FLAT|NRC_LED|NRC_WFSC|"
    dist.meta.exposure.type = "NRC_IMAGE"
    dist.meta.psubarray = "FULL|SUB64P|SUB160)|SUB160P|SUB320|SUB400P|SUB640|"
    dist.meta.subarray.name = "FULL"

    # Populate the following so that no validation warnings or errors happen
    dist.meta.instrument.module = "A"
    dist.meta.instrument.channel = "SHORT"
    dist.meta.input_units = u.degree
    dist.meta.output_units = u.degree
    dist.meta.description = "NIRCam distortion reference file"
    dist.meta.author = "Hank the Septopus"
    dist.meta.pedigree = "Cleveland"
    dist.meta.useafter = "2000-01-01T00:00:00"

    return dist


def test_distortion_schema(distortion_model, tmpdir):
    """Make sure DistortionModel roundtrips"""
    path = str(tmpdir.join("test_dist.asdf"))
    dist = distortion_model
    dist.save(path)

    with pytest.warns(None) as report:
        with DistortionModel(path) as dist1:
            assert dist1.meta.instrument.p_pupil == dist.meta.instrument.p_pupil
            assert dist1.meta.instrument.pupil == dist.meta.instrument.pupil
            assert dist1.meta.exposure.p_exptype == dist.meta.exposure.p_exptype
            assert dist1.meta.exposure.type == dist.meta.exposure.type
            assert dist1.meta.psubarray == dist.meta.psubarray
            assert dist1.meta.subarray.name == dist.meta.subarray.name
        assert len(report) == 0


def test_distortion_strict_validation(distortion_model):
    """Make sure strict validation works"""
    distortion_model.validate()


def test_distortion_schema_bad_valueerror(distortion_model):
    """Check that ValueError is raised for ReferenceFile missing items"""
    dist = DistortionModel(distortion_model, strict_validation=True)
    dist.meta.author = None

    with pytest.raises(ValueError):
        dist.validate()


def test_distortion_schema_bad_assertionerror(distortion_model):
    """Check that AssertionError is raised for distortion-specific missing items"""
    dist = DistortionModel(distortion_model, strict_validation=True)
    dist.meta.instrument.channel = None

    with pytest.raises(AssertionError):
        dist.validate()


@pytest.mark.parametrize("cls", find_all_wcs_ref_models_classes())
def test_simplemodel_subclasses(cls):
    """Test that expected validation errors are raised"""
    model = cls()
    with pytest.warns(None) as report:
        model.validate()
    assert len(report) >= 1

    model = cls(strict_validation=True)
    with pytest.raises((ValueError, KeyError)):
        model.validate()
