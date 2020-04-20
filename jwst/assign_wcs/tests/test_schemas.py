from astropy.modeling import models
from astropy import units as u
from jwst.datamodels import DistortionModel


def test_distortion_schema(tmpdir):
    """Make sure DistortionModel roundtrips"""
    m = models.Shift(1) & models.Shift(2)
    dist = DistortionModel(model=m, input_units=u.pixel, output_units=u.arcsec)

    dist.meta.instrument.name = "NIRCAM"
    dist.meta.instrument.detector = "NRCA1"
    dist.meta.instrument.p_pupil = "F162M|F164N|CLEAR|"
    dist.meta.instrument.pupil = "F162M"
    dist.meta.exposure.p_exptype = "NRC_IMAGE|NRC_TSIMAGE|NRC_FLAT|NRC_LED|NRC_WFSC|"
    dist.meta.exposure.type = "NRC_IMAGE"
    dist.meta.psubarray = "FULL|SUB64P|SUB160)|SUB160P|SUB320|SUB400P|SUB640|"
    dist.meta.subarray.name = "FULL"
    path = str(tmpdir.join("test_dist.asdf"))
    dist.save(path)

    with DistortionModel(path) as dist1:
        assert dist1.meta.instrument.p_pupil == dist.meta.instrument.p_pupil
        assert dist1.meta.instrument.pupil == dist.meta.instrument.pupil
        assert dist1.meta.exposure.p_exptype == dist.meta.exposure.p_exptype
        assert dist1.meta.exposure.type == dist.meta.exposure.type
        assert dist1.meta.psubarray == dist.meta.psubarray
        assert dist1.meta.subarray.name == dist.meta.subarray.name
