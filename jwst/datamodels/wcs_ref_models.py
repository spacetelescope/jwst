import numpy as np
from astropy.modeling.core import Model
from astropy import units as u
from . import model_base

from .extension import BaseExtension
from jwst.transforms.jwextension import JWSTExtension
from gwcs.extension import GWCSExtension
from .reference import ReferenceFileModel

__all__ = ['DistortionModel', 'DistortionMRSModel', 'SpecwcsModel', 'RegionsModel',
           'WavelengthrangeModel', 'CameraModel', 'CollimatorModel', 'OTEModel',
           'FOREModel', "FPAModel", 'IFUPostModel', 'IFUFOREModel', 'IFUSlicerModel',
           'MSAModel', 'FilteroffsetModel', 'DisperserModel',
           'NIRCAMGrismModel', 'NIRISSGrismModel', 'WaveCorrModel']


class _SimpleModel(ReferenceFileModel):
    """
    A model for a reference file of type "distortion".
    """
    schema_url = None
    reftype = None

    def __init__(self, init=None, model=None, input_units=None, output_units=None, **kwargs):

        super(_SimpleModel, self).__init__(init=init, **kwargs)
        if model is not None:
            self.model = model
        if input_units is not None:
            self.meta.input_units = input_units
        if output_units is not None:
            self.meta.output_units = output_units
        if init is None:
            try:
                self.populate_meta()
            except NotImplementedError:
                pass

    def on_save(self, path=None):
        self.meta.reftype = self.reftype

    def populate_meta(self):
        """
        Subclasses can overwrite this to populate specific meta keywords.
        """
        raise NotImplementedError

    def to_fits(self):
        raise NotImplementedError("FITS format is not supported for this file.")

    def validate(self):
        super(_SimpleModel, self).validate()
        assert isinstance(self.model, Model)
        assert self.meta.instrument.name in ["NIRCAM", "NIRSPEC", "MIRI", "TFI", "FGS", "NIRISS"]


class DistortionModel(_SimpleModel):
    """
    A model for a reference file of type "distortion".
    """
    schema_url = "distortion.schema.yaml"
    reftype = "distortion"

    def validate(self):
        super(DistortionModel, self).validate()
        assert isinstance(self.meta.input_units, (str, u.NamedUnit))
        assert isinstance(self.meta.output_units, (str, u.NamedUnit))
        if self.meta.instrument.name == 'NIRCAM':
            assert self.meta.instrument.module is not None
            assert self.meta.instrument.channel is not None
            assert self.meta.instrument.p_pupil is not None


class DistortionMRSModel(ReferenceFileModel):
    """
    A model for a reference file of type "distortion" for the MIRI MRS.
    """
    schema_url = "distortion_mrs.schema.yaml"
    reftype = "distortion"

    def __init__(self, init=None, x_model=None, y_model=None, alpha_model=None, beta_model=None,
                 bzero=None, bdel=None, input_units=None, output_units=None, **kwargs):

        super(DistortionMRSModel, self).__init__(init=init, **kwargs)

        if x_model is not None:
            self.x_model = x_model
        if y_model is not None:
            self.y_model = y_model
        if alpha_model is not None:
            self.alpha_model = alpha_model
        if beta_model is not None:
            self.beta_model = beta_model
        if bzero is not None:
            self.bzero = bzero
        if bdel is not None:
            self.bdel = bdel
        if input_units is not None:
            self.meta.input_units = input_units
        if output_units is not None:
            self.meta.output_units = output_units
        if init is None:
            try:
                self.populate_meta()
            except NotImplementedError:
                pass

    def on_save(self, path=None):
        self.meta.reftype = self.reftype

    def populate_meta(self):
        self.meta.instrument.name = "MIRI"
        self.meta.exposure.type = "MIR_MRS"
        self.meta.input_units = u.pix
        self.meta.output_units = u.arcsec

    def to_fits(self):
        raise NotImplementedError("FITS format is not supported for this file.")

    def validate(self):
        super(DistortionMRSModel, self).validate()
        assert isinstance(self.meta.input_units, (str, u.NamedUnit))
        assert isinstance(self.meta.output_units, (str, u.NamedUnit))
        assert self.meta.instrument.name == "MIRI"
        assert self.meta.exposure.type == "MIR_MRS"
        assert self.meta.instrument.channel in ("12", "34", "1", "2", "3", "4")
        assert self.meta.instrument.band in ("SHORT", "LONG", "MEDIUM")
        assert self.meta.instrument.detector in ("MIRIFUSHORT", "MIRIFULONG")
        assert all([isinstance(m, Model) for m in self.x_model])
        assert all([isinstance(m, Model) for m in self.y_model])
        assert all([isinstance(m, Model) for m in self.alpha_model])
        assert all([isinstance(m, Model) for m in self.beta_model])
        assert len(self.abv2v3_model.model) == 2
        assert len(self.abv2v3_model.channel_band) == 2


class SpecwcsModel(_SimpleModel):
    """
    A model for a reference file of type "specwcs".
    """
    schema_url = "specwcs.schema.yaml"
    reftype = "specwcs"

    def validate(self):
        assert isinstance(self.meta.input_units, (str, u.NamedUnit))
        assert isinstance(self.meta.output_units, (str, u.NamedUnit))
        assert self.meta.instrument.name in ["NIRCAM",
                                             "NIRSPEC",
                                             "MIRI",
                                             "TFI",
                                             "FGS",
                                             "NIRISS"]


class NIRCAMGrismModel(ReferenceFileModel):
    """
    A model for a reference file of type "specwcs" for NIRCAM grisms.

    This reference file contains the models for wave, x, and y polynomial
    solutions that describe dispersion through the grism
    """
    schema_url = "specwcs_nircam_grism.schema.yaml"
    reftype = "specwcs"

    def __init__(self, init=None,
                       displ=None,
                       dispx=None,
                       dispy=None,
                       invdispl=None,
                       invdispx=None,
                       invdispy=None,
                       orders=None,
                       **kwargs):
        super(NIRCAMGrismModel, self).__init__(init=init, **kwargs)

        if init is None:
            self.populate_meta()
        if displ is not None:
            self.displ = displ
        if dispx is not None:
            self.dispx = dispx
        if dispy is not None:
            self.dispy = dispy
        if invdispl is not None:
            self.invdispl = invdispl
        if invdispx is not None:
            self.invdispx = invdispx
        if invdispy is not None:
            self.invdispy = invdispy
        if orders is not None:
            self.orders = orders

    def populate_meta(self):
        self.meta.instrument.name = "NIRCAM"
        self.meta.exposure.type = "NRC_WFSS"
        self.meta.reftype = self.reftype

    def validate(self):
        assert isinstance(self.meta.input_units, (str, u.NamedUnit))
        assert isinstance(self.meta.output_units, (str, u.NamedUnit))
        assert self.meta.instrument.name == "NIRCAM"
        assert self.meta.exposure.type == "NRC_WFSS"
        assert self.meta.reftype == self.reftype

    def to_fits(self):
        raise NotImplementedError("FITS format is not supported for this file.")


class NIRISSGrismModel(ReferenceFileModel):
    """
    A model for a reference file of type "specwcs" for NIRISS grisms.
    """
    schema_url = "specwcs_niriss_grism.schema.yaml"
    reftype = "specwcs"

    def __init__(self, init=None,
                       displ=None,
                       dispx=None,
                       dispy=None,
                       invdispl=None,
                       orders=None,
                       fwcpos_ref=None,
                       **kwargs):
        super(NIRISSGrismModel, self).__init__(init=init, **kwargs)

        if init is None:
            self.populate_meta()
        if displ is not None:
            self.displ = displ
        if dispx is not None:
            self.dispx = dispx
        if dispy is not None:
            self.dispy = dispy
        if invdispl is not None:
            self.invdispl = invdispl
        if orders is not None:
            self.orders = orders
        if fwcpos_ref is not None:
            self.fwcpos_ref = fwcpos_ref

    def populate_meta(self):
        self.meta.instrument.name = "NIRISS"
        self.meta.instrument.detector = "NIS"
        self.meta.exposure.type = "NIS_WFSS"
        self.meta.reftype = self.reftype

    def validate(self):
        assert isinstance(self.meta.input_units, (str, u.NamedUnit))
        assert isinstance(self.meta.output_units, (str, u.NamedUnit))
        assert self.meta.instrument.name == "NIRISS"
        assert self.meta.exposure.type == "NIS_WFSS"
        assert self.meta.instrument.detector == "NIS"
        assert self.meta.reftype == self.reftype

    def to_fits(self):
        raise NotImplementedError("FITS format is not supported for this file.")


class RegionsModel(ReferenceFileModel):
    """
    A model for a reference file of type "regions".
    """
    schema_url = "regions.schema.yaml"
    reftype = "regions"

    def __init__(self, init=None, regions=None, **kwargs):
        super(RegionsModel, self).__init__(init=init, **kwargs)
        if regions is not None:
            self.regions = regions
        if init is None:
            self.populate_meta()

    def on_save(self, path=None):
        self.meta.reftype = self.reftype

    def populate_meta(self):
        self.meta.instrument.name = "MIRI"
        self.meta.exposure.type = "MIR_MRS"

    def to_fits(self):
        raise NotImplementedError("FITS format is not supported for this file.")

    def validate(self):
        super(RegionsModel, self).validate()
        assert isinstance(self.regions.copy(), np.ndarray)
        assert self.meta.instrument.name == "MIRI"
        assert self.meta.exposure.type == "MIR_MRS"
        assert self.meta.instrument.channel in ("12", "34", "1", "2", "3", "4")
        assert self.meta.instrument.band in ("SHORT", "LONG")
        assert self.meta.instrument.detector in ("MIRIFUSHORT", "MIRIFULONG")


class WavelengthrangeModel(ReferenceFileModel):
    """
    A model for a reference file of type "wavelengthrange".
    The model is used by MIRI, NIRSPEC, NIRCAM, and NIRISS
    """
    schema_url = "wavelengthrange.schema.yaml"
    reftype = "wavelengthrange"

    def __init__(self, init=None, wrange_selector=None, wrange=None, order=None, wunits=None, **kwargs):

        super(WavelengthrangeModel, self).__init__(init=init, **kwargs)
        if wrange_selector is not None:
            self.waverange_selector = wrange_selector
        if wrange is not None:
            self.wavelengthrange = wrange
        if order is not None:
            self.order = order
        if wunits is not None:
            self.meta.wavelength_units = wunits

    def on_save(self, path=None):
        self.meta.reftype = self.reftype

    def to_fits(self):
        raise NotImplementedError("FITS format is not supported for this file")

    def validate(self):
        super(WavelengthrangeModel, self).validate()
        assert self.meta.instrument.name in ("MIRI", "NIRSPEC", "NIRCAM", "NIRISS")


class FPAModel(ReferenceFileModel):
    """
    A model for a NIRSPEC reference file of type "fpa".
    """
    schema_url = "fpa.schema.yaml"
    reftype = "fpa"

    def __init__(self, init=None, nrs1_model=None, nrs2_model=None, **kwargs):

        super(FPAModel, self).__init__(init=init, **kwargs)
        if nrs1_model is not None:
            self.nrs1_model = nrs1_model
        if nrs2_model is not None:
            self.nrs2_model = nrs2_model
        if init is None:
            self.populate_meta()

    def on_save(self, path=None):
        self.meta.reftype = self.reftype

    def populate_meta(self):
        self.meta.instrument.name = "NIRSPEC"
        self.meta.instrument.p_detector = "NRS1|NRS2|"
        self.meta.exposure.p_exptype = "NRS_TACQ|NRS_TASLIT|NRS_TACONFIRM|\
        NRS_CONFIRM|NRS_FIXEDSLIT|NRS_IFU|NRS_MSASPEC|NRS_IMAGE|NRS_FOCUS|\
        NRS_MIMF|NRS_BOTA|NRS_LAMP|NRS_BRIGHTOBJ|"
        self.meta.exposure.type = "N/A"

    def to_fits(self):
        raise NotImplementedError("FITS format is not supported for this file.")

    def validate(self):
        super(FPAModel, self).validate()
        assert isinstance(self.nrs1_model, Model)
        assert isinstance(self.nrs2_model, Model)


class IFUPostModel(ReferenceFileModel):
    """
    A model for a NIRSPEC reference file of type "ifupost".

    Parameters
    ----------
    init : str
        A file name.
    slice_models : dict
        A dictionary with slice transforms with the following entries
        {"slice_N":
        {'linear': astropy.modeling.Model,
        'xpoly': astropy.modeling.Model,
        'xpoly_distortion': astropy.modeling.Model,
        'ypoly': astropy.modeling.Model,
        'ypoly_distortion': astropy.modeling.Model,
        }
        }
    """
    schema_url = "ifupost.schema.yaml"
    reftype = "ifupost"

    def __init__(self, init=None, slice_models=None, **kwargs):

        super(IFUPostModel, self).__init__(init=init, **kwargs)
        if slice_models is not None:
            if len(slice_models) != 30:
                raise ValueError("Expected 30 slice models, got {0}".format(len(slice_models)))
            else:
                for key, val in slice_models.items():
                    setattr(self, key, val)
        if init is None:
            self.populate_meta()

    def on_save(self, path=None):
        self.meta.reftype = self.reftype

    def populate_meta(self):
        self.meta.instrument.name = "NIRSPEC"
        self.meta.instrument.p_detector = "NRS1|NRS2|"
        self.meta.exposure.type = "NRS_IFU"
        self.meta.exposure.p_exptype = "NRS_IFU"

    def to_fits(self):
        raise NotImplementedError("FITS format is not supported for this file.")

    def validate(self):
        super(IFUPostModel, self).validate()


class IFUSlicerModel(ReferenceFileModel):
    """
    A model for a NIRSPEC reference file of type "ifuslicer".
    """
    schema_url = "ifuslicer.schema.yaml"
    reftype = "ifuslicer"

    def __init__(self, init=None, model=None, data=None, **kwargs):

        super(IFUSlicerModel, self).__init__(init=init, **kwargs)
        if model is not None:
            self.model = model
        if data is not None:
            seld.data = data
        if init is None:
            self.populate_meta()

    def on_save(self, path=None):
        self.meta.reftype = self.reftype

    def populate_meta(self):
        self.meta.instrument.name = "NIRSPEC"
        self.meta.instrument.p_detector = "NRS1|NRS2|"
        self.meta.exposure.type = "NRS_IFU"
        self.meta.exposure.p_exptype = "NRS_IFU"

    def to_fits(self):
        raise NotImplementedError("FITS format is not supported for this file.")

    def validate(self):
        super(IFUSlicerModel, self).validate()


class MSAModel(ReferenceFileModel):
    """
    A model for a NIRSPEC reference file of type "msa".
    """
    schema_url = "msa.schema.yaml"
    reftype = "msa"

    def __init__(self, init=None, models=None, data=None, **kwargs):
        super(MSAModel, self).__init__(init=init, **kwargs)
        if models is not None and data is not None:
            self.Q1 = {'model': models['Q1'], 'data': data['Q1']}
            self.Q2 = {'model': models['Q2'], 'data': data['Q2']}
            self.Q3 = {'model': models['Q3'], 'data': data['Q3']}
            self.Q4 = {'model': models['Q4'], 'data': data['Q4']}
            self.Q5 = {'model': models['Q5'], 'data': data['Q5']}
        if init is None:
            self.populate_meta()

    def on_save(self, path=None):
        self.meta.reftype = self.reftype

    def populate_meta(self):
        self.meta.instrument.name = "NIRSPEC"
        self.meta.instrument.p_detector = "NRS1|NRS2|"
        self.meta.exposure.p_exptype = "NRS_TACQ|NRS_TASLIT|NRS_TACONFIRM|\
        NRS_CONFIRM|NRS_FIXEDSLIT|NRS_IFU|NRS_MSASPEC|NRS_IMAGE|NRS_FOCUS|\
        NRS_MIMF|NRS_BOTA|NRS_LAMP|NRS_BRIGHTOBJ|"
        self.meta.exposure.type = "N/A"

    def to_fits(self):
        raise NotImplementedError("FITS format is not supported for this file.")

    def validate(self):
        super(MSAModel, self).validate()


class DisperserModel(ReferenceFileModel):
    """
    A model for a NIRSPEC reference file of type "disperser".
    """
    schema_url = "disperser.schema.yaml"
    reftype = "disperser"

    def __init__(self, init=None, angle=None, gwa_tiltx=None, gwa_tilty=None,
                 kcoef=None, lcoef=None, tcoef=None, pref=None, tref=None,
                 theta_x=None, theta_y=None,theta_z=None,
                 groovedensity=None, **kwargs):
        super(DisperserModel, self).__init__(init=init, **kwargs)
        if groovedensity is not None:
            self.groovedensity = groovedensity
        if angle is not None:
            self.angle = angle
        if gwa_tiltx is not None:
            self.gwa_tiltx = gwa_tiltx
        if gwa_tilty is not None:
            self.gwa_tilty = gwa_tilty
        if kcoef is not None:
            self.kcoef = kcoef
        if lcoef is not None:
            self.lcoef = lcoef
        if tcoef is not None:
            self.tcoef = tcoef
        if pref is not None:
            self.pref = pref
        if tref is not None:
            self.tref = tref
        if theta_x is not None:
            self.theta_x = theta_x
        if theta_y is not None:
            self.theta_y = theta_y
        if theta_z is not None:
            self.theta_z = theta_z
        if gwa_tiltx is not None:
            self.gwa_tiltx = gwa_tiltx
        if gwa_tilty is not None:
            self.gwa_tilty = gwa_tilty
        if init is None:
            self.populate_meta()

    def on_save(self, path=None):
        self.meta.reftype = self.reftype

    def populate_meta(self):
        self.meta.instrument.name = "NIRSPEC"
        self.meta.instrument.p_detector = "NRS1|NRS2|"
        self.meta.exposure.p_exptype = "NRS_TACQ|NRS_TASLIT|NRS_TACONFIRM|\
        NRS_CONFIRM|NRS_FIXEDSLIT|NRS_IFU|NRS_MSASPEC|NRS_IMAGE|NRS_FOCUS|\
        NRS_MIMF|NRS_BOTA|NRS_LAMP|NRS_BRIGHTOBJ|"
        self.meta.exposure.type = "N/A"

    def to_fits(self):
        raise NotImplementedError("FITS format is not supported for this file.")

    def validate(self):
        super(DisperserModel, self).validate()
        assert self.meta.instrument.grating in ["G140H", "G140M", "G235H", "G235M",
                                                "G395H", "G395M", "MIRROR", "PRISM"]


class FilteroffsetModel(ReferenceFileModel):
    """
    A model for a NIRSPEC reference file of type "disperser".
    """
    schema_url = "filteroffset.schema.yaml"
    reftype = "filteroffset"

    def __init__(self, init=None, filters=None, **kwargs):
        super(FilteroffsetModel, self).__init__(init, **kwargs)
        if filters is not None:
            self.filters = filters

    def populate_meta(self):
        self.meta.instrument.name = "MIRI"
        self.meta.instrument.detector = "MIRIMAGE"
        self.meta.instrument.pfilter = "F1130W|F1140C|F2300C|F2100W|F1800W|\
        F1550C|F560W|F2550WR|FND|F2550W|F1500W|F1000W|F1065C|F770W|F1280W|"

    def on_save(self, path=None):
        self.meta.reftype = self.reftype

    def validate(self):
        assert self.meta.instrument.name == "MIRI"
        assert self.meta.instrument.detector == "MIRIMAGE"
        super(FilteroffsetModel, self).validate()


class IFUFOREModel(_SimpleModel):
    """
    A model for a NIRSPEC reference file of type "ifufore".
    """
    schema_url = "ifufore.schema.yaml"
    reftype = "ifufore"

    def populate_meta(self):
        self.meta.instrument.name = "NIRSPEC"
        self.meta.instrument.p_detector = "NRS1|NRS2|"
        self.meta.exposure.type = "NRS_IFU"
        self.meta.exposure.p_exptype = "NRS_IFU"


class CameraModel(_SimpleModel):
    """
    A model for a reference file of type "camera".
    """
    schema_url = "camera.schema.yaml"
    reftype = 'camera'

    def populate_meta(self):
        self.meta.instrument.name = "NIRSPEC"
        self.meta.instrument.p_detector = "NRS1|NRS2|"
        self.meta.exposure.p_exptype = "NRS_TACQ|NRS_TASLIT|NRS_TACONFIRM|\
        NRS_CONFIRM|NRS_FIXEDSLIT|NRS_IFU|NRS_MSASPEC|NRS_IMAGE|NRS_FOCUS|\
        NRS_MIMF|NRS_BOTA|NRS_LAMP|NRS_BRIGHTOBJ|"
        self.meta.exposure.type = "N/A"


class CollimatorModel(_SimpleModel):
    """
    A model for a reference file of type "collimator".
    """
    schema_url = "collimator.schema.yaml"
    reftype = 'collimator'

    def populate_meta(self):
        self.meta.instrument.name = "NIRSPEC"
        self.meta.instrument.p_detector = "NRS1|NRS2|"
        self.meta.exposure.p_exptype = "NRS_TACQ|NRS_TASLIT|NRS_TACONFIRM|\
        NRS_CONFIRM|NRS_FIXEDSLIT|NRS_IFU|NRS_MSASPEC|NRS_IMAGE|NRS_FOCUS|\
        NRS_MIMF|NRS_BOTA|NRS_LAMP|NRS_BRIGHTOBJ|"
        self.meta.exposure.type = "N/A"


class OTEModel(_SimpleModel):
    """
    A model for a reference file of type "ote".
    """
    schema_url = "ote.schema.yaml"
    reftype = 'ote'

    def populate_meta(self):
        self.meta.instrument.name = "NIRSPEC"
        self.meta.instrument.p_detector = "NRS1|NRS2|"
        self.meta.exposure.p_exptype = "NRS_TACQ|NRS_TASLIT|NRS_TACONFIRM|\
        NRS_CONFIRM|NRS_FIXEDSLIT|NRS_IFU|NRS_MSASPEC|NRS_IMAGE|NRS_FOCUS|\
        NRS_MIMF|NRS_BOTA|NRS_LAMP|NRS_BRIGHTOBJ|"
        self.meta.exposure.type = "N/A"


class FOREModel(_SimpleModel):
    """
    A model for a reference file of type "fore".
    """
    schema_url = "fore.schema.yaml"
    reftype = 'fore'

    def populate_meta(self):
        self.meta.instrument.name = "NIRSPEC"
        self.meta.instrument.p_detector = "NRS1|NRS2|"
        self.meta.exposure.p_exptype = "NRS_TACQ|NRS_TASLIT|NRS_TACONFIRM|\
        NRS_CONFIRM|NRS_FIXEDSLIT|NRS_IFU|NRS_MSASPEC|NRS_IMAGE|NRS_FOCUS|\
        NRS_MIMF|NRS_BOTA|NRS_LAMP|NRS_BRIGHTOBJ|"
        self.meta.exposure.type = "N/A"

    def on_save(self, path=None):
        self.meta.reftype = self.reftype

    def validate(self):
        super(FOREModel, self).validate()
        assert self.meta.instrument.filter in ["CLEAR", "F070LP", "F100LP", "F110W",
                                               "F140X", "F170LP", "F290LP"]


class WaveCorrModel(ReferenceFileModel):

    reftype = "wavecorr"
    schema_url = "wavecorr.schema.yaml"

    def __init__(self, init=None, apertures=None, **kwargs):
        super(WaveCorrModel, self).__init__(init, **kwargs)
        if apertures is not None:
            self.apertures = apertures
        if init is None:
            self.populate_meta()

    @property
    def aperture_names(self):
        return [getattr(ap, 'aperture_name') for ap in self.apertures]

    def populate_meta(self):
        self.meta.instrument.name = "NIRSPEC"

    def on_save(self, path=None):
        self.meta.reftype = self.reftype

    def validate(self):
        super(WaveCorrModel, self).validate()
        assert self.aperture_names is not None
        assert self.apertures is not None
