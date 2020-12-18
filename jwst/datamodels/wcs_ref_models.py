import traceback
import warnings

import numpy as np
from astropy.modeling.core import Model
from astropy import units as u
from stdatamodels.validate import ValidationWarning

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
        super().validate()
        try:
            assert isinstance(self.model, Model)
            assert self.meta.instrument.name in ["NIRCAM", "NIRSPEC", "MIRI", "TFI", "FGS", "NIRISS"]
        except AssertionError:
            if self._strict_validation:
                raise
            else:
                warnings.warn(traceback.format_exc(), ValidationWarning)

class DistortionModel(_SimpleModel):
    """
    A model for a reference file of type "distortion".
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/distortion.schema"
    reftype = "distortion"

    def validate(self):
        super().validate()
        try:
            assert isinstance(self.meta.input_units, (str, u.NamedUnit))
            assert isinstance(self.meta.output_units, (str, u.NamedUnit))
            if self.meta.instrument.name == 'NIRCAM':
                assert self.meta.instrument.module is not None
                assert self.meta.instrument.channel is not None
                assert self.meta.instrument.p_pupil is not None
        except AssertionError:
            if self._strict_validation:
                raise
            else:
                warnings.warn(traceback.format_exc(), ValidationWarning)

class DistortionMRSModel(ReferenceFileModel):
    """
    A model for a reference file of type "distortion" for the MIRI MRS.
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/distortion_mrs.schema"
    reftype = "distortion"

    def __init__(self, init=None, x_model=None, y_model=None, alpha_model=None, beta_model=None,
                 bzero=None, bdel=None, input_units=None, output_units=None, **kwargs):

        super().__init__(init=init, **kwargs)

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
        super().validate()
        try:
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
        except AssertionError:
            if self._strict_validation:
                raise
            else:
                warnings.warn(traceback.format_exc(), ValidationWarning)


class SpecwcsModel(_SimpleModel):
    """
    A model for a reference file of type "specwcs".

    Notes
    -----
    For NIRISS and NIRCAM WFSS modes the specwcs file is
    used during extract_2D. See NIRCAMGrismModel and
    NIRISSGrismModel.
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/specwcs.schema"
    reftype = "specwcs"

    def validate(self):
        super().validate()
        try:
            assert isinstance(self.meta.input_units, (str, u.NamedUnit))
            assert isinstance(self.meta.output_units, (str, u.NamedUnit))
            assert self.meta.instrument.name in ["NIRCAM",
                                                 "NIRSPEC",
                                                 "MIRI",
                                                 "TFI",
                                                 "FGS",
                                                 "NIRISS"]
        except AssertionError:
            if self._strict_validation:
                raise
            else:
                warnings.warn(traceback.format_exc(), ValidationWarning)


class NIRCAMGrismModel(ReferenceFileModel):
    """
    A model for a reference file of type "specwcs" for NIRCAM WFSS.

    This reference file contains the models for wave, x, and y polynomial
    solutions that describe dispersion through the grism.

    Parameters
    ----------
    displ: `~astropy.modeling.Model`
          Nircam Grism wavelength dispersion model
    dispx : `~astropy.modeling.Model`
          Nircam Grism row dispersion model
    dispy : `~astropy.modeling.Model`
          Nircam Grism column dispersion model
    invdispl : `~astropy.modeling.Model`
          Nircam Grism inverse wavelength dispersion model
    invdispx : `~astropy.modeling.Model`
          Nircam Grism inverse row dispersion model
    invdispy : `~astropy.modeling.Model`
          Nircam Grism inverse column dispersion model
    orders : `~astropy.modeling.Model`
          NIRCAM Grism orders, matched to the array locations of the
          dispersion models

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/specwcs_nircam_grism.schema"
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
        super().__init__(init=init, **kwargs)

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
        super().validate()
        try:
            assert isinstance(self.meta.input_units, (str, u.NamedUnit))
            assert isinstance(self.meta.output_units, (str, u.NamedUnit))
            assert self.meta.instrument.name == "NIRCAM"
            assert self.meta.exposure.type == "NRC_WFSS"
            assert self.meta.reftype == self.reftype
        except AssertionError:
            if self._strict_validation:
                raise
            else:
                warnings.warn(traceback.format_exc(), ValidationWarning)

    def to_fits(self):
        raise NotImplementedError("FITS format is not supported for this file.")


class NIRISSGrismModel(ReferenceFileModel):
    """
    A model for a reference file of type "specwcs" for NIRISS grisms.

    Parameters
    ----------
    displ: `~astropy.modeling.Model`
          NIRISS Grism wavelength dispersion model
    dispx : `~astropy.modeling.Model`
          NIRISS Grism row dispersion model
    dispy : `~astropy.modeling.Model`
          NIRISS Grism column dispersion model
    invdispl : `~astropy.modeling.Model`
          NIRISS Grism inverse wavelength dispersion model
    invdispx : `~astropy.modeling.Model`
          NIRISS Grism inverse row dispersion model
    invdispy : `~astropy.modeling.Model`
          NIRISS Grism inverse column dispersion model
    orders : `~astropy.modeling.Model`
          NIRISS Grism orders, matched to the array locations of the
          dispersion models
    fwcpos_ref : float
        The reference value for the filter wheel position
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/specwcs_niriss_grism.schema"
    reftype = "specwcs"

    def __init__(self, init=None,
                       displ=None,
                       dispx=None,
                       dispy=None,
                       invdispl=None,
                       orders=None,
                       fwcpos_ref=None,
                       **kwargs):
        super().__init__(init=init, **kwargs)

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
        super(NIRISSGrismModel, self).validate()
        try:
            assert isinstance(self.meta.input_units, (str, u.NamedUnit))
            assert isinstance(self.meta.output_units, (str, u.NamedUnit))
            assert self.meta.instrument.name == "NIRISS"
            assert self.meta.exposure.type == "NIS_WFSS"
            assert self.meta.instrument.detector == "NIS"
            assert self.meta.reftype == self.reftype
        except AssertionError:
            if self._strict_validation:
                raise
            else:
                warnings.warn(traceback.format_exc(), ValidationWarning)

    def to_fits(self):
        raise NotImplementedError("FITS format is not supported for this file.")


class RegionsModel(ReferenceFileModel):
    """
    A model for a reference file of type "regions".
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/regions.schema"
    reftype = "regions"

    def __init__(self, init=None, regions=None, **kwargs):
        super().__init__(init=init, **kwargs)
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
        super().validate()
        try:
            assert isinstance(self.regions, np.ndarray)
            assert self.meta.instrument.name == "MIRI"
            assert self.meta.exposure.type == "MIR_MRS"
            assert self.meta.instrument.channel in ("12", "34", "1", "2", "3", "4")
            assert self.meta.instrument.band in ("SHORT", "LONG")
            assert self.meta.instrument.detector in ("MIRIFUSHORT", "MIRIFULONG")
        except AssertionError:
            if self._strict_validation:
                raise
            else:
                warnings.warn(traceback.format_exc(), ValidationWarning)


class WavelengthrangeModel(ReferenceFileModel):
    """
    A model for a reference file of type "wavelengthrange".

    The model is used by MIRI, NIRSPEC, NIRCAM, and NIRISS.


    Parameters
    ----------
    wrange : list
        Contains a list of [order, filter, min wave, max wave]
    order : list
        A list of orders that are available and described in the file
    extract_orders : list
        A list of filters and the orders that should be extracted by default
    wunits : `~astropy.units`
        The units for the wavelength data

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/wavelengthrange.schema"
    reftype = "wavelengthrange"

    def __init__(self, init=None, wrange_selector=None, wrange=None,
                 order=None, extract_orders=None, wunits=None, **kwargs):

        super().__init__(init=init, **kwargs)
        if wrange_selector is not None:
            self.waverange_selector = wrange_selector
        if wrange is not None:
            self.wavelengthrange = wrange
        if order is not None:
            self.order = order
        if extract_orders is not None:
            self.extract_orders = extract_orders
        if wunits is not None:
            self.meta.wavelength_units = wunits

    def on_save(self, path=None):
        self.meta.reftype = self.reftype

    def to_fits(self):
        raise NotImplementedError("FITS format is not supported for this file")

    def validate(self):
        super().validate()
        try:
            assert self.meta.instrument.name in ("MIRI", "NIRSPEC", "NIRCAM", "NIRISS")
        except AssertionError:
            if self._strict_validation:
                raise
            else:
                warnings.warn(traceback.format_exc(), ValidationWarning)

    def get_wfss_wavelength_range(self, filter, orders):
        """ Retrieve the wavelength range for a WFSS observation.

        Parameters
        ----------
        filter : str
            Filter for which to retrieve the wavelength range.
        orders : list
            List of spectral orders

        Returns
        -------
        wave_range : dict
            Pairs of {order: (wave_min, wave_max)} for each order and the specific filter.
        """
        wave_range = {}
        for order in orders:
            range_select = [(x[2], x[3]) for x in self.wavelengthrange if (x[0] == order and x[1] == filter)]
            wave_range[order] = range_select[0]
        return wave_range


class FPAModel(ReferenceFileModel):
    """
    A model for a NIRSPEC reference file of type "fpa".
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/fpa.schema"
    reftype = "fpa"

    def __init__(self, init=None, nrs1_model=None, nrs2_model=None, **kwargs):

        super().__init__(init=init, **kwargs)
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
        NRS_MIMF|NRS_MSATA|NRS_WATA|NRS_LAMP|NRS_BRIGHTOBJ|"
        self.meta.exposure.type = "N/A"

    def to_fits(self):
        raise NotImplementedError("FITS format is not supported for this file.")

    def validate(self):
        super().validate()
        try:
            assert isinstance(self.nrs1_model, Model)
            assert isinstance(self.nrs2_model, Model)
        except AssertionError:
            if self._strict_validation:
                raise
            else:
                warnings.warn(traceback.format_exc(), ValidationWarning)


class IFUPostModel(ReferenceFileModel):
    """
    A model for a NIRSPEC reference file of type "ifupost".

    Parameters
    ----------
    init : str
        A file name.
    slice_models : dict
        A dictionary with slice transforms with the following entries
        {"slice_N": {'linear': astropy.modeling.Model,
        ...         'xpoly': astropy.modeling.Model,
        ...         'xpoly_distortion': astropy.modeling.Model,
        ...         'ypoly': astropy.modeling.Model,
        ...         'ypoly_distortion': astropy.modeling.Model,
        ...         }}

    """

    schema_url = "http://stsci.edu/schemas/jwst_datamodel/ifupost.schema"
    reftype = "ifupost"

    def __init__(self, init=None, slice_models=None, **kwargs):

        super().__init__(init=init, **kwargs)
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
        super().validate()


class IFUSlicerModel(ReferenceFileModel):
    """
    A model for a NIRSPEC reference file of type "ifuslicer".
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/ifuslicer.schema"
    reftype = "ifuslicer"

    def __init__(self, init=None, model=None, data=None, **kwargs):

        super().__init__(init=init, **kwargs)
        if model is not None:
            self.model = model
        if data is not None:
            self.data = data
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
        super().validate()


class MSAModel(ReferenceFileModel):
    """
    A model for a NIRSPEC reference file of type "msa".
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/msa.schema"
    reftype = "msa"

    def __init__(self, init=None, models=None, data=None, **kwargs):
        super().__init__(init=init, **kwargs)
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
        NRS_MIMF|NRS_MSATA|NRS_WATA|NRS_LAMP|NRS_BRIGHTOBJ|"
        self.meta.exposure.type = "N/A"

    def to_fits(self):
        raise NotImplementedError("FITS format is not supported for this file.")

    def validate(self):
        super().validate()


class DisperserModel(ReferenceFileModel):
    """
    A model for a NIRSPEC reference file of type "disperser".
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/disperser.schema"
    reftype = "disperser"

    def __init__(self, init=None, angle=None, gwa_tiltx=None, gwa_tilty=None,
                 kcoef=None, lcoef=None, tcoef=None, pref=None, tref=None,
                 theta_x=None, theta_y=None,theta_z=None,
                 groovedensity=None, **kwargs):
        super().__init__(init=init, **kwargs)
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
        NRS_MIMF|NRS_MSATA|NRS_WATA|NRS_LAMP|NRS_BRIGHTOBJ|"
        self.meta.exposure.type = "N/A"

    def to_fits(self):
        raise NotImplementedError("FITS format is not supported for this file.")

    def validate(self):
        super().validate()
        try:
            assert self.meta.instrument.grating in ["G140H", "G140M", "G235H", "G235M",
                                                    "G395H", "G395M", "MIRROR", "PRISM"]
        except AssertionError:
            if self._strict_validation:
                raise
            else:
                warnings.warn(traceback.format_exc(), ValidationWarning)


class FilteroffsetModel(ReferenceFileModel):
    """
    A model for a NIRSPEC reference file of type "disperser".
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/filteroffset.schema"
    reftype = "filteroffset"

    def __init__(self, init=None, filters=None, instrument=None, **kwargs):
        super().__init__(init, **kwargs)
        if filters is not None:
            self.filters = filters
            if instrument is None or instrument not in ("NIRCAM", "MIRI", "NIRISS"):
                raise ValueError("Please specify the instrument")
            self.meta.instrument.name = instrument

    def populate_meta(self):
        self.meta.reftype = self.reftype
        if self.meta.instrument.name == "MIRI":
            self.meta.instrument.detector = "MIRIMAGE"
            self.meta.instrument.pfilter = "F1130W|F1140C|F2300C|F2100W|F1800W|\
            F1550C|F560W|F2550WR|FND|F2550W|F1500W|F1000W|F1065C|F770W|F1280W|"
        elif self.meta.instrument.name == "NIRCAM":
            self.meta.instrument.pfilter = "F070W|F090W|F115W|F140M|F150W|\
            F162M|F164N|F150W2|F182M|F187N|F200W|F210M|F212N|"
        elif self.meta.instrument.name == "NIRISS":
            self.meta.instrument.pfilter = "F070W|F115W|F140M|F150W|F158M|\
            F200W|F277W|F356W|F380M|F430M|F444W|F480M|"
        else:
            raise ValueError(f"Unsupported instrument: {self.meta.instrument.name}")

    def validate(self):
        super().validate()

        instrument_name = self.meta.instrument.name
        nircam_channels = ["SHORT", "LONG"]
        nircam_module = ["A", "B"]
        if instrument_name not in ['MIRI', 'NIRCAM', 'NIRISS']:
            self.print_err('Expected "meta.instrument.name" to be one of "NIRCAM, "MIRI" or "NIRISS"')
        if instrument_name == "MIRI" and self.meta.instrument.detector != "MIRIMAGE":
            self.print_err('Expected detector to be MIRIMAGE for instrument MIRI')
        elif instrument_name == "NIRCAM":
            if self.meta.instrument.channel not in nircam_channels:
                self.print_err(f'Expected meta.instrument.channel for instrument NIRCAM to be one of {nircam_channels}')
            if self.meta.instrument.module not in nircam_module:
                self.print_err(f'Expected meta.instrument.module for instrument NIRCAM to be one of {nircam_module}')


class IFUFOREModel(_SimpleModel):
    """
    A model for a NIRSPEC reference file of type "ifufore".
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/ifufore.schema"
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
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/camera.schema"
    reftype = 'camera'

    def populate_meta(self):
        self.meta.instrument.name = "NIRSPEC"
        self.meta.instrument.p_detector = "NRS1|NRS2|"
        self.meta.exposure.p_exptype = "NRS_TACQ|NRS_TASLIT|NRS_TACONFIRM|\
        NRS_CONFIRM|NRS_FIXEDSLIT|NRS_IFU|NRS_MSASPEC|NRS_IMAGE|NRS_FOCUS|\
        NRS_MIMF|NRS_MSATA|NRS_WATA|NRS_LAMP|NRS_BRIGHTOBJ|"
        self.meta.exposure.type = "N/A"


class CollimatorModel(_SimpleModel):
    """
    A model for a reference file of type "collimator".
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/collimator.schema"
    reftype = 'collimator'

    def populate_meta(self):
        self.meta.instrument.name = "NIRSPEC"
        self.meta.instrument.p_detector = "NRS1|NRS2|"
        self.meta.exposure.p_exptype = "NRS_TACQ|NRS_TASLIT|NRS_TACONFIRM|\
        NRS_CONFIRM|NRS_FIXEDSLIT|NRS_IFU|NRS_MSASPEC|NRS_IMAGE|NRS_FOCUS|\
        NRS_MIMF|NRS_MSATA|NRS_WATA|NRS_LAMP|NRS_BRIGHTOBJ|"
        self.meta.exposure.type = "N/A"


class OTEModel(_SimpleModel):
    """
    A model for a reference file of type "ote".
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/ote.schema"
    reftype = 'ote'

    def populate_meta(self):
        self.meta.instrument.name = "NIRSPEC"
        self.meta.instrument.p_detector = "NRS1|NRS2|"
        self.meta.exposure.p_exptype = "NRS_TACQ|NRS_TASLIT|NRS_TACONFIRM|\
        NRS_CONFIRM|NRS_FIXEDSLIT|NRS_IFU|NRS_MSASPEC|NRS_IMAGE|NRS_FOCUS|\
        NRS_MIMF|NRS_MSATA|NRS_WATA|NRS_LAMP|NRS_BRIGHTOBJ|"
        self.meta.exposure.type = "N/A"


class FOREModel(_SimpleModel):
    """
    A model for a reference file of type "fore".
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/fore.schema"
    reftype = 'fore'

    def populate_meta(self):
        self.meta.instrument.name = "NIRSPEC"
        self.meta.instrument.p_detector = "NRS1|NRS2|"
        self.meta.exposure.p_exptype = "NRS_TACQ|NRS_TASLIT|NRS_TACONFIRM|\
        NRS_CONFIRM|NRS_FIXEDSLIT|NRS_IFU|NRS_MSASPEC|NRS_IMAGE|NRS_FOCUS|\
        NRS_MIMF|NRS_MSATA|NRS_WATA|NRS_LAMP|NRS_BRIGHTOBJ|"
        self.meta.exposure.type = "N/A"

    def on_save(self, path=None):
        self.meta.reftype = self.reftype

    def validate(self):
        super().validate()
        try:
            assert self.meta.instrument.filter in ["CLEAR", "F070LP", "F100LP", "F110W",
                                                   "F140X", "F170LP", "F290LP"]
        except AssertionError:
            if self._strict_validation:
                raise
            else:
                warnings.warn(traceback.format_exc(), ValidationWarning)


class WaveCorrModel(ReferenceFileModel):

    reftype = "wavecorr"
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/wavecorr.schema"

    def __init__(self, init=None, apertures=None, **kwargs):
        super().__init__(init, **kwargs)
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
        super().validate()
        try:
            assert self.aperture_names is not None
            assert self.apertures is not None
        except AssertionError:
            if self._strict_validation:
                raise
            else:
                warnings.warn(traceback.format_exc(), ValidationWarning)
