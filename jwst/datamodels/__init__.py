from astropy.io import registry

from . import ndmodel

from .model_base import DataModel
from .amilg import AmiLgModel
from .asn import AsnModel
from .barshadow import BarshadowModel
from .combinedspec import CombinedSpecModel
from .container import ModelContainer
from .contrast import ContrastModel
from .cube import CubeModel
from .dark import DarkModel
from .darkMIRI import DarkMIRIModel
from .drizpars import DrizParsModel
from .drizproduct import DrizProductModel
from .extract1dimage import Extract1dImageModel
from .flat import FlatModel
from .fringe import FringeModel
from .gain import GainModel
from .gls_rampfit import GLS_RampFitModel
from .guiderraw import GuiderRawModel
from .guidercal import GuiderCalModel
from .ifucube import IFUCubeModel
from .ifucubepars import IFUCubeParsModel, NirspecIFUCubeParsModel, MiriIFUCubeParsModel
from .ifuimage import IFUImageModel
from .image import ImageModel
from .ipc import IPCModel
from .irs2 import IRS2Model
from .lastframe import LastFrameModel
from .level1b import Level1bModel
from .linearity import LinearityModel
from .mask import MaskModel
from .miri_ramp import MIRIRampModel
from .multiexposure import MultiExposureModel
from .multiextract1d import MultiExtract1dImageModel
from .multiprod import MultiProductModel
from .multislit import MultiSlitModel
from .multispec import MultiSpecModel
from .nirspec_flat import NirspecFlatModel, NirspecQuadFlatModel
from .outlierpars import OutlierParsModel
from .pathloss import PathlossModel
from .persat import PersistenceSatModel
from .photom import PhotomModel, FgsPhotomModel, FgsImgPhotomModel, NircamPhotomModel, NrcImgPhotomModel
from .photom import NirissPhotomModel, NisImgPhotomModel, NirspecPhotomModel, NirspecFSPhotomModel
from .photom import MiriImgPhotomModel, MirImgPhotomModel, MiriMrsPhotomModel
from .pixelarea import PixelAreaModel, NirspecSlitAreaModel, NirspecMosAreaModel, NirspecIfuAreaModel
from .psfmask import PsfMaskModel
from .quad import QuadModel
from .ramp import RampModel
from .rampfitoutput import RampFitOutputModel
from .readnoise import ReadnoiseModel
from .reference import ReferenceFileModel, ReferenceImageModel, ReferenceCubeModel, ReferenceQuadModel
from .reset import ResetModel
from .resolution import ResolutionModel, MiriResolutionModel
from .rscd import RSCDModel
from .saturation import SaturationModel
from .slit import SlitModel, SlitDataModel
from .source_container import SourceModelContainer
from .spec import SpecModel
from .steppars import StepParsModel
from .straylight import StrayLightModel
from .superbias import SuperBiasModel
from .throughput import ThroughputModel
from .trapdensity import TrapDensityModel
from .trappars import TrapParsModel
from .trapsfilled import TrapsFilledModel
from .tsophot import TsoPhotModel
from .wcs_ref_models import (DistortionModel, DistortionMRSModel, SpecwcsModel,
    RegionsModel, WavelengthrangeModel, CameraModel, CollimatorModel, OTEModel,
    FOREModel, FPAModel, IFUPostModel, IFUFOREModel, IFUSlicerModel, MSAModel,
    FilteroffsetModel, DisperserModel, NIRCAMGrismModel, NIRISSGrismModel,
    WaveCorrModel)
from .wfssbkg import WfssBkgModel
from .util import open



__all__ = [
    'open',
    'DataModel',
    'AmiLgModel', 'AsnModel',
    'BarshadowModel', 'CameraModel', 'CollimatorModel',
    'CombinedSpecModel', 'ContrastModel', 'CubeModel',
    'DarkModel', 'DarkMIRIModel',
    'DisperserModel', 'DistortionModel', 'DistortionMRSModel',
    'DrizProductModel',
    'DrizParsModel',
    'Extract1dImageModel',
    'FilteroffsetModel',
    'FlatModel', 'NirspecFlatModel', 'NirspecQuadFlatModel',
    'FOREModel', 'FPAModel',
    'FringeModel', 'GainModel', 'GLS_RampFitModel',
    'GuiderRawModel', 'GuiderCalModel',
    'IFUCubeModel',
    'IFUCubeParsModel', 'NirspecIFUCubeParsModel', 'MiriIFUCubeParsModel',
    'IFUFOREModel', 'IFUImageModel', 'IFUPostModel', 'IFUSlicerModel',
    'ImageModel', 'IPCModel', 'IRS2Model', 'LastFrameModel', 'Level1bModel',
    'LinearityModel', 'MaskModel', 'ModelContainer', 'MSAModel',
    'MultiExposureModel', 'MultiExtract1dImageModel', 'MultiProductModel', 'MultiSlitModel',
    'MultiSpecModel', 'OTEModel',
    'NIRCAMGrismModel','NIRISSGrismModel',
    'OutlierParsModel',
    'PathlossModel',
    'PersistenceSatModel',
    'PixelAreaModel', 'NirspecSlitAreaModel', 'NirspecMosAreaModel', 'NirspecIfuAreaModel',
    'PhotomModel', 'FgsPhotomModel', 'FgsImgPhotomModel', 'MiriImgPhotomModel', 'MirImgPhotomModel', 'MiriMrsPhotomModel',
    'NircamPhotomModel', 'NrcImgPhotomModel', 'NirissPhotomModel', 'NisImgPhotomModel',
    'NirspecPhotomModel', 'NirspecFSPhotomModel',
    'PsfMaskModel',
    'QuadModel', 'RampModel', 'MIRIRampModel',
    'RampFitOutputModel', 'ReadnoiseModel',
    'ReferenceFileModel', 'ReferenceCubeModel', 'ReferenceImageModel', 'ReferenceQuadModel',
    'RegionsModel', 'ResetModel',
    'ResolutionModel', 'MiriResolutionModel',
    'RSCDModel', 'SaturationModel', 'SlitDataModel', 'SlitModel', 'SpecModel',
    'SourceModelContainer', 'StepParsModel',
    'StrayLightModel', 'SuperBiasModel', 'SpecwcsModel',
    'ThroughputModel',
    'TrapDensityModel', 'TrapParsModel', 'TrapsFilledModel',
    'TsoPhotModel',
    'WavelengthrangeModel', 'WaveCorrModel', 'WfssBkgModel']

# Initialize the astropy.io registry,
# but only the first time this module is called

try:
    _defined_models
except NameError:
    with registry.delay_doc_updates(DataModel):
        registry.register_reader('datamodel', DataModel, ndmodel.read)
        registry.register_writer('datamodel', DataModel, ndmodel.write)
        registry.register_identifier('datamodel', DataModel, ndmodel.identify)

_all_models = __all__[1:]
_local_dict = locals()
_defined_models = {k:_local_dict[k] for k in _all_models}
