from astropy.io import registry

from stdatamodels import ndmodel

from .model_base import JwstDataModel, DataModel
from .abvega_offset import ABVegaOffsetModel
from .amilg import AmiLgModel
from .apcorr import FgsImgApcorrModel, MirImgApcorrModel
from .apcorr import NrcImgApcorrModel, NisImgApcorrModel
from .apcorr import MirLrsApcorrModel, MirMrsApcorrModel
from .apcorr import NrcWfssApcorrModel, NisWfssApcorrModel
from .apcorr import NrsMosApcorrModel, NrsIfuApcorrModel, NrsFsApcorrModel
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
from .extract1d_spec import Extract1dIFUModel
from .flat import FlatModel
from .fringe import FringeModel
from .gain import GainModel
from .gls_rampfit import GLS_RampFitModel
from .guider import GuiderRawModel, GuiderCalModel
from .ifucube import IFUCubeModel
from .ifucubepars import NirspecIFUCubeParsModel, MiriIFUCubeParsModel
from .ifuimage import IFUImageModel
from .image import ImageModel
from .ipc import IPCModel
from .irs2 import IRS2Model
from .lastframe import LastFrameModel
from .level1b import Level1bModel
from .linearity import LinearityModel
from .mask import MaskModel
from .ramp import MIRIRampModel
from .multicombinedspec import MultiCombinedSpecModel
from .multiexposure import MultiExposureModel
from .multiextract1d import MultiExtract1dImageModel
from .multiprod import MultiProductModel
from .multislit import MultiSlitModel
from .multispec import MultiSpecModel
from .nirspec_flat import NirspecFlatModel, NirspecQuadFlatModel
from .outlierpars import OutlierParsModel
from .pathloss import PathlossModel
from .persat import PersistenceSatModel
from .photom import FgsImgPhotomModel
from .photom import MirImgPhotomModel, MirLrsPhotomModel, MirMrsPhotomModel
from .photom import NrcImgPhotomModel, NrcWfssPhotomModel
from .photom import NisImgPhotomModel, NisSossPhotomModel, NisWfssPhotomModel
from .photom import NrsFsPhotomModel, NrsMosPhotomModel
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
    'DataModel', 'JwstDataModel',
    'ABVegaOffsetModel',
    'AmiLgModel',
    'FgsImgApcorrModel', 'MirImgApcorrModel', 'NrcImgApcorrModel', 'NisImgApcorrModel',
    'MirLrsApcorrModel', 'MirMrsApcorrModel', 'NrcWfssApcorrModel', 'NisWfssApcorrModel',
    'NrsMosApcorrModel', 'NrsFsApcorrModel','NrsIfuApcorrModel',
    'AsnModel',
    'BarshadowModel', 'CameraModel', 'CollimatorModel',
    'CombinedSpecModel', 'ContrastModel', 'CubeModel',
    'DarkModel', 'DarkMIRIModel',
    'DisperserModel', 'DistortionModel', 'DistortionMRSModel',
    'DrizParsModel',
    'DrizProductModel',
    'Extract1dImageModel',
    'Extract1dIFUModel',
    'FilteroffsetModel',
    'FlatModel', 'NirspecFlatModel', 'NirspecQuadFlatModel',
    'FOREModel', 'FPAModel',
    'FringeModel', 'GainModel', 'GLS_RampFitModel',
    'GuiderRawModel', 'GuiderCalModel',
    'IFUCubeModel',
    'NirspecIFUCubeParsModel', 'MiriIFUCubeParsModel',
    'IFUFOREModel', 'IFUImageModel', 'IFUPostModel', 'IFUSlicerModel',
    'ImageModel', 'IPCModel', 'IRS2Model', 'LastFrameModel', 'Level1bModel',
    'LinearityModel', 'MaskModel', 'ModelContainer', 'MSAModel',
    'MultiCombinedSpecModel','MultiExposureModel',
    'MultiExtract1dImageModel', 'MultiSlitModel',
    'MultiProductModel',
    'MultiSpecModel', 'OTEModel',
    'NIRCAMGrismModel','NIRISSGrismModel',
    'OutlierParsModel',
    'PathlossModel',
    'PersistenceSatModel',
    'PixelAreaModel', 'NirspecSlitAreaModel', 'NirspecMosAreaModel', 'NirspecIfuAreaModel',
    'FgsImgPhotomModel',
    'MirImgPhotomModel', 'MirLrsPhotomModel', 'MirMrsPhotomModel',
    'NrcImgPhotomModel', 'NrcWfssPhotomModel',
    'NisImgPhotomModel', 'NisSossPhotomModel', 'NisWfssPhotomModel',
    'NrsFsPhotomModel', 'NrsMosPhotomModel',
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
    with registry.delay_doc_updates(JwstDataModel):
        registry.register_reader('datamodel', JwstDataModel, ndmodel.read)
        registry.register_writer('datamodel', JwstDataModel, ndmodel.write)
        registry.register_identifier('datamodel', JwstDataModel, ndmodel.identify)

_all_models = __all__[1:]
_local_dict = locals()
_defined_models = {k:_local_dict[k] for k in _all_models}
