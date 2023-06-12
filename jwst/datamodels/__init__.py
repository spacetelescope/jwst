# The jwst.datamodels submodule was moved to stdatamodels.jwst.datamodels
# https://github.com/spacetelescope/jwst/pull/7439

import importlib
from inspect import ismodule
import sys

from stdatamodels.jwst.datamodels.model_base import JwstDataModel, DataModel
from stdatamodels.jwst.datamodels.abvega_offset import ABVegaOffsetModel
from stdatamodels.jwst.datamodels.amilg import AmiLgModel
from stdatamodels.jwst.datamodels.apcorr import FgsImgApcorrModel, MirImgApcorrModel
from stdatamodels.jwst.datamodels.apcorr import NrcImgApcorrModel, NisImgApcorrModel
from stdatamodels.jwst.datamodels.apcorr import MirLrsApcorrModel, MirMrsApcorrModel
from stdatamodels.jwst.datamodels.apcorr import NrcWfssApcorrModel, NisWfssApcorrModel
from stdatamodels.jwst.datamodels.apcorr import NrsMosApcorrModel, NrsIfuApcorrModel, NrsFsApcorrModel
from stdatamodels.jwst.datamodels.asn import AsnModel
from stdatamodels.jwst.datamodels.barshadow import BarshadowModel
from stdatamodels.jwst.datamodels.combinedspec import CombinedSpecModel
from stdatamodels.jwst.datamodels.contrast import ContrastModel
from stdatamodels.jwst.datamodels.cube import CubeModel
from stdatamodels.jwst.datamodels.dark import DarkModel
from stdatamodels.jwst.datamodels.darkMIRI import DarkMIRIModel
from stdatamodels.jwst.datamodels.drizpars import DrizParsModel
from stdatamodels.jwst.datamodels.extract1dimage import Extract1dImageModel
from stdatamodels.jwst.datamodels.extract1d_spec import Extract1dIFUModel
from stdatamodels.jwst.datamodels.flat import FlatModel
from stdatamodels.jwst.datamodels.fringe import FringeModel
from stdatamodels.jwst.datamodels.fringefreq import FringeFreqModel
from stdatamodels.jwst.datamodels.gain import GainModel
from stdatamodels.jwst.datamodels.gls_rampfit import GLS_RampFitModel
from stdatamodels.jwst.datamodels.guider import GuiderRawModel, GuiderCalModel
from stdatamodels.jwst.datamodels.ifucube import IFUCubeModel
from stdatamodels.jwst.datamodels.ifucubepars import NirspecIFUCubeParsModel, MiriIFUCubeParsModel
from stdatamodels.jwst.datamodels.ifuimage import IFUImageModel
from stdatamodels.jwst.datamodels.mrsptcorr import MirMrsPtCorrModel
from stdatamodels.jwst.datamodels.image import ImageModel
from stdatamodels.jwst.datamodels.ipc import IPCModel
from stdatamodels.jwst.datamodels.irs2 import IRS2Model
from stdatamodels.jwst.datamodels.lastframe import LastFrameModel
from stdatamodels.jwst.datamodels.level1b import Level1bModel
from stdatamodels.jwst.datamodels.linearity import LinearityModel
from stdatamodels.jwst.datamodels.mask import MaskModel
from stdatamodels.jwst.datamodels.mrsxartcorr import MirMrsXArtCorrModel
from stdatamodels.jwst.datamodels.multicombinedspec import MultiCombinedSpecModel
from stdatamodels.jwst.datamodels.multiexposure import MultiExposureModel
from stdatamodels.jwst.datamodels.multiextract1d import MultiExtract1dImageModel
from stdatamodels.jwst.datamodels.multislit import MultiSlitModel
from stdatamodels.jwst.datamodels.multispec import MultiSpecModel
from stdatamodels.jwst.datamodels.nirspec_flat import NirspecFlatModel, NirspecQuadFlatModel
from stdatamodels.jwst.datamodels.outlierpars import OutlierParsModel
from stdatamodels.jwst.datamodels.outlierifuoutput import OutlierIFUOutputModel
from stdatamodels.jwst.datamodels.pathloss import PathlossModel, MirLrsPathlossModel
from stdatamodels.jwst.datamodels.persat import PersistenceSatModel
from stdatamodels.jwst.datamodels.photom import FgsImgPhotomModel
from stdatamodels.jwst.datamodels.photom import MirImgPhotomModel, MirLrsPhotomModel, MirMrsPhotomModel
from stdatamodels.jwst.datamodels.photom import NrcImgPhotomModel, NrcWfssPhotomModel
from stdatamodels.jwst.datamodels.photom import NisImgPhotomModel, NisSossPhotomModel, NisWfssPhotomModel
from stdatamodels.jwst.datamodels.photom import NrsFsPhotomModel, NrsMosPhotomModel
from stdatamodels.jwst.datamodels.pixelarea import PixelAreaModel, NirspecSlitAreaModel, NirspecMosAreaModel, NirspecIfuAreaModel
from stdatamodels.jwst.datamodels.psfmask import PsfMaskModel
from stdatamodels.jwst.datamodels.quad import QuadModel
from stdatamodels.jwst.datamodels.ramp import RampModel
from stdatamodels.jwst.datamodels.rampfitoutput import RampFitOutputModel
from stdatamodels.jwst.datamodels.readnoise import ReadnoiseModel
from stdatamodels.jwst.datamodels.reference import ReferenceFileModel, ReferenceImageModel, ReferenceCubeModel, ReferenceQuadModel
from stdatamodels.jwst.datamodels.reset import ResetModel
from stdatamodels.jwst.datamodels.resolution import ResolutionModel, MiriResolutionModel
from stdatamodels.jwst.datamodels.rscd import RSCDModel
from stdatamodels.jwst.datamodels.saturation import SaturationModel
from stdatamodels.jwst.datamodels.segmap import SegmentationMapModel
from stdatamodels.jwst.datamodels.slit import SlitModel, SlitDataModel
from stdatamodels.jwst.datamodels.sossextractmodel import SossExtractModel
from stdatamodels.jwst.datamodels.sosswavegrid import SossWaveGridModel
from stdatamodels.jwst.datamodels.spec import SpecModel
from stdatamodels.jwst.datamodels.speckernel import SpecKernelModel
from stdatamodels.jwst.datamodels.specprofile import SpecProfileModel, SpecProfileSingleModel
from stdatamodels.jwst.datamodels.spectrace import SpecTraceModel, SpecTraceSingleModel
from stdatamodels.jwst.datamodels.straylight import StrayLightModel
from stdatamodels.jwst.datamodels.superbias import SuperBiasModel
from stdatamodels.jwst.datamodels.throughput import ThroughputModel
from stdatamodels.jwst.datamodels.trapdensity import TrapDensityModel
from stdatamodels.jwst.datamodels.trappars import TrapParsModel
from stdatamodels.jwst.datamodels.trapsfilled import TrapsFilledModel
from stdatamodels.jwst.datamodels.tsophot import TsoPhotModel
from stdatamodels.jwst.datamodels.wavemap import WaveMapModel, WaveMapSingleModel
from stdatamodels.jwst.datamodels.wcs_ref_models import (DistortionModel, DistortionMRSModel, SpecwcsModel,
                             RegionsModel, WavelengthrangeModel, CameraModel, CollimatorModel, OTEModel,
                             FOREModel, FPAModel, IFUPostModel, IFUFOREModel, IFUSlicerModel, MSAModel,
                             FilteroffsetModel, DisperserModel, NIRCAMGrismModel, NIRISSGrismModel,
                             WaveCorrModel)
from stdatamodels.jwst.datamodels.wfssbkg import WfssBkgModel
>>>>>>> b709e8466 (updates from review)
from stdatamodels.jwst.datamodels.util import open

from .container import ModelContainer
from .source_container import SourceModelContainer

import stdatamodels.jwst.datamodels

# Import everything defined in stdatamodels.jwst.datamodels.__all__
from stdatamodels.jwst.datamodels import * # noqa: F403

# Define __all__ to include stdatamodels.jwst.datamodels.__all__
__all__ = [
    'open',
    'ModelContainer', 'SourceModelContainer',
] + stdatamodels.jwst.datamodels.__all__


# Modules that are not part of stdatamodels
_jwst_modules = ["container", "source_container"]
<<<<<<< HEAD

# Models that are not part of stdatamodels
_jwst_models = ["ModelContainer", "SourceModelContainer"]

# Deprecated modules in stdatamodels
_deprecated_modules = ['drizproduct', 'multiprod']

# Deprecated models in stdatamodels
_deprecated_models = ['DrizProductModel', 'MultiProductModel', 'MIRIRampModel']

# Import all submodules from stdatamodels.jwst.datamodels
for attr in dir(stdatamodels.jwst.datamodels):
    if attr[0] == '_':
        continue
    if attr in _jwst_models or attr in _deprecated_modules or attr in _deprecated_models:
        continue
    obj = getattr(stdatamodels.jwst.datamodels, attr)
    if ismodule(obj):
        # Make the submodule available locally
        locals()[attr] = obj
        # Add the submodule to sys.modules so that a call
        # to jwst.datamodels.dqflags will return the submodule
        # stdatamodels.jwst.datamodels.dqflags
        sys.modules[f"jwst.datamodels.{attr}"] = obj

# Add a few submodules to sys.modules without exposing them locally
for _submodule_name in ['schema', 'schema_editor', 'validate']:
    _submodule = importlib.import_module(f"stdatamodels.jwst.datamodels.{_submodule_name}")
    sys.modules[f"jwst.datamodels.{_submodule_name}"] = _submodule
=======
>>>>>>> b709e8466 (updates from review)
