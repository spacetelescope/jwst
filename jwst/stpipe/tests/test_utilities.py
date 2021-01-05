"""Test utility funcs"""
from jwst.stpipe.utilities import all_steps


# Snapshot of known steps
KNOWN_STEPS = set([
    'AlignRefsStep','Ami3Pipeline', 'AmiAnalyzeStep', 'AmiAverageStep',
    'AmiNormalizeStep', 'AssignMTWcsStep', 'AssignWcsStep',
    'BackgroundStep', 'BarShadowStep',
    'Combine1dStep', 'Coron3Pipeline', 'CubeBuildStep', 'CubeSkyMatchStep',
    'DQInitStep', 'DarkCurrentStep', 'DarkPipeline', 'Detector1Pipeline',
    'Extract1dStep', 'Extract2dStep',
    'FirstFrameStep', 'FlatFieldStep', 'FringeStep',
    'GainScaleStep', 'GroupScaleStep', 'GuiderCdsStep', 'GuiderPipeline',
    'HlspStep',
    'IPCStep', 'Image2Pipeline', 'Image3Pipeline', 'ImprintStep',
    'JumpStep',
    'KlipStep',
    'LastFrameStep', 'LinearityStep',
    'MRSIMatchStep', 'MSAFlagOpenStep', 'MasterBackgroundNrsSlitsStep', 'MasterBackgroundStep',
    'OutlierDetectionScaledStep', 'OutlierDetectionStackStep', 'OutlierDetectionStep',
    'PathLossStep', 'PersistenceStep', 'PhotomStep',
    'RscdStep', 'RSCD_Step', 'RampFitStep', 'RefPixStep', 'ResampleSpecStep', 'ResampleStep', 'ResetStep',
    'SaturationStep', 'SkyMatchStep', 'SourceCatalogStep', 'SourceTypeStep', 'Spec2Pipeline', 'Spec3Pipeline',
    'StackRefsStep', 'StraylightStep', 'SubtractImagesStep', 'SuperBiasStep',
    'TSOPhotometryStep', 'Tso3Pipeline', 'TweakRegStep',
    'WavecorrStep', 'WfsCombineStep', 'WhiteLightStep',
])

def test_all_steps():
    """Test finding all defined steps and pipelines"""
    found_steps = all_steps()
    diff = KNOWN_STEPS.symmetric_difference(found_steps)
    assert not diff, f'Steps not accounted for. Confirm and check suffix and CRDS calpars.\n{diff}'
