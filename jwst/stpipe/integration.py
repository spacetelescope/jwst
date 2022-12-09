"""
Entry point implementations.
"""


def get_steps():
    """
    Return tuples describing the stpipe.Step subclasses provided
    by this package.  This method is registered with the stpipe.steps
    entry point.

    Returns
    -------
    list of tuple (str, str, bool)
        The first element each tuple is a fully-qualified Step
        subclass name.  The second element is an optional class
        alias.  The third element indicates that the class
        is a subclass of Pipeline.
    """
    # Unit tests ensure that this list is kept in sync with the actual
    # class definitions.  We need to avoid importing jwst.pipeline and
    # jwst.step to keep the CLI snappy.
    return [
        ("jwst.pipeline.Ami3Pipeline", 'calwebb_ami3', True),
        ("jwst.pipeline.Coron3Pipeline", 'calwebb_coron3', True),
        ("jwst.pipeline.DarkPipeline", 'calwebb_dark', True),
        ("jwst.pipeline.Detector1Pipeline", 'calwebb_detector1', True),
        ("jwst.pipeline.GuiderPipeline", 'calwebb_guider', True),
        ("jwst.pipeline.Image2Pipeline", 'calwebb_image2', True),
        ("jwst.pipeline.Image3Pipeline", 'calwebb_image3', True),
        ("jwst.pipeline.Spec2Pipeline", 'calwebb_spec2', True),
        ("jwst.pipeline.Spec3Pipeline", 'calwebb_spec3', True),
        ("jwst.pipeline.Tso3Pipeline", 'calwebb_tso3', True),
        ("jwst.step.AmiAnalyzeStep", 'ami_analyze', False),
        ("jwst.step.AmiAverageStep", 'ami_average', False),
        ("jwst.step.AmiNormalizeStep", 'ami_normalize', False),
        ("jwst.step.AssignMTWcsStep", 'assign_mtwcs', False),
        ("jwst.step.AssignWcsStep", 'assign_wcs', False),
        ("jwst.step.BackgroundStep", 'background', False),
        ("jwst.step.BarShadowStep", 'barshadow', False),
        ("jwst.step.Combine1dStep", 'combine_1d', False),
        ("jwst.step.StackRefsStep", 'stack_refs', False),
        ("jwst.step.AlignRefsStep", 'align_refs', False),
        ("jwst.step.KlipStep", 'klip', False),
        ("jwst.step.HlspStep", 'hlsp', False),
        ("jwst.step.CubeBuildStep", 'cube_build', False),
        ("jwst.step.CubeSkyMatchStep", 'cube_skymatch', False),
        ("jwst.step.DarkCurrentStep", 'dark_current', False),
        ("jwst.step.DQInitStep", 'dq_init', False),
        ("jwst.step.Extract1dStep", 'extract_1d', False),
        ("jwst.step.Extract2dStep", 'extract_2d', False),
        ("jwst.step.FirstFrameStep", 'firstframe', False),
        ("jwst.step.FlatFieldStep", 'flat_field', False),
        ("jwst.step.FringeStep", 'fringe', False),
        ("jwst.step.GainScaleStep", 'gain_scale', False),
        ("jwst.step.GroupScaleStep", 'group_scale', False),
        ("jwst.step.GuiderCdsStep", 'guider_cds', False),
        ("jwst.step.ImprintStep", 'imprint', False),
        ("jwst.step.IPCStep", 'ipc', False),
        ("jwst.step.JumpStep", 'jump', False),
        ("jwst.step.LastFrameStep", 'lastframe', False),
        ("jwst.step.LinearityStep", 'linearity', False),
        ("jwst.step.MasterBackgroundStep", 'master_background', False),
        ("jwst.step.MasterBackgroundMosStep", 'master_background_mos', False),
        ("jwst.step.MRSIMatchStep", 'mrs_imatch', False),
        ("jwst.step.MSAFlagOpenStep", 'msa_flagging', False),
        ("jwst.step.OutlierDetectionStep", 'outlier_detection', False),
        ("jwst.step.OutlierDetectionScaledStep", 'outlier_detection_scaled', False),
        ("jwst.step.OutlierDetectionStackStep", 'outlier_detection_stack', False),
        ("jwst.step.PathLossStep", 'pathloss', False),
        ("jwst.step.PersistenceStep", 'persistence', False),
        ("jwst.step.PhotomStep", 'photom', False),
        ("jwst.step.PixelReplaceStep", 'pixel_replace', False),
        ("jwst.step.RampFitStep", 'ramp_fit', False),
        ("jwst.step.RefPixStep", 'refpix', False),
        ("jwst.step.ResampleStep", 'resample', False),
        ("jwst.step.ResampleSpecStep", 'resample_spec', False),
        ("jwst.step.ResetStep", 'reset', False),
        ("jwst.step.ResidualFringeStep", 'residual_fringe', False),
        ("jwst.step.RscdStep", 'rscd', False),
        ("jwst.step.SaturationStep", 'saturation', False),
        ("jwst.step.SkyMatchStep", 'skymatch', False),
        ("jwst.step.SourceCatalogStep", 'source_catalog', False),
        ("jwst.step.SourceTypeStep", 'srctype', False),
        ("jwst.step.StraylightStep", 'straylight', False),
        ("jwst.step.SuperBiasStep", 'superbias', False),
        ("jwst.step.TSOPhotometryStep", 'tso_photometry', False),
        ("jwst.step.TweakRegStep", 'tweakreg', False),
        ("jwst.step.UndersamplingCorrectionStep", 'undersampling_correction', False),
        ("jwst.step.WavecorrStep", 'wavecorr', False),
        ("jwst.step.WfsCombineStep", 'calwebb_wfs-image3', False),
        ("jwst.step.WfssContamStep", 'wfss_contam', False),
        ("jwst.step.WhiteLightStep", 'white_light', False),
    ]
