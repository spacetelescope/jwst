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
        ("jwst.pipeline.Ami3Pipeline", ['calwebb_ami3', 'jw_ami3'], True),
        ("jwst.pipeline.Coron3Pipeline", ['calwebb_coron3', 'jw_coron3'], True),
        ("jwst.pipeline.DarkPipeline", ['calwebb_dark', 'jw_dark'], True),
        ("jwst.pipeline.Detector1Pipeline", ['calwebb_detector1', 'jw_detector1'], True),
        ("jwst.pipeline.GuiderPipeline", ['calwebb_guider', 'jw_guider'], True),
        ("jwst.pipeline.Image2Pipeline", ['calwebb_image2', 'jw_image2'], True),
        ("jwst.pipeline.Image3Pipeline", ['calwebb_image3', 'jw_image3'], True),
        ("jwst.pipeline.Spec2Pipeline", ['calwebb_spec2', 'jw_spec2'], True),
        ("jwst.pipeline.Spec3Pipeline", ['calwebb_spec3', 'jw_spec3'], True),
        ("jwst.pipeline.Tso3Pipeline", ['calwebb_tso3', 'jw_tso3'], True),
        ("jwst.step.AmiAnalyzeStep", 'jw_ami_analyze', False),
        ("jwst.step.AmiAverageStep", 'jw_ami_average', False),
        ("jwst.step.AmiNormalizeStep", 'jw_ami_normalize', False),
        ("jwst.step.AssignMTWcsStep", 'jw_assign_mtwcs', False),
        ("jwst.step.AssignWcsStep", 'jw_assign_wcs', False),
        ("jwst.step.BackgroundStep", 'jw_background', False),
        ("jwst.step.BarShadowStep", 'jw_barshadow', False),
        ("jwst.step.Combine1dStep", 'jw_combine_1d', False),
        ("jwst.step.StackRefsStep", 'jw_stack_refs', False),
        ("jwst.step.AlignRefsStep", 'jw_align_refs', False),
        ("jwst.step.KlipStep", 'jw_klip', False),
        ("jwst.step.HlspStep", 'jw_hlsp', False),
        ("jwst.step.CubeBuildStep", 'jw_cube_build', False),
        ("jwst.step.CubeSkyMatchStep", 'jw_cube_skymatch', False),
        ("jwst.step.DarkCurrentStep", 'jw_dark_current', False),
        ("jwst.step.DQInitStep", 'jw_dq_init', False),
        ("jwst.step.Extract1dStep", 'jw_extract_1d', False),
        ("jwst.step.Extract2dStep", 'jw_extract_2d', False),
        ("jwst.step.FirstFrameStep", 'jw_firstframe', False),
        ("jwst.step.FlatFieldStep", 'jw_flat_field', False),
        ("jwst.step.FringeStep", 'jw_fringe', False),
        ("jwst.step.GainScaleStep", 'jw_gain_scale', False),
        ("jwst.step.GroupScaleStep", 'jw_group_scale', False),
        ("jwst.step.GuiderCdsStep", 'jw_guider_cds', False),  #GuiderCDS name?
        ("jwst.step.ImprintStep", 'jw_imprint', False),
        ("jwst.step.IPCStep", 'jw_ipc', False),
        ("jwst.step.JumpStep", 'jw_jump', False),
        ("jwst.step.LastFrameStep", 'jw_lastframe', False),
        ("jwst.step.LinearityStep", 'jw_linearity', False),
        ("jwst.step.MasterBackgroundStep", 'jw_master_background', False),
        ("jwst.step.MasterBackgroundNrsSlitsStep", 'jw_master_background_nrs', False),
        ("jwst.step.MRSIMatchStep", 'jw_mrs_imatch', False),
        ("jwst.step.MSAFlagOpenStep", 'jw_msa_flagging', False),
        ("jwst.step.OutlierDetectionStep", 'jw_outlier_detection', False),
        ("jwst.step.OutlierDetectionScaledStep", 'jw_outlier_detection_scaled', False),
        ("jwst.step.OutlierDetectionStackStep", 'jw_outlier_detection_stack', False),
        ("jwst.step.PathLossStep", 'jw_pathloss', False),
        ("jwst.step.PersistenceStep", 'jw_persistence', False),
        ("jwst.step.PhotomStep", 'jw_photom', False),
        ("jwst.step.RampFitStep", 'jw_ramp_fit', False),  #RampFit name?
        ("jwst.step.RefPixStep", 'jw_refpix', False),
        ("jwst.step.ResampleStep", 'jw_resample', False),
        ("jwst.step.ResampleSpecStep", 'jw_resample_spec', False),
        ("jwst.step.ResetStep", 'jw_reset', False),
        ("jwst.step.RscdStep", 'jw_rscd', False),
        ("jwst.step.RSCD_Step", None, False),
        ("jwst.step.SaturationStep", 'jw_saturation', False),
        ("jwst.step.SkyMatchStep", 'jw_skymatch', False),
        ("jwst.step.SourceCatalogStep", 'jw_source_catalog', False),
        ("jwst.step.SourceTypeStep", 'jw_srctype', False),
        ("jwst.step.StraylightStep", 'jw_straylight', False),
        ("jwst.step.SuperBiasStep", 'jw_superbias', False),
        ("jwst.step.TSOPhotometryStep", 'jw_tso_photometry', False),
        ("jwst.step.TweakRegStep", 'jw_tweakreg', False),
        ("jwst.step.WavecorrStep", 'jw_wavecorr', False),
        ("jwst.step.WfsCombineStep", ['calwebb_wfs-image3', 'jw_wfs-image3'], False),  #WfsCombine name?
        ("jwst.step.WfssContamStep", 'jw_wfss_contam', False),
        ("jwst.step.WhiteLightStep", 'jw_white_light', False),
    ]
