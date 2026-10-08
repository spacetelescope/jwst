from jwst.assign_wcs.assign_wcs_step import AssignWcsStep

__all__ = ["get_reference_files"]


def get_reference_files(datamodel):
    """
    Get the reference files associated with the assign_wcs step.

    Parameters
    ----------
    datamodel : `~stdatamodels.jwst.datamodels.JwstDataModel`
        Input datamodel.

    Returns
    -------
    dict
        Keys are reference file type. Values are reference file paths.
    """
    refs = {}
    step = AssignWcsStep()
    for reftype in AssignWcsStep.reference_file_types:
        val = step.get_reference_file(datamodel, reftype)
        if val == "N/A":
            refs[reftype] = None
        else:
            refs[reftype] = val

    return refs
