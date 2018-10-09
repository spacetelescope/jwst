from .model_base import DataModel


__all__ = ['RampFitOutputModel']


class RampFitOutputModel(DataModel):
    """
    A data model for the optional output of the ramp fitting step.

    In the parameter definitions below, `n_int` is the number of
    integrations, `max_seg` is the maximum number of segments that
    were fit, `nreads` is the number of reads in an integration, and
    `ny` and `nx` are the height and width of the image.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    slope : numpy array (n_int, max_seg, ny, nx)

    sigslope : numpy array (n_int, max_seg, ny, nx)

    var_poisson : numpy array (n_int, max_seg, ny, nx)

    var_rnoise : numpy array (n_int, max_seg, ny, nx)

    yint : numpy array (n_int, max_seg, ny, nx)

    sigyint : numpy array (n_int, max_seg, ny, nx)

    pedestal : numpy array (n_int, max_seg, ny, nx)

    weights : numpy array (n_int, max_seg, ny, nx)

    crmag : numpy array (n_int, max_seg, ny, nx)
    """
    schema_url = "rampfitoutput.schema.yaml"
