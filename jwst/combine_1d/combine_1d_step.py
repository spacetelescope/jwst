from stdatamodels.jwst import datamodels

from ..stpipe import Step
from . import combine1d


__all__ = ["Combine1dStep"]


class Combine1dStep(Step):
    """
    Combine 1D spectra.

    Attributes
    ----------
    exptime_key : str
        A case-insensitive string that identifies the metadata element
        (or FITS keyword) for the weight to apply to the input data.  The default
        is "exposure_time".  If the string is "effinttm" or starts with
        "integration", the integration time (FITS keyword EFFINTTM) is used
        as the weight.  If the string is "effexptm" or starts with "exposure",
        the exposure time (FITS keyword EFFEXPTM) is used as the weight.  If
        the string is "unit_weight" or "unit weight", the same weight (1) will
        be used for all input spectra.  If the string is anything else, a warning
        will be logged and unit weight will be used.

    sigma_clip : float or None
        Optional factor for sigma clipping outliers when combining spectra. If
        a floating point value is provided for ``sigma_clip``, this value will be
        used to set an outlier threshold for any pixels in the input spectra that
        deviate from the median and median absolute deviation of the inputs.
        Defaults to None (such that no clipping is performed).
    """

    class_alias = "combine_1d"

    spec = """
    exptime_key = string(default="exposure_time") # Metadata key to use for weighting
    sigma_clip = float(default=None) # Factor for clipping outliers
    """  # noqa: E501

    def process(self, input_data):
        """
        Combine the input data.

        Parameters
        ----------
        input_data : str or ModelContainer or MultiSpecModel
            Input is expected to be an association file name, ModelContainer,
            or MultiSpecModel containing multiple spectra to be combined.
            Individual members of the association or container are expected
            to be MultiSpecModel instances.

        Returns
        -------
        output_spectrum : MultiCombinedSpecModel
            A single combined 1D spectrum.
        """
        with datamodels.open(input_data) as input_model:
            try:
                result = combine1d.combine_1d_spectra(
                    input_model, self.exptime_key, sigma_clip=self.sigma_clip
                )
            except TypeError:
                self.log.error("Invalid input model for combine_1d; skipping.")
                result = input_model.copy()
                result.meta.cal_step.combine_1d = "SKIPPED"

        return result
