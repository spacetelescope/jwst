import logging

from jwst.coron import klip
from jwst.stpipe import Step

__all__ = ["KlipStep"]

log = logging.getLogger(__name__)


class KlipStep(Step):
    """
    Performs KLIP processing on a science target coronagraphic exposure.

    The input science exposure is assumed to be a fully calibrated
    level-2b image. The processing is performed using a set of reference PSF
    images observed in the same coronagraphic mode.
    """

    class_alias = "klip"

    spec = """
        truncate = integer(default=50,min=0) # The number of KL transform rows to keep
    """  # noqa: E501

    def process(self, target, psfrefs):
        """
        Execute the KLIP calibration step.

        Parameters
        ----------
        target : str or CubeModel
            CubeModel or file containing science target exposure
        psfrefs : str or CubeModel
            CubeModel or file containing PSF Reference exposures

        Returns
        -------
        psf_sub : CubeModel
            Science target CubeModel with the PSF subtracted
        """
        target_model = self.prepare_output(target)
        refs_model = self.prepare_output(psfrefs)

        # Retrieve the parameter values
        truncate = self.truncate
        log.info("KL transform truncation = %d", truncate)

        # Call the KLIP routine
        target_model = klip.klip(target_model, refs_model, truncate, return_psf=False)

        # Update the step completion status
        target_model.meta.cal_step.klip = "COMPLETE"

        # Return the target model
        return target_model
