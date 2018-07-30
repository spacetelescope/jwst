from ..stpipe import Step
from .. import datamodels
from . import path_loss

__all__ = ["PathLossStep"]


class PathLossStep(Step):
    """
    PathLossStep: Inserts the pathloss and wavelength arrays
    into the data.

    Pathloss depends on the centering of the source in the aperture if the
    source is a point source. This step fills the following attributes in the
    datamodel:

    - for exposure type `NRS_IFU`, the 1-d arrays .wavelength_pointsource,
      .pathloss_pointsource, .wavelength_uniformsource and
      .pathloss_uniformsource

    - for exposure types NRS_FIXEDSLIT, NRS_BRIGHTOBJ, and NRS_MSASPEC, the
      1-d arrays .slits[n].wavelength_pointsource,
      .slits[n].pathloss_pointsource, .slits[n].wavelength_uniformsource and
      .slits[n].pathloss_uniformsource

    In all of these `EXP_TYPES`, these arrays are added to each
    member of the `slits` attribute.
    """

    spec = """
    """

    reference_file_types = ['pathloss']

    def process(self, input):

        # Open the input data model
        with datamodels.open(input) as input_model:

            # Get the name of the path loss reference file to use
            self.pathloss_name = self.get_reference_file(input_model,
                                                         'pathloss')
            self.log.info('Using PATHLOSS reference file %s',
                          self.pathloss_name)

            # Check for a valid reference file
            if self.pathloss_name == 'N/A':
                self.log.warning('No PATHLOSS reference file found')
                self.log.warning('Path loss step will be skipped')
                result = input_model.copy()
                result.meta.cal_step.pathloss = 'SKIPPED'
                return result

            instrument = input_model.meta.instrument.name
            # Open the pathloss ref file data model
            pathloss_model = datamodels.PathlossModel(self.pathloss_name)

            # Do the pathloss correction
            result = path_loss.do_correction(input_model, pathloss_model)

            pathloss_model.close()

        return result
