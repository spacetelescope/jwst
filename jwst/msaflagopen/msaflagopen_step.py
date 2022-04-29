from ..stpipe import Step
from .. import datamodels
from . import msaflag_open

__all__ = ["MSAFlagOpenStep"]


class MSAFlagOpenStep(Step):
    """
    MSAFlagOpenStep: Flags pixels affected by MSA failed open shutters
    """

    class_alias = "msa_flagging"

    spec = """

    """

    reference_file_types = ['msaoper']

    def process(self, input):

        # Open the input data model
        with datamodels.open(input) as input_model:
            self.reference_name = self.get_reference_file(input_model,
                                                          'msaoper')
            self.log.info('Using reference file %s', self.reference_name)

            # Check for a valid reference file
            if self.reference_name == 'N/A':
                self.log.warning('No reference file found')
                self.log.warning('Step will be skipped')
                result = input_model.copy()
                result.meta.cal_step.msa_flagging = 'SKIPPED'
                return result

            #
            # Get the reference file names for constructing the WCS pipeline
            wcs_reffile_names = create_reference_filename_dictionary(input_model)
            # Do the DQ flagging
            result = msaflag_open.do_correction(input_model,
                                                self.reference_name,
                                                wcs_reffile_names)

            # set the step status to complete
            result.meta.cal_step.msa_flagging = 'COMPLETE'

        return result


def create_reference_filename_dictionary(input_model):
    reffiles = {}
    a = Step()
    reffiles['distortion'] = input_model.meta.ref_file.distortion.name
    reffiles['filteroffset'] = input_model.meta.ref_file.filteroffset.name
    reffiles['specwcs'] = input_model.meta.ref_file.filteroffset.name
    reffiles['regions'] = input_model.meta.ref_file.regions.name
    reffiles['wavelengthrange'] = input_model.meta.ref_file.wavelengthrange.name
    if input_model.meta.ref_file.v2v3.name is not None:
        reffiles['v2v3'] = input_model.meta.ref_file.v2v3.name
    reffiles['camera'] = input_model.meta.ref_file.camera.name
    reffiles['collimator'] = input_model.meta.ref_file.collimator.name
    reffiles['disperser'] = input_model.meta.ref_file.disperser.name
    reffiles['fore'] = input_model.meta.ref_file.fore.name
    reffiles['fpa'] = input_model.meta.ref_file.fpa.name
    reffiles['msa'] = input_model.meta.ref_file.msa.name
    reffiles['ote'] = input_model.meta.ref_file.ote.name
    reffiles['ifupost'] = input_model.meta.ref_file.ifupost.name
    reffiles['ifufore'] = input_model.meta.ref_file.ifufore.name
    reffiles['ifuslicer'] = input_model.meta.ref_file.ifuslicer.name
    # Convert from crds protocol to absolute filenames
    for key in reffiles.keys():
        if reffiles[key].startswith('crds://'):
            reffiles[key] = a.reference_uri_to_cache_path(reffiles[key], input_model.crds_observatory)
    return reffiles
