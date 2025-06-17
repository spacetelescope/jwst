#! /usr/bin/env python
from stdatamodels.jwst import datamodels

from jwst.stpipe import Step
from . import extract_2d


__all__ = ["Extract2dStep"]


class Extract2dStep(Step):
    """Class that provides method to perform a 2D extraction of spectra."""

    class_alias = "extract_2d"

    spec = """
        slit_names = force_list(default=None)   # slits to be extracted
        source_ids = force_list(default=None)     # source ids to be extracted
        extract_orders = int_list(default=None)  # list of orders to extract
        grism_objects = list(default=None)  # list of grism objects to use
        tsgrism_extract_height =  integer(default=None)  # extraction height in pixels, TSGRISM mode
        wfss_extract_half_height =  integer(default=5)  # extraction half height in pixels, WFSS mode
        wfss_mmag_extract = float(default=None)  # minimum abmag to extract, WFSS mode
        wfss_nbright = integer(default=1000)  # number of brightest objects to extract, WFSS mode
    """  # noqa: E501

    reference_file_types = ["wavelengthrange"]

    def process(self, input_model):
        """
        Perform the extract_2d calibration step.

        Parameters
        ----------
        input_model : DataModel
            DataModel to be processed

        Returns
        -------
        output_model : DataModel
            The resulting DataModel of the extract_2d step
        """
        reference_file_names = {}
        if input_model.meta.exposure.type in extract_2d.slitless_modes:
            # The wavelengthrange file is used only by the WFSS modes.
            # If retrieved by a Nirspec mode, it would override the name of
            # the file in meta.ref_file if a custom file was used.
            for reftype in self.reference_file_types:
                reffile = self.get_reference_file(input_model, reftype)
                reference_file_names[reftype] = reffile if reffile else ""
        with datamodels.open(input_model) as dm:
            output_model = extract_2d.extract2d(
                dm,
                self.slit_names,
                self.source_ids,
                reference_files=reference_file_names,
                grism_objects=self.grism_objects,
                tsgrism_extract_height=self.tsgrism_extract_height,
                wfss_extract_half_height=self.wfss_extract_half_height,
                extract_orders=self.extract_orders,
                mmag_extract=self.wfss_mmag_extract,
                nbright=self.wfss_nbright,
            )

        return output_model
