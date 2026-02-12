from jwst.extract_2d import extract_2d
from jwst.stpipe import Step

__all__ = ["Extract2dStep"]


class Extract2dStep(Step):
    """Extract 2D spectral cutouts."""

    class_alias = "extract_2d"

    spec = """
        slit_names = force_list(default=None)   # slits to be extracted
        source_ids = force_list(default=None)     # source ids to be extracted
        source_ra = force_list(default=None)  # source RAs to be extracted, WFSS only
        source_dec = force_list(default=None)  # source DECs to be extracted, WFSS only
        source_max_sep = float(default=2.0)  # maximum separation in arcseconds around source_ra, source_dec, WFSS only
        extract_orders = int_list(default=None)  # list of orders to extract
        grism_objects = list(default=None)  # list of grism objects to use
        tsgrism_extract_height =  integer(default=None)  # extraction height in pixels, TSGRISM mode
        wfss_extract_half_height =  integer(default=5)  # extraction half height in pixels, WFSS mode
        wfss_mmag_extract = float(default=None)  # minimum abmag to extract, WFSS mode
        wfss_nbright = integer(default=1000)  # number of brightest objects to extract, WFSS mode
    """  # noqa: E501

    reference_file_types = ["wavelengthrange"]

    def process(self, input_data):
        """
        Extract 2D cutouts from a spectral image.

        Parameters
        ----------
        input_data : str, `~stdatamodels.jwst.datamodels.ImageModel`, or \
                     `~stdatamodels.jwst.datamodels.CubeModel`
            Input datamodel or file name.

        Returns
        -------
        result : `~stdatamodels.jwst.datamodels.MultiSlitModel`, or \
                 `~stdatamodels.jwst.datamodels.SlitModel`
            Output datamodel containing spectral cutouts in separate extensions.
        """
        output_model = self.prepare_output(input_data)

        reference_file_names = {}
        if output_model.meta.exposure.type in extract_2d.slitless_modes:
            # The wavelengthrange file is used only by the WFSS modes.
            # If retrieved by a Nirspec mode, it would override the name of
            # the file in meta.ref_file if a custom file was used.
            for reftype in self.reference_file_types:
                reffile = self.get_reference_file(output_model, reftype)
                reference_file_names[reftype] = reffile if reffile else ""

        result = extract_2d.extract2d(
            output_model,
            self.slit_names,
            self.source_ids,
            reference_files=reference_file_names,
            grism_objects=self.grism_objects,
            tsgrism_extract_height=self.tsgrism_extract_height,
            wfss_extract_half_height=self.wfss_extract_half_height,
            extract_orders=self.extract_orders,
            mmag_extract=self.wfss_mmag_extract,
            nbright=self.wfss_nbright,
            source_ra=self.source_ra,
            source_dec=self.source_dec,
            max_sep=self.source_max_sep,
        )

        # Result is a new model if the step succeeded,
        # so close the input in that case if it was opened here
        if output_model is not input_data and output_model is not result:
            output_model.close()

        return result
