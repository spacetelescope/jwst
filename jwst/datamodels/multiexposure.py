from copy import deepcopy

from asdf import schema as asdf_schema
from asdf import treeutil, AsdfFile

from .model_base import JwstDataModel
from .image import ImageModel
from .slit import SlitModel, SlitDataModel


__all__ = ['MultiExposureModel']


class MultiExposureModel(JwstDataModel):
    """
    A data model for multi-slit images derived from
    numerous exposures. The intent is that all slits
    in this model are of the same source, with each slit
    representing a separate exposure of that source.

    This model has a special member `exposures` that can be used to
    deal with an entire slit at a time.  It behaves like a list::

       >>> from .image import ImageModel
       >>> multiexposure_model = MultiExposureModel()
       >>> multiexposure_model.exposures.append(ImageModel())
       >>> multiexposure_model.exposures[0]      # doctest: +SKIP
       <ImageModel>

    Also, there is an extra attribute, `meta`. This will contain the
    meta attribute from the exposure from which each slit has been taken.

    See the module `exp_to_source` for the initial creation of these
    models. This is part of the Level 3 processing of multi-objection
    observations.

    Parameters
    __________
    exposures.items.data : numpy float32 array

    exposures.items.dq : numpy uint32 array

    exposures.items.err : numpy float32 array

    exposures.items.area : numpy float32 array
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/multiexposure.schema"
    core_schema_url = 'http://stsci.edu/schemas/jwst_datamodel/core.schema'

    def __init__(self, init=None, **kwargs):

        # Lets create a schema
        schema = self._build_schema()

        if isinstance(init, (SlitModel, SlitDataModel, ImageModel)):
            super(MultiExposureModel, self).__init__(
                init=None,
                schema=schema,
                **kwargs
            )
            self.update(init)
            self.exposures.append(self.exposures.item())
            self.exposures[0].data = init.data
            self.exposures[0].dq = init.dq
            self.exposures[0].err = init.err
            self.exposures[0].var_poisson = init.var_poisson
            self.exposures[0].var_rnoise = init.var_rnoise
            self.exposures[0].var_flat = init.var_flat
            self.exposures[0].wavelength = init.wavelength
            self.exposures[0].barshadow = init.barshadow
            self.exposures[0].flatfield_point = init.flatfield_point
            self.exposures[0].flatfield_uniform = init.flatfield_uniform
            self.exposures[0].pathloss_point = init.pathloss_point
            self.exposures[0].pathloss_uniform = init.pathloss_uniform
            self.exposures[0].photom_point = init.photom_point
            self.exposures[0].photom_uniform = init.photom_uniform
            self.exposures[0].area = init.area
            return

        super(MultiExposureModel, self).__init__(
            init=init,
            schema=schema,
            **kwargs
        )

    def _build_schema(self):
        """Build the schema, incorporating the core."""
        # Get the schemas
        schema = asdf_schema.load_schema(
            self.schema_url,
            resolver=AsdfFile().resolver,
            resolve_references=True
        )
        core_schema = asdf_schema.load_schema(
            self.core_schema_url,
            resolver=AsdfFile().resolver,
            resolve_references=True
        )

        # Create a new core.meta that will co-locate
        # with each exposure entry. This is done
        # by saving the meta information in a separate
        # FITS HDU.
        core_meta_schema = deepcopy(core_schema['properties']['meta'])
        treeutil.walk(core_meta_schema, remove_fits)
        exposure_schema = schema['allOf'][1]['properties']['exposures']['items']
        exposure_meta_schema = exposure_schema['allOf'][1]['properties']['meta']
        exposure_meta_schema.update(core_meta_schema)

        # That's all folks
        return schema


# Utilities
def set_hdu(obj, hdu_id='EXP'):
    """Add fits_hdu specification to fits-connected properties"""
    try:
        if 'fits_keyword' in obj.keys():
            obj['fits_hdu'] = hdu_id
    except AttributeError:
        pass


def remove_fits(obj):
    try:
        obj.pop('fits_keyword', None)
    except (AttributeError, KeyError, TypeError):
        pass
    try:
        obj.pop('fits_hdu', None)
    except (AttributeError, KeyError, TypeError):
        pass
