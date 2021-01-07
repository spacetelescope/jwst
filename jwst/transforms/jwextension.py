# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os.path
from asdf.extension import AsdfExtension
from asdf import util
from .tags import (GratingEquationType, CoordsType, RotationSequenceType,
                   Gwa2SlitType, Slit2MsaType, LogicalType, NirissSOSSType, V23ToSkyType,
                   RefractionIndexType, SnellType, MIRI_AB2SliceType, NIRCAMGrismDispersionType,
                   NIRISSGrismDispersionType)

from .jwst_types import _jwst_types


__all__ = ['GratingEquationType', 'CoordsType', 'RotationSequenceType',
           'Gwa2SlitType', 'Slit2MsaType', 'LogicalType', 'NirissSOSSType', 'V23ToSkyType',
           'RefractionIndexType', 'SnellType', 'MIRI_AB2SliceType', 'NIRCAMGrismDispersionType',
           'NIRISSGrismDispersionType']


SCHEMA_PATH = os.path.abspath(
    os.path.join(os.path.dirname(__file__), 'schemas'))


class JWSTExtension(AsdfExtension):
    @property
    def types(self):
        return _jwst_types

    @property
    def tag_mapping(self):
        return [('tag:stsci.edu:jwst_pipeline',
                 'http://stsci.edu/schemas/jwst_pipeline{tag_suffix}')]

    @property
    def url_mapping(self):
        return [('http://stsci.edu/schemas/jwst_pipeline/',
                 util.filepath_to_url(os.path.join(SCHEMA_PATH, "stsci.edu")) +
                 '/jwst_pipeline/{url_suffix}.yaml')]
