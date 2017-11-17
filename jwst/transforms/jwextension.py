from __future__ import absolute_import, division, unicode_literals, print_function
import os.path
from asdf.extension import AsdfExtension
from asdf import util
from .tags import *

SCHEMA_PATH = os.path.abspath(
    os.path.join(os.path.dirname(__file__), 'schemas'))


class JWSTExtension(AsdfExtension):
    @property
    def types(self):
        return [GratingEquationType,
                CoordsType,
                RotationSequenceType,
                LRSWavelengthType,
                Gwa2SlitType,
                Slit2MsaType,
                MIRI_AB2SliceType,
                SnellType,
                NIRCAMGrismDispersionType,
                NIRISSGrismDispersionType,
                LogicalType,
                TPCorrType,
                ]

    @property
    def tag_mapping(self):
        # return [('tag:stsci.edu:jwst',
                 # 'http://stsci.edu/schemas/jwst{tag_suffix}')]
        return [('tag:stsci.edu:jwst_pipeline',
                 'http://stsci.edu/schemas/jwst_pipeline{tag_suffix}')]

    @property
    def url_mapping(self):
        return [('http://stsci.edu/schemas/jwst_pipeline/',
                 util.filepath_to_url(SCHEMA_PATH) +
                 '/{url_suffix}.yaml')]
