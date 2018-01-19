import os.path
from asdf.extension import AsdfExtension
from asdf import util
from .tags import *
from .jwst_types import _jwst_types


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
