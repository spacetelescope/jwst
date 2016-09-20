from __future__ import absolute_import, division, unicode_literals, print_function
import os.path
from asdf.extension import AsdfExtension
from asdf import util

SCHEMA_PATH = os.path.abspath(
    os.path.join(os.path.dirname(__file__), 'schemas'))


class BaseExtension(AsdfExtension):
    @property
    def types(self):
        return []

    @property
    def tag_mapping(self):
        return []

    @property
    def url_mapping(self):
        return [('http://jwst.stsci.edu/schemas/',
                 util.filepath_to_url(SCHEMA_PATH) +
                 '/{url_suffix}')]
