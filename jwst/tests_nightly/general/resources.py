from ci_watson.base_classes import BaseTest

class BaseJWSTTest(BaseTest):
    ignore_hdus = ['ASDF']
    ignore_keywords = ['DATE', 'CAL_VER', 'CAL_VCS', 'CRDS_VER', 'CRDS_CTX']
    input_repo = 'scsb-jwst-pipeline'
    copy_local = False  # Do not make additional copy by default
    rtol = 0.00001

    def set_environ(self):
        self.tree = ''

class NIRCamTest(BaseJWSTTest):
    input_loc = 'nircam'

class MIRITest(BaseJWSTTest):
    input_loc = 'miri'

class NIRISSTest(BaseJWSTTest):
    input_loc = 'niriss'

class NIRSpecTest(BaseJWSTTest):
    input_loc = 'nirspec'

class FGSTest(BaseJWSTTest):
    input_loc = 'fgs'
