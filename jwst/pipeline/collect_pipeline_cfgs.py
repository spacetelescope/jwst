# Collect the pipeline configurations

import os
import shutil
from glob import glob
from importlib.util import find_spec


def collect_pipeline_cfgs(dst='./'):
    """Copy step and pipeline .cfg files to destination"""
    os.makedirs(dst, exist_ok=True)

    cfg_dir = os.path.join(find_spec('jwst').submodule_search_locations[0], 'pipeline')
    for cfg in glob(os.path.join(cfg_dir, "*.cfg")):
        shutil.copy(cfg, dst)
