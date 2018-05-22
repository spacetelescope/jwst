# Collect the pipeline configurations

import os
import shutil


def collect_pipeline_cfgs(dst='./'):
    # find installed pipeline *.cfg files
    import jwst.pipeline as p

    if not os.path.exists(dst):
        os.makedirs(dst)

    cfg_dir = p.__path__[0]
    for i in os.listdir(cfg_dir):
        if i.endswith('.cfg'):
            cfg = os.path.join(cfg_dir, i)
            shutil.copy(cfg, dst)
