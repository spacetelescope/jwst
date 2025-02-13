# Collect the pipeline configurations

from pathlib import Path
import shutil
from importlib.util import find_spec


def collect_pipeline_cfgs(dst="./"):
    """
    Copy step and pipeline .cfg files to destination.

    Parameters
    ----------
    dst : str, optional
        Destination directory. Default is current directory.
    """
    Path(dst).mkdir(parents=True, exist_ok=True)

    cfg_dir = Path(find_spec("jwst").submodule_search_locations[0]) / "pipeline"
    for cfg in cfg_dir.glob("*.cfg"):
        shutil.copy(cfg, dst)
