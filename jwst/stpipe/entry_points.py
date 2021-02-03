from pkg_resources import iter_entry_points
from collections import namedtuple
import warnings


STEPS_GROUP = "stpipe.steps"


StepInfo = namedtuple("StepInfo", ["class_name", "class_alias", "is_pipeline", "package_name", "package_version"])


def get_steps():
    """
    Get the list of steps registered with stpipe's entry point group.  Each entry
    point is expected to return a list of tuples, where the first tuple element
    is a fully-qualified Step subclass name, the second element is an optional
    class alias, and the third is a bool indicating whether the class is to be
    listed as a pipeline in the CLI output.

    Returns
    -------
    list of StepInfo
    """
    steps = []

    for entry_point in iter_entry_points(group=STEPS_GROUP):
        package_name = entry_point.dist.project_name
        package_version = entry_point.dist.version
        package_steps = []

        try:
            elements = entry_point.load()()

            for element in elements:
                package_steps.append(StepInfo(*element, package_name, package_version))
        except Exception as e:
            warnings.warn(
                f"{STEPS_GROUP} plugin from package {package_name}=={package_version} "
                f"failed to load:\n\n"
                f"{e.__class__.__name__}: {e}"
            )

        steps.extend(package_steps)

    return steps
