import logging

log = logging.getLogger(__name__)

__all__ = ["determine_ncores"]


def determine_ncores(max_cores, num_cores):
    """
    Determine the number of cores to use for multiprocessing.

    Parameters
    ----------
    max_cores : str or int
        Number of cores to use for multiprocessing. If set to 'none'
        (the default), then no multiprocessing will be done. The other
        allowable string values are 'quarter', 'half', and 'all', which indicate
        the fraction of cores to use for multi-proc. The total number of
        cores includes the SMT cores (Hyper Threading for Intel).
        If an integer is provided, it will be the exact number of cores used.
    num_cores : int
        Number of cores available on the machine

    Returns
    -------
    ncpus : int
        Number of cores to use for multiprocessing
    """
    match max_cores:
        case "none":
            return 1
        case None:
            return 1
        case "quarter":
            return num_cores // 4 or 1
        case "half":
            return num_cores // 2 or 1
        case "all":
            return num_cores
        case int():
            if max_cores <= num_cores and max_cores > 0:
                return max_cores
            log.warning(
                f"Requested {max_cores} cores exceeds the number of cores available "
                "on this machine ({num_cores}). Using all available cores."
            )
            return num_cores
        case _:
            raise ValueError(f"Invalid value for max_cores: {max_cores}")
