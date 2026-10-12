import logging

log = logging.getLogger(__name__)


def check_cube_size(naxis1, naxis2, naxis3, limit):
    """Check the size of cube and throw an exception if too large."""
    total = naxis1 * naxis2 * naxis3
    if total > limit:
        msg = f"naxis1 * naxis2 * naxis3 = {total} exceeds limit of {limit}"
        log.error(msg)
        raise ValueError(msg)
