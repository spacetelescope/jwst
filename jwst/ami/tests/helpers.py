"""Useful constants and functions for AMI tests."""

import numpy as np

PXSC_DEG = 65.6 / (60.0 * 60.0 * 1000)
PXSC_RAD = PXSC_DEG * np.pi / (180)
PXSC_MAS = PXSC_DEG * 3600 * 1000

__all__ = []  # type: ignore[var-annotated]
