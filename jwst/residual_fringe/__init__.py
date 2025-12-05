"""Correct residual fringes in MIRI MRS data."""

from .residual_fringe_step import ResidualFringeStep
from .utils import fit_residual_fringes_1d

__all__ = ["ResidualFringeStep", "fit_residual_fringes_1d"]
