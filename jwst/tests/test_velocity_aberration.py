"""
Test script for set_velocity_aberration.py
"""
from numpy import isclose
import os
import sys

sys.path.insert(
    0,
    os.path.join(os.path.dirname(__file__), '../../scripts')
)

import set_velocity_aberration as sva # noqa: E402

# Testing constants
GOOD_VELOCITY = (100.0, 100.0, 100.0)
GOOD_POS = (0., 0.)
GOOD_SCALE_FACTOR = 1.000333731048419
GOOD_OFFSET_X = 0.00033356409519815205
GOOD_OFFSET_Y = 0.00033356409519815205

ZERO_VELOCITY = 0.
ZERO_SCALE_FACTOR = 1.0
ZERO_OFFSET_X = 0.
ZERO_OFFSET_Y = 0.


def test_scale_factor_valid():
    scale_factor = sva.aberration_scale(
        GOOD_VELOCITY[0], GOOD_VELOCITY[1], GOOD_VELOCITY[2],
        GOOD_POS[0], GOOD_POS[1]
    )
    assert isclose(scale_factor, GOOD_SCALE_FACTOR)


def test_scale_factor_zero_velocity():
    scale_factor = sva.aberration_scale(
        ZERO_VELOCITY, ZERO_VELOCITY, ZERO_VELOCITY,
        GOOD_POS[0], GOOD_POS[1]
    )
    assert isclose(scale_factor, ZERO_SCALE_FACTOR)


def test_offset_valid():
    delta_x, delta_y = sva.aberration_offset(
        GOOD_VELOCITY[0], GOOD_VELOCITY[1], GOOD_VELOCITY[2],
        GOOD_POS[0], GOOD_POS[1]
    )
    assert isclose(delta_x, GOOD_OFFSET_X)
    assert isclose(delta_y, GOOD_OFFSET_Y)


def test_offset_zero_velocity():
    delta_x, delta_y = sva.aberration_offset(
        ZERO_VELOCITY, ZERO_VELOCITY, ZERO_VELOCITY,
        GOOD_POS[0], GOOD_POS[1]
    )
    assert isclose(delta_x, ZERO_OFFSET_X)
    assert isclose(delta_y, ZERO_OFFSET_Y)
