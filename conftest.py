"""Project default for pytest"""
import pytest


# Add option to run slow tests
def pytest_addoption(parser):
    parser.addoption(
        "--runslow",
        action="store_true",
        help="run slow tests"
    )
