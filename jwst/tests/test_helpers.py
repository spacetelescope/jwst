import logging
import warnings

import pytest

from jwst.tests.helpers import _help_pytest_warns


def test_warning_inside_capture_logging():
    """
    Using _help_pytest_warns with pytest.warns doesn't prevent catching
    expected warnings.
    """
    with _help_pytest_warns(), pytest.warns(UserWarning, match="expected"):
        warnings.warn("expected", UserWarning)


@pytest.mark.xfail
def test_unexpected_warning_inside_capture_logging():
    """
    This should fail (see xfail) since the warning doesn't match.
    We're using an xfail here because we're testing that pytest.warns
    fails to catch a warnings.
    """
    with _help_pytest_warns(), pytest.warns(UserWarning, match="expected"):
        warnings.warn("not the warning you're looking for", UserWarning)


@pytest.mark.filterwarnings("ignore:foo")
def test_python_warning_logging_incompatibility():
    """
    An example and smoke test for python showing the incompatibility
    of logging.captureWarnings used within a warnings.catch_warnings
    context.
    """
    with warnings.catch_warnings(record=True) as rec:
        logging.captureWarnings(True)
        warnings.warn("foo", UserWarning)
        logging.captureWarnings(False)

    assert len(rec) == 0  # note this is 0 and not 1
