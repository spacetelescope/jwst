import pytest

from jwst.stpipe import integration


@pytest.fixture
def mock_stpipe_entry_points(monkeypatch):
    """
    Patch get_steps to return a single mock step.

    Parameters
    ----------
    monkeypatch : _pytest.monkeypatch.MonkeyPatch
        Pytest fixture to patch objects.
    """

    def get_steps():
        return [
            ("jwst.stpipe.tests.steps.AnotherDummyStep", "stpipe_dummy", False),
        ]

    monkeypatch.setattr(integration, "get_steps", get_steps)
