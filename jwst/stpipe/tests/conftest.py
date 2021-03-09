import pytest

from jwst.stpipe import integration


@pytest.fixture
def mock_stpipe_entry_points(monkeypatch):
    def get_steps():
        return [
            ("jwst.stpipe.tests.steps.AnotherDummyStep", "stpipe_dummy", False),
        ]
    monkeypatch.setattr(integration, "get_steps", get_steps)
