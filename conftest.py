from stpipe.crds_client import get_context_used


def pytest_report_header(config):
    """Add CRDS_CONTEXT to pytest report header"""
    return f"crds_context: {get_context_used('jwst')}"
