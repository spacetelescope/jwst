def pytest_addoption(parser):
    parser.addoption("--report-crds-context", action="store_true",
                     help="Report CRDS context in test suite header")


def pytest_report_header(config):
    """Add CRDS_CONTEXT to pytest report header"""

    if config.getoption("report_crds_context"):
        from stpipe.crds_client import get_context_used
        return f"crds_context: {get_context_used('jwst')}"
    else:
        return []
