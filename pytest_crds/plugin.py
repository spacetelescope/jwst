def pytest_addoption(parser):
    parser.addoption("--report-crds-context", action="store_true",
                     help="Report CRDS context in test suite header")


def pytest_report_header(config):
    """Add CRDS_CONTEXT to pytest report header"""

    if config.getoption("report_crds_context"):
        from stpipe.crds_client import get_context_used
        try:
            context = get_context_used('jwst')
        except SystemExit:
            context = 'unknown'
        return f"crds_context: {context}"
    else:
        return []
