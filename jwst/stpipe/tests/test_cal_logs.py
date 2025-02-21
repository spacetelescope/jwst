from jwst.stpipe.tests.steps import CalLogsStep, CalLogsPipeline


def test_cal_logs_step():
    m = CalLogsStep().run("foo")
    assert any(("foo" in l for l in m.cal_logs.cal_logs_step))


def test_cal_logs_pipeline():
    m = CalLogsPipeline().run("foo")
    assert not hasattr(m.cal_logs, "cal_logs_step")
    assert any(("foo" in l for l in m.cal_logs.cal_logs_pipeline))
