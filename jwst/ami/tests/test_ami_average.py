from jwst.ami.ami_average_step import AmiAverageStep


def test_step_init():
    """
    Just a simple smoke test to make sure the step
    can be created with the default spec.
    """
    step = AmiAverageStep()
    assert step
