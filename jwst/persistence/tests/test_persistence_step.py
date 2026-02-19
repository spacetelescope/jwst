import asdf
import os
import numpy as np
import pytest

from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels import dqflags

from jwst.persistence import persistence
from jwst.persistence.persistence_step import PersistenceStep


@pytest.mark.skip(reason="testing int_times usage needs to be completed.")
def test_persistence_int_times(create_sci_model):
    nints, ngroups, nrows, ncols = 2, 7, 1, 2
    model = create_sci_model(nints=nints, ngroups=ngroups, nrows=nrows, ncols=ncols)
    model.groupdq[0, 5:, 0, 1] |= dqflags.group["SATURATED"]

    step = PersistenceStep(persistence_time=70)
    res = step.run(model)

    # XXX Add an int_times to the model.


def test_persistence_time_none_keeps_groupdq_unchanged(create_sci_model):
    """Test that when persistence_time is None, no PERSISTENCE flags are added to groupdq."""
    model = create_sci_model()
    # Pre-mark a pixel as SATURATED to prove no extra flagging occurs when disabled
    y, x = 1, 2
    model.groupdq[0, 0, y, x] |= dqflags.group["SATURATED"]
    step = PersistenceStep(persistence_time=None)
    res = step.run(model)
    assert np.array_equal(res.groupdq, model.groupdq)


def test_persistence_time_nonneg_sec(create_sci_model):
    """Test persitence flag gets set inside of persistence window"""
    nints, ngroups, nrows, ncols = 2, 7, 1, 2
    model = create_sci_model(nints=nints, ngroups=ngroups, nrows=nrows, ncols=ncols)
    model.groupdq[0, 5:, 0, 1] |= dqflags.group["SATURATED"]

    step = PersistenceStep(persistence_time=70)
    res = step.run(model)

    # With a persistence window of 70 seconds and group time of 21.47354 seconds, the 5th
    # group of the first integration for pixel (0, 1) and the following three groups
    # (crossing the integration) will be flagged as persistent.

    check1 = np.array([0, 0, 0, 0, 0, 0, 0], dtype=np.uint8)
    check2 = np.array([0, 0, 0, 0, 0, 0, 0], dtype=np.uint8)
    check3 = np.array([0, 0, 0, 0, 0, 34, 34], dtype=np.uint8)
    check4 = np.array([32, 32, 0, 0, 0, 0, 0], dtype=np.uint8)
    np.testing.assert_equal(res.groupdq[0, :, 0, 0], check1)
    np.testing.assert_equal(res.groupdq[1, :, 0, 0], check2)
    np.testing.assert_equal(res.groupdq[0, :, 0, 1], check3)
    np.testing.assert_equal(res.groupdq[1, :, 0, 1], check4)


def test_persistence_time_0_sec(create_sci_model):
    """Test persitence flag gets set inside of persistence window of 0 seconds"""
    nints, ngroups, nrows, ncols = 2, 7, 1, 2
    model = create_sci_model(nints=nints, ngroups=ngroups, nrows=nrows, ncols=ncols)
    model.groupdq[0, 5:, 0, 1] |= dqflags.group["SATURATED"]

    step = PersistenceStep(persistence_time=0)
    res = step.run(model)

    # With a zero persistence window, no groups should be flagged as persistent and
    # persistence_time gets set to None.
    assert step.persistence_time is None
    check1 = np.array([0, 0, 0, 0, 0, 0, 0], dtype=np.uint8)
    check2 = np.array([0, 0, 0, 0, 0, 0, 0], dtype=np.uint8)
    check3 = np.array([0, 0, 0, 0, 0, 2, 2], dtype=np.uint8)
    check4 = np.array([0, 0, 0, 0, 0, 0, 0], dtype=np.uint8)
    np.testing.assert_equal(res.groupdq[0, :, 0, 0], check1)
    np.testing.assert_equal(res.groupdq[1, :, 0, 0], check2)
    np.testing.assert_equal(res.groupdq[0, :, 0, 1], check3)
    np.testing.assert_equal(res.groupdq[1, :, 0, 1], check4)


def test_persistence_time_neg_sec(create_sci_model):
    """Test persitence flag gets set inside of negative persistence window."""
    nints, ngroups, nrows, ncols = 2, 7, 1, 2
    model = create_sci_model(nints=nints, ngroups=ngroups, nrows=nrows, ncols=ncols)
    model.groupdq[0, 5:, 0, 1] |= dqflags.group["SATURATED"]

    step = PersistenceStep(persistence_time=-100)

    res = step.run(model)

    # With a negative persistence window, no groups should be flagged as persistent and
    # persistence_time gets set to None.
    assert step.persistence_time is None
    check1 = np.array([0, 0, 0, 0, 0, 0, 0], dtype=np.uint8)
    check2 = np.array([0, 0, 0, 0, 0, 0, 0], dtype=np.uint8)
    check3 = np.array([0, 0, 0, 0, 0, 2, 2], dtype=np.uint8)
    check4 = np.array([0, 0, 0, 0, 0, 0, 0], dtype=np.uint8)
    np.testing.assert_equal(res.groupdq[0, :, 0, 0], check1)
    np.testing.assert_equal(res.groupdq[1, :, 0, 0], check2)
    np.testing.assert_equal(res.groupdq[0, :, 0, 1], check3)
    np.testing.assert_equal(res.groupdq[1, :, 0, 1], check4)


def test_persistence_time_with_array(create_sci_model, tmp_path):
    """Test persitence flag gets set using a persistence_array"""
    nints, ngroups, nrows, ncols = 2, 7, 1, 3
    model = create_sci_model(nints=nints, ngroups=ngroups, nrows=nrows, ncols=ncols)
    model.groupdq[0, 5:, 0, 1] |= dqflags.group["SATURATED"]

    asdf_file = os.path.join(tmp_path,  "persistence_array.asdf")
    persistence_array = np.array([[0.0, 0.0, 1648999345.9470801]], dtype=np.float64)
    tree = {"persistence_data": persistence_array}
    with asdf.AsdfFile(tree) as af:
        af.write_to(asdf_file)

    step = PersistenceStep(persistence_array_file=asdf_file, persistence_time=70)
    res = step.run(model)

    # With a non-zero persistence_array entry for pixel (0, 2), that pixel should have
    # persistence flagging applied based on the time in the array.
    checkz = np.array([0, 0, 0, 0, 0, 0, 0], dtype=np.uint8)
    check1 = np.array([0, 0, 0, 0, 0, 34, 34], dtype=np.uint8)
    check2 = np.array([32, 32, 0, 0, 0, 0, 0], dtype=np.uint8)
    np.testing.assert_equal(res.groupdq[0, :, 0, 0], checkz)
    np.testing.assert_equal(res.groupdq[1, :, 0, 0], checkz)
    np.testing.assert_equal(res.groupdq[0, :, 0, 1], check1)
    np.testing.assert_equal(res.groupdq[1, :, 0, 1], check2)
    np.testing.assert_equal(res.groupdq[0, :, 0, 2], check2)
    np.testing.assert_equal(res.groupdq[1, :, 0, 2], checkz)


def test_persistence_time_dnu(create_sci_model):
    """Test persitence flag gets set with DO_NOT_USE"""
    nints, ngroups, nrows, ncols = 2, 7, 1, 2
    model = create_sci_model(nints=nints, ngroups=ngroups, nrows=nrows, ncols=ncols)
    model.groupdq[0, 5:, 0, 1] |= dqflags.group["SATURATED"]

    step = PersistenceStep(persistence_dnu=True, persistence_time=70)
    res = step.run(model)

    # With persistence_dnu=True, the DO_NOT_USE flag should be set along with PERSISTENCE
    checkz = np.array([0, 0, 0, 0, 0, 0, 0], dtype=np.uint8)
    check1 = np.array([0, 0, 0, 0, 0, 35, 35], dtype=np.uint8)
    check2 = np.array([33, 33, 0, 0, 0, 0, 0], dtype=np.uint8)
    np.testing.assert_equal(res.groupdq[0, :, 0, 0], checkz)
    np.testing.assert_equal(res.groupdq[1, :, 0, 0], checkz)
    np.testing.assert_equal(res.groupdq[0, :, 0, 1], check1)
    np.testing.assert_equal(res.groupdq[1, :, 0, 1], check2)


def test_step_persistence_fails(monkeypatch, create_sci_model):
    # mock a known error condition in the persistence call
    monkeypatch.setattr(
        persistence.DataSet, "do_all", lambda *args: (datamodels.RampModel(), None, None, True)
    )

    sci = create_sci_model()
    result = PersistenceStep.call(sci)

    # Step is skipped
    assert result.meta.cal_step.persistence == "SKIPPED"

    # Input is not modified
    assert result is not sci
    assert sci.meta.cal_step.persistence is None
