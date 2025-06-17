"""Test for duplication/missing associations"""

import logging

from astropy.utils.data import get_pkg_data_filename

from jwst.associations.tests.helpers import (
    level3_rule_path,
    registry_level3_only,
)
from jwst.associations import AssociationPool, generate
from jwst.associations.main import Main
from jwst.associations.lib.utilities import constrain_on_candidates
from jwst.tests.helpers import LogWatcher


def test_duplicate_names(monkeypatch):
    """
    For Level 3 association, there should be no association
    with the same product name.
    """
    pool = AssociationPool.read(
        get_pkg_data_filename("data/jw00632_dups.csv", package="jwst.associations.tests")
    )
    constrain_all_candidates = constrain_on_candidates(None)
    rules = registry_level3_only(global_constraints=constrain_all_candidates)

    # watch for an error log message, we don't use caplog here because
    # something in the test suite messes up the logging during some runs
    # (probably stpipe or the association generator) and causes caplog
    # to sometimes miss the message
    watcher = LogWatcher(
        "Following associations have the same product name but significant differences"
    )
    monkeypatch.setattr(logging.getLogger("jwst.associations.lib.prune"), "warning", watcher)
    asns = generate(pool, rules)

    watcher.assert_seen()


def test_duplicate_generate():
    """Test for duplicate/overwrite association

    The pool has two exposures, one without a valid `asn_candidate`,
    and one with a valid observation `asn_candidate`.
    When set with the "all candidates" constraint, only one association
    should be made.

    The prompt for this test was that three associations were being created,
    two of which were the observation candidate, with the second
    being a duplicate of the first. The third was an extraneous
    discovered candidate.
    """
    pool = AssociationPool.read(
        get_pkg_data_filename("data/pool_duplicate.csv", package="jwst.associations.tests")
    )
    constrain_all_candidates = constrain_on_candidates(None)
    rules = registry_level3_only(global_constraints=constrain_all_candidates)
    asns = generate(pool, rules)
    assert len(asns) == 1
    asn = asns[0]
    assert asn["asn_type"] == "image3"
    assert asn["asn_id"] == "o029"


def test_duplicate_main():
    """Test the same but with Main."""
    cmd_args = [
        get_pkg_data_filename("data/pool_duplicate.csv", package="jwst.associations.tests"),
        "--dry-run",
        "-r",
        level3_rule_path(),
        "--ignore-default",
    ]

    generated = Main.cli(cmd_args)
    asns = generated.associations
    assert len(asns) == 2
    asn_types = set([asn["asn_type"] for asn in asns])
    assert len(asn_types) == 1
    assert "image3" in asn_types
    asn_ids = set([asn["asn_id"] for asn in asns])
    assert asn_ids == set(("a3001", "o029"))
