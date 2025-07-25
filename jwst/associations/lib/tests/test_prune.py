from copy import deepcopy

import pytest

from jwst.associations.asn_from_list import asn_from_list
from jwst.associations.lib.prune import prune_remove
from jwst.associations.lib.rules_level2b import Asn_Lv2Image


def test_prune_duplicate():
    """
    Test pruning when two equivalent associations are in a list.

    Here, the second rateints association is the duplicate and should
    be removed. The others should stay.
    """
    rate = asn_from_list(["jw00016006001_04105_00001_nrca4_rate.fits"], rule=Asn_Lv2Image)
    rateints = asn_from_list(["jw00016006001_04105_00001_nrca4_rateints.fits"], rule=Asn_Lv2Image)
    rateints_duplicate = deepcopy(rateints)

    # all three asns are equal, according to comparison rules
    assert rate == rateints
    assert rateints == rateints_duplicate
    assert rateints is not rateints_duplicate

    asn_list = [rate, rateints, rateints_duplicate]
    to_remove = [rateints_duplicate]
    known_dups = []

    prune_remove(asn_list, to_remove, known_dups)

    assert len(asn_list) == 2
    assert asn_list[0] is rate
    assert asn_list[1] is rateints


def test_prune_duplicate_not_found():
    """Test pruning error condition."""
    rate = asn_from_list(["jw00016006001_04105_00001_nrca4_rate.fits"], rule=Asn_Lv2Image)
    rateints = asn_from_list(["jw00016006001_04105_00001_nrca4_rateints.fits"], rule=Asn_Lv2Image)
    rateints_duplicate = deepcopy(rateints)
    asn_list = [rate, rateints, rateints_duplicate]

    # Specify an equivalent association to remove, but it
    # is not the exact same object as the one in the list.
    # This should raise an error.
    to_remove = [deepcopy(rateints_duplicate)]
    known_dups = []

    with pytest.raises(IndexError, match="Expected duplicate association not found"):
        prune_remove(asn_list, to_remove, known_dups)
