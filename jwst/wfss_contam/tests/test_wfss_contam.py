import pytest
from jwst.wfss_contam.wfss_contam import _determine_multiprocessing_ncores

@pytest.mark.parametrize("max_cores, num_cores, expected", 
                         [("none", 4, 1), 
                          ("quarter", 4, 1), 
                          ("half", 4, 2), 
                          ("all", 4, 4), 
                          ("none", 1, 1),])
def test_determine_multiprocessing_ncores(max_cores, num_cores, expected):
    assert _determine_multiprocessing_ncores(max_cores, num_cores) == expected  


