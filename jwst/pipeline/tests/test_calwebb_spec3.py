import pytest
import numpy as np
from astropy.io import fits
import stdatamodels.jwst.datamodels as dm
from jwst.pipeline.calwebb_spec3 import _save_wfss_x1d, _save_wfss_c1d


# TODO: decrease amount of code duplication in these two tests

def test_save_wfss_x1d(tmp_cwd):
    """Test flat file format saving of WFSS x1d data."""

    # set up a list of MultiSlitModel objects that look like outputs from extract_1d
    n_sources = 5
    n_exposures = 4
    n_rows = 3
    exposure_filenames = [f"exposure_{i}.fits" for i in range(n_exposures)]

    results_list = []
    for i in range(n_sources):
        multi = dm.MultiSpecModel()
        for j in range(n_exposures):
            # remove a single exposure, source pair to ensure filename matching works right
            if i == 1 and j == 2:
                continue
            # create a new SlitModel for each exposure
            spec = dm.SpecModel()
            spectable_dtype = spec.schema["properties"]["spec_table"]["datatype"]
            recarray_dtype = [(d["name"], d["datatype"]) for d in spectable_dtype]
            spec.meta.filename = exposure_filenames[j]
            spec.source_id = n_sources - i # reverse the order to test sorting
            spec.name = str(spec.source_id)
            spec_table = np.recarray((n_rows,), dtype=recarray_dtype)
            spec_table["WAVELENGTH"] = [1.0, 2.0, 3.0]
            spec_table["FLUX"] = [1.0, 2.0, 3.0]
            spec.spec_table = spec_table
            multi.spec.append(spec)
        results_list.append(multi)
    
    # save
    _save_wfss_x1d(results_list, "test_x1d.fits")

    # re-open as FITS and check the data
    hdul = fits.open("test_x1d.fits")
    assert len(hdul) == n_exposures + 2  # one primary HDU, one for each source, one ASDF extension
    for j in range(n_exposures):
        bintable = hdul[j + 1].data
        if j == 2:
            # Source 1 not observed in this exposure
            # Ensure the table is therefore one row shorter
            assert bintable.shape == (n_sources - 1,)
        else:
            assert bintable.shape == (n_sources,)
        
        # Ensure the source_ids are sorted
        source_ids = bintable["SOURCE_ID"]
        names = bintable["NAME"]
        assert np.all(np.diff(source_ids) >= 0)
        # Ensure names match the source_ids
        assert np.all([name == str(source_id) for name, source_id in zip(names, source_ids)])

    hdul.close()


def test_save_wfss_c1d(tmp_cwd):
    """Test flat file format saving of WFSS c1d data."""

    # Set up a list of MultiSlitModel objects that look like outputs from combine_1d
    n_sources = 5
    n_rows = 3

    results_list = []
    for i in range(n_sources):
        multi = dm.MultiCombinedSpecModel()
        # Only one spec per model in this case
        spec = dm.CombinedSpecModel()
        spectable_dtype = spec.schema["properties"]["spec_table"]["datatype"]
        recarray_dtype = [(d["name"], d["datatype"]) for d in spectable_dtype]
        spec.source_id = n_sources - i # reverse the order to test sorting
        spec_table = np.recarray((n_rows,), dtype=recarray_dtype)
        spec_table["WAVELENGTH"] = [1.0, 2.0, 3.0]
        spec_table["FLUX"] = [1.0, 2.0, 3.0]
        spec.spec_table = spec_table
        multi.spec.append(spec)
        results_list.append(multi)

    # Save
    _save_wfss_c1d(results_list, "test_c1d.fits")

    # Re-open as FITS and check the data
    hdul = fits.open("test_c1d.fits")
    assert len(hdul) == 3  # one primary HDU, one BinTableHDU, one ASDF extension
    bintable = hdul[1].data
    assert bintable.shape == (n_sources,)

    # Ensure the source_ids are sorted
    source_ids = bintable["SOURCE_ID"]
    assert np.all(np.diff(source_ids) >= 0)

    hdul.close()
