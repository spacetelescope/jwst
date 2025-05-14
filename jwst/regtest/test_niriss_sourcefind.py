import pytest
import numpy as np

from stdatamodels.jwst import datamodels
from jwst.tweakreg import tweakreg_catalog


@pytest.mark.bigdata
def test_tweakreg_catalog_starfinder_overlap(rtdata):
    '''
    Test that the IRAF and segmentation star finders find overlapping stars
    within a given tolerance for xcentroid and ycentroid.
    '''

    stem = "jw01088003001_01101_00005"
    rtdata.get_data(f"niriss/imaging/{stem}_nis_cal.fits")
    model = datamodels.ImageModel(rtdata.input)

    # Generate catalogs for both starfinder techniques
    catalogs = {}
    for starfinder in ["iraf", "segmentation"]:
        catalogs[starfinder] = tweakreg_catalog.make_tweakreg_catalog(
            model, 2.5, kernel_fwhm=2.5, bkg_boxsize=400.0, starfinder_name=starfinder, starfinder_kwargs={
                'brightest': None,
                'sharphi': 3.0,
                'minsep_fwhm': 2.5,
                'sigma_radius': 2.5,
            })

    # Compare the two catalogs
    iraf_x = np.array(catalogs["iraf"]["xcentroid"])
    iraf_y = np.array(catalogs["iraf"]["ycentroid"])
    seg_x = np.array(catalogs["segmentation"]["xcentroid"])
    seg_y = np.array(catalogs["segmentation"]["ycentroid"])


    tolerance = 0.5  # pixel tolerance for matching centroids
    matches = 0

    # Compute pairwise differences. dx, dy, and distances have shape (n_iraf, n_segmentation)
    dx = iraf_x[:, None] - seg_x[None, :]
    dy = iraf_y[:, None] - seg_y[None, :]

    # Calculate distances and find matches within the tolerance
    distances = np.sqrt(dx**2 + dy**2)
    matches = np.sum(np.any(distances <= tolerance, axis=1))

    # With these parameters, IRAF finds ~4300 stars and segmentation finds ~5300 stars
    assert matches > 4000
