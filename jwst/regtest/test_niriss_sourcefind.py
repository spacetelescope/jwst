import pytest
from astropy.io import ascii
from numpy.testing import assert_allclose

from stdatamodels.jwst import datamodels
from jwst.tweakreg import tweakreg_catalog


@pytest.mark.bigdata
@pytest.mark.parametrize(
    "starfinder",
    ["iraf", "segmentation"],
)
def test_tweakreg_catalog_starfinder_alternatives(rtdata, starfinder):
    """
    Test that the IRAF and segmentation star finders give expected results for undersampled NIRISS data
    It is well known that DAOStarFinder gives bad results so is not included in this test
    """

    stem = "jw01088003001_01101_00005"
    rtdata.get_data(f"niriss/imaging/{stem}_nis_cal.fits")
    model = datamodels.ImageModel(rtdata.input)
    catalog = tweakreg_catalog.make_tweakreg_catalog(
        model,
        2.5,
        2.5,
        10.0,
        starfinder_name=starfinder,
        starfinder_kwargs={
            "brightest": None,
            "sharphi": 3.0,
            "minsep_fwhm": 2.5,
            "sigma_radius": 2.5,
        },
    )

    output_name = f"{stem}_{starfinder}_cat.ecsv"
    catalog.write(output_name, format="ascii.ecsv")
    rtdata.output = output_name
    rtdata.get_truth(f"truth/test_niriss_sourcefind/{stem}_{starfinder}_cat.ecsv")
    catalog_truth = ascii.read(rtdata.truth)

    # rtol is larger than default because of numerical differences on Linux vs MacOS
    assert_allclose(catalog["xcentroid"], catalog_truth["xcentroid"], rtol=1e-3)
    assert_allclose(catalog["ycentroid"], catalog_truth["ycentroid"], rtol=1e-3)
