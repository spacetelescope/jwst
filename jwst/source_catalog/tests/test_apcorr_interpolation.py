from jwst.source_catalog import SourceCatalogStep
from jwst.source_catalog.tests.helpers import make_nircam_model, mock_apcorr


def test_nondefault_aperture_ee_is_interpolated():
    model = make_nircam_model()
    apcorr = mock_apcorr()

    step = SourceCatalogStep(
        aperture_ee1=35,
        aperture_ee2=50,
        aperture_ee3=70,
        snr_threshold=0.5,
        npixels=5,
        bkg_boxsize=50,
        kernel_fwhm=2.0,
        save_results=False,
        override_apcorr=apcorr,
    )

    catalog = step.run(model)

    assert catalog is not None
    assert "aper35_flux" in catalog.colnames
    assert "aper50_flux" in catalog.colnames
    assert "aper70_flux" in catalog.colnames

    model.close()
