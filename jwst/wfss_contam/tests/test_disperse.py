import numpy as np
import pytest
from numpy.testing import assert_allclose

from jwst.wfss_contam.disperse import disperse


@pytest.mark.parametrize("phot_per_lam", [True, False])
def test_disperse_oversample_same_result(grism_wcs, segmentation_map, phot_per_lam):
    """Coverage for bug where wavelength oversampling led to double-counted fluxes."""
    x0 = np.array([300.5])
    y0 = np.array([300.5])
    order = 1
    flxs = np.array([1.0])
    source_id = np.array([0])
    naxis = (300, 500)
    sens_waves = np.linspace(1.708, 2.28, 100)
    wmin, wmax = np.min(sens_waves), np.max(sens_waves)
    sens_resp = np.ones(100)
    seg_wcs = segmentation_map.meta.wcs
    (
        0,
        (300, 500),
        2,
        False,
    )
    xoffset = 2200
    yoffset = 1000

    output_images = []
    for os in [2, 3]:
        src = disperse(
            x0,
            y0,
            flxs,
            source_id,
            order,
            wmin,
            wmax,
            sens_waves,
            sens_resp,
            seg_wcs,
            grism_wcs,
            naxis,
            oversample_factor=os,
            xoffset=xoffset,
            yoffset=yoffset,
            phot_per_lam=phot_per_lam,
        )
        output_images.append(src[source_id[0]]["image"])

    # different oversampling gives different effects at the ends
    # unsure if this is a bug or not, but the middle should definitely be the same
    assert_allclose(output_images[0][2:-2, :], output_images[1][2:-2], rtol=1e-5)
