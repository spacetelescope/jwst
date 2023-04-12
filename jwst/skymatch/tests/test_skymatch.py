from itertools import product
from copy import deepcopy

import pytest
import numpy as np
from astropy.modeling.models import Shift
from gwcs import WCS

from stdatamodels.jwst.datamodels import ImageModel, dqflags

from jwst.datamodels import ModelContainer
from jwst.assign_wcs import AssignWcsStep
from jwst.skymatch import SkyMatchStep
from jwst.tweakreg.utils import adjust_wcs
from jwst.associations.asn_from_list import asn_from_list
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base


DO_NOT_USE = dqflags.pixel['DO_NOT_USE']
SATURATED = dqflags.pixel['SATURATED']


@pytest.fixture
def nircam_rate():
    # Theoretical NIRCAM image size
    xsize0 = 2040
    ysize0 = 2040

    # Actual image size for this simulation:
    xsize = 96
    ysize = 96

    shape = (ysize, xsize)
    im = ImageModel(shape)
    im.var_rnoise += 0

    im.meta.wcsinfo = {
        'ctype1': 'RA---TAN',
        'ctype2': 'DEC--TAN',
        'dec_ref': 11.99875540218638,
        'ra_ref': 22.02351763251896,
        'roll_ref': 0.005076934167039675,
        'v2_ref': 86.039011,
        'v3_ref': -493.385704,
        'v3yangle': -0.07385127,
        'vparity': -1,
        'wcsaxes': 2
    }
    im.meta.instrument = {
        'channel': 'LONG',
        'detector': 'NRCALONG',
        'filter': 'F444W',
        'lamp_mode': 'NONE',
        'module': 'A',
        'name': 'NIRCAM',
        'pupil': 'CLEAR'
    }
    im.meta.subarray = {
        'fastaxis': -1,
        'name': 'FULL',
        'slowaxis': 2,
        'xsize': xsize,
        'xstart': 1,
        'ysize': ysize,
        'ystart': 1
    }
    im.meta.observation = {
        'activity_id': '01',
        'date': '2021-10-25',
        'exposure_number': '00001',
        'obs_id': 'V42424001001P0000000001101',
        'observation_label': 'nircam_ptsrc_only',
        'observation_number': '001',
        'program_number': '42424',
        'sequence_id': '1',
        'time': '16:58:27.258',
        'visit_group': '01',
        'visit_id': '42424001001',
        'visit_number': '001'
    }
    im.meta.exposure = {
        'duration': 161.05155,
        'end_time': 59512.70899968495,
        'exposure_time': 150.31478,
        'frame_time': 10.73677,
        'group_time': 21.47354,
        'groupgap': 1,
        'integration_time': 150.31478,
        'mid_time': 59512.70812980775,
        'nframes': 1,
        'ngroups': 7,
        'nints': 1,
        'nresets_at_start': 1,
        'nresets_between_ints': 1,
        'readpatt': 'BRIGHT1',
        'sample_time': 10,
        'start_time': 59512.70725993055,
        'type': 'NRC_IMAGE'
    }
    im.meta.photometry = {
        'pixelarea_steradians': 1e-13,
        'pixelarea_arcsecsq': 4e-3,
    }

    im = AssignWcsStep.call(im, sip_approx=False)

    # Adjust WCS to account for reduced image size (for performance reasons):
    wcs = im.meta.wcs
    pipeline = deepcopy(wcs.pipeline)

    s = pipeline[0]
    shiftx = (xsize0 - xsize) / 2.0
    shifty = (ysize0 - ysize) / 2.0
    s.transform = (Shift(shiftx) & Shift(shifty)) | s.transform

    wcs = WCS(pipeline)
    wcs.bounding_box = ((-0.5, xsize - 0.5), (-0.5, ysize - 0.5))
    im.meta.wcs = wcs

    return im


def _add_bad_pixels(im, sat_val, dont_use_val):
    # Add two types of "bad" pixels: 1) in the center of the image that will
    # lie in the intersection region common to all images (we assume images
    # are rotated around a common center) and 2) bad pixels at the corners
    # of the images that have a different flag and will be excluded from the
    # analysis of rotated images because they will lie outside of the
    # intersection region common to all images (and not because they are
    # masked out based on DQ array).
    mask = np.ones(im.data.shape, dtype=bool)
    # Add some "bad" pixels:
    # corners
    im.data[:5, :5] = sat_val
    im.data[-5:, :5] = sat_val
    im.data[-5:, -5:] = sat_val
    im.data[:5, -5:] = sat_val

    im.dq[:5, :5] = SATURATED
    im.dq[-5:, :5] = SATURATED
    im.dq[-5:, -5:] = SATURATED
    im.dq[:5, -5:] = SATURATED

    mask[:5, :5] = False
    mask[-5:, :5] = False
    mask[-5:, -5:] = False
    mask[:5, -5:] = False

    cy, cx = (x // 2 for x in im.data.shape)
    cx -= 5
    cy -= 5

    # center
    im.data[cx:cx + 10, cy:cy + 10] = dont_use_val
    im.dq[cx:cx + 10, cy:cy + 10] = DO_NOT_USE
    mask[cx:cx + 10, cy:cy + 10] = False

    return im, mask


@pytest.mark.parametrize(
    'skymethod, subtract, skystat, match_down, grouped',
    tuple(
        product(
            ['local', 'match', 'global', 'global+match'],
            [False, True],
            ['median', 'mean', 'midpt', 'mode'],
            [False, True],
            [False, True]
        )
    )
)
def test_skymatch(nircam_rate, skymethod, subtract, skystat, match_down,
                  grouped):
    # test basic functionality and correctness of sky computations
    np.random.seed(1)
    im1 = nircam_rate.copy()
    im2 = im1.copy()
    im3 = im1.copy()

    # add "bad" data
    im1, dq_mask = _add_bad_pixels(im1, 1e6, 1e9)
    im2, _ = _add_bad_pixels(im2, 5e6, 3e9)
    im3, _ = _add_bad_pixels(im3, 7e6, 1e8)

    if grouped:
        im2.meta.observation.sequence_id = "1"
        im3.meta.observation.sequence_id = "1"
    else:
        im2.meta.observation.sequence_id = "2"
        im3.meta.observation.sequence_id = "3"

    container = ModelContainer([im1, im2, im3])

    # define some background:
    levels = [9.12, 8.28, 2.56]

    for im, lev in zip(container, levels):
        im.data += np.random.normal(
            loc=lev,
            scale=0.1,
            size=im.data.shape
        )

    # exclude central DO_NOT_USE and corner SATURATED pixels
    result = SkyMatchStep.call(
        container,
        skymethod=skymethod,
        match_down=match_down,
        subtract=subtract,
        skystat=skystat,
        binwidth=0.2,
        nclip=0,
        dqbits='~DO_NOT_USE+SATURATED'
    )

    if skymethod == 'match' and grouped:
        # nothing to "match" when there is only one group:
        assert im.meta.background.method is None
        assert im.meta.background.subtracted is None

        # test that output models have original sky levels on failure:
        for im, lev in zip(result, levels):
            assert abs(np.mean(im.data[dq_mask]) - lev) < 0.01

        return

    if skymethod in ['local', 'global+match']:
        if grouped:
            ref_levels = len(levels) * [min(levels)]
        else:
            ref_levels = levels

    elif skymethod == 'match':
        if grouped:
            lev0 = 0
        else:
            lev0 = min(levels) if match_down else max(levels)
        ref_levels = np.subtract(levels, lev0)

    elif skymethod == 'global':
        ref_levels = len(levels) * [min(levels)]

    sub_levels = np.subtract(levels, ref_levels)

    for im, lev, rlev, slev in zip(result, levels, ref_levels, sub_levels):
        # check that meta was set correctly:
        assert im.meta.background.method == skymethod
        assert im.meta.background.subtracted == subtract

        # test computed/measured sky values:
        assert abs(im.meta.background.level - rlev) < 0.01

        # test
        if subtract:
            assert abs(np.mean(im.data[dq_mask]) - slev) < 0.01
        else:
            assert abs(np.mean(im.data[dq_mask]) - lev) < 0.01


@pytest.mark.parametrize(
    'skymethod, subtract, skystat',
    tuple(
        product(
            ['local', 'match', 'global'],
            [False, True],
            ['mean']
        )
    )
)
def test_skymatch_overlap(nircam_rate, skymethod, subtract, skystat):
    # test that computations are performed only in the area of overlap
    # between images (bad pixels in the corners of rotated images are ignored).
    # Set 'nclip' to 0 in order to not clip bad pixels in computing mean.
    np.random.seed(1)
    im1a = nircam_rate.copy()
    im1b = im1a.copy()
    im2a = im1a.copy()
    im2b = im1a.copy()
    im3 = im1a.copy()

    # add "bad" data
    im1a, dq_mask = _add_bad_pixels(im1a, 1e6, 1e9)
    im1b, _ = _add_bad_pixels(im1b, 1e6, 1e9)
    im2a, _ = _add_bad_pixels(im2a, 5e6, 3e9)
    im2b, _ = _add_bad_pixels(im2b, 5e6, 3e9)
    im3, _ = _add_bad_pixels(im3, 7e6, 1e8)

    # rotate images so that corner SATURATED pixels are not in the overlap
    # region:
    im1a.meta.wcs = adjust_wcs(im2a.meta.wcs, delta_ra=0.1, delta_dec=0.1)
    im2a.meta.wcs = adjust_wcs(im2a.meta.wcs, delta_roll=30)
    im2b.meta.wcs = adjust_wcs(
        im2a.meta.wcs,
        delta_ra=0.1,
        delta_dec=0.1,
        delta_roll=30
    )
    im3.meta.wcs = adjust_wcs(im3.meta.wcs, delta_roll=60)

    im2a.meta.observation.sequence_id = "2"
    im2b.meta.observation.sequence_id = "2"
    im3.meta.observation.sequence_id = "3"

    container = ModelContainer([im1a, im1b, im2a, im2b, im3])

    # define some background:
    levels = [9.12, 9.12, 8.28, 8.28, 2.56]

    for im, lev in zip(container, levels):
        im.data += np.random.normal(loc=lev, scale=0.1, size=im.data.shape)

    # We do not exclude SATURATED pixels. They should be ignored because
    # images are rotated and SATURATED pixels in the corners are not in the
    # common intersection of all input images. This is the purpose of this test
    result = SkyMatchStep.call(
        container,
        skymethod=skymethod,
        match_down=True,
        subtract=subtract,
        skystat=skystat,
        nclip=0,
        dqbits='~DO_NOT_USE'  # specifically DO NOT add 'SATURATED' flag
    )

    if skymethod in ['local', 'global+match']:
        ref_levels = levels

    elif skymethod == 'match':
        ref_levels = np.subtract(levels, min(levels))

    elif skymethod == 'global':
        ref_levels = len(levels) * [min(levels)]

    sub_levels = np.subtract(levels, ref_levels)

    for im, lev, rlev, slev in zip(result, levels, ref_levels, sub_levels):
        # check that meta was set correctly:
        assert im.meta.background.method == skymethod
        assert im.meta.background.subtracted == subtract

        if skymethod in ['local', 'global']:
            # These two sky methods must fail because they do not take
            # into account (do not compute) overlap regions and use
            # entire images:

            # test computed/measured sky values:
            assert abs(im.meta.background.level - rlev) > 1000  # FAIL

            # test
            if subtract:
                assert abs(np.mean(im.data[dq_mask]) - slev) > 1000  # FAIL
            else:
                assert abs(np.mean(im.data[dq_mask]) - lev) < 0.01  # PASS
        else:
            # test computed/measured sky values:
            assert abs(im.meta.background.level - rlev) < 0.01

            # test
            if subtract:
                assert abs(np.mean(im.data[dq_mask]) - slev) < 0.01
            else:
                assert abs(np.mean(im.data[dq_mask]) - lev) < 0.01


def test_asn_input(nircam_rate, tmpdir):
    # This is the same test as 'test_skymatch_overlap' with
    # skymethod='match', subtract=True, skystat='mean' and with memory saving
    # feature enabled (data loaded from files as needed).
    np.random.seed(1)
    im1 = nircam_rate.copy()
    im2 = im1.copy()
    im3 = im1.copy()

    # add "bad" data
    im1, dq_mask = _add_bad_pixels(im1, 1e6, 1e9)
    im2, _ = _add_bad_pixels(im2, 5e6, 3e9)
    im3, _ = _add_bad_pixels(im3, 7e6, 1e8)

    # rotate images so that corner SATURATED pixels are not in the overlap
    # region:
    im2.meta.wcs = adjust_wcs(im2.meta.wcs, delta_roll=30)
    im3.meta.wcs = adjust_wcs(im3.meta.wcs, delta_roll=60)

    im2.meta.observation.sequence_id = "2"
    im3.meta.observation.sequence_id = "3"

    container = ModelContainer([im1, im2, im3])

    # define some background:
    levels = [9.12, 8.28, 2.56]

    for im, lev in zip(container, levels):
        im.data += np.random.normal(loc=lev, scale=0.1, size=im.data.shape)

    im1_path = str(tmpdir / "skymatch_im1.fits")
    im2_path = str(tmpdir / "skymatch_im2.fits")
    im3_path = str(tmpdir / "skymatch_im3.fits")

    im1.write(im1_path)
    im2.write(im2_path)
    im3.write(im3_path)

    assoc_out = asn_from_list(
        [im1_path, im2_path, im3_path],
        rule=DMS_Level3_Base,
        product_name='skymatch'
    )
    asn_out_fname, out_serialized = assoc_out.dump(format='json')
    asn_out_fname = str(tmpdir / asn_out_fname)
    with open(asn_out_fname, "w") as asn_out:
        asn_out.write(out_serialized)

    # We do not exclude SATURATED pixels. They should be ignored because
    # images are rotated and SATURATED pixels in the corners are not in the
    # common intersection of all input images. This is the purpose of this test
    step = SkyMatchStep(
        minimize_memory=True,
        skymethod='match',
        match_down=True,
        subtract=True,
        skystat='mean',
        nclip=0,
        dqbits='~DO_NOT_USE'  # specifically DO NOT add 'SATURATED' flag
    )

    result = step.run(asn_out_fname)

    assert isinstance(result, str)

    ref_levels = np.subtract(levels, min(levels))
    sub_levels = np.subtract(levels, ref_levels)

    result = ModelContainer(result)

    for im, lev, rlev, slev in zip(result, levels, ref_levels, sub_levels):
        # check that meta was set correctly:
        assert im.meta.background.method == 'match'
        assert im.meta.background.subtracted is True

        # test computed/measured sky values:
        assert abs(im.meta.background.level - rlev) < 0.01

        # test
        assert abs(np.mean(im.data[dq_mask]) - slev) < 0.01


@pytest.mark.parametrize(
    'skymethod, subtract',
    tuple(
        product(
            ['local', 'match', 'global', 'global+match'],
            [False, True]
        )
    )
)
def test_skymatch_2x(nircam_rate, tmpdir, skymethod, subtract):
    # Test that repetitive applications of skymatch produce the same results
    np.random.seed(1)
    im1 = nircam_rate.copy()
    im2 = im1.copy()
    im3 = im1.copy()

    # add "bad" data
    im1, dq_mask = _add_bad_pixels(im1, 1e6, 1e9)
    im2, _ = _add_bad_pixels(im2, 5e6, 3e9)
    im3, _ = _add_bad_pixels(im3, 7e6, 1e8)

    im2.meta.observation.sequence_id = "2"
    im3.meta.observation.sequence_id = "3"

    container = ModelContainer([im1, im2, im3])

    # define some background:
    levels = [9.12, 8.28, 2.56]

    for im, lev in zip(container, levels):
        im.data += np.random.normal(loc=lev, scale=0.1, size=im.data.shape)

    im1_path = str(tmpdir / "skymatch_im1.fits")
    im2_path = str(tmpdir / "skymatch_im2.fits")
    im3_path = str(tmpdir / "skymatch_im3.fits")

    im1.write(im1_path)
    im2.write(im2_path)
    im3.write(im3_path)

    assoc_out = asn_from_list(
        [im1_path, im2_path, im3_path],
        rule=DMS_Level3_Base,
        product_name='skymatch'
    )
    asn_out_fname, out_serialized = assoc_out.dump(format='json')
    asn_out_fname = str(tmpdir / asn_out_fname)
    with open(asn_out_fname, "w") as asn_out:
        asn_out.write(out_serialized)

    # We do not exclude SATURATED pixels. They should be ignored because
    # images are rotated and SATURATED pixels in the corners are not in the
    # common intersection of all input images. This is the purpose of this test
    step = SkyMatchStep(
        minimize_memory=True,
        skymethod=skymethod,
        match_down=True,
        subtract=subtract,
        skystat='mean',
        nclip=0,
        dqbits='~DO_NOT_USE+SATURATED'
    )

    result = step.run(asn_out_fname)

    with pytest.raises(ValueError):
        step.subtract = not subtract
        step.run(result)

    # 2nd run.
    step.subtract = subtract
    result2 = step.run(result)

    # compute expected levels
    if skymethod in ['local', 'global+match']:
        ref_levels = levels

    elif skymethod == 'match':
        ref_levels = np.subtract(levels, min(levels))

    elif skymethod == 'global':
        ref_levels = len(levels) * [min(levels)]

    sub_levels = np.subtract(levels, ref_levels)

    result2 = ModelContainer(result2)

    # compare results
    for im2, lev, rlev, slev in zip(result2, levels, ref_levels, sub_levels):
        # check that meta was set correctly:
        assert im2.meta.background.method == skymethod
        assert im2.meta.background.subtracted == subtract

        # test computed/measured sky values:
        assert abs(im2.meta.background.level - rlev) < 0.01

        # test
        if subtract:
            assert abs(np.mean(im2.data[dq_mask]) - slev) < 0.01
        else:
            assert abs(np.mean(im2.data[dq_mask]) - lev) < 0.01
