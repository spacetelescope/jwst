import numpy as np

from jwst.background.weighted_sigma_clip import WeightedSigmaClip


def test_weighted_sigma_clip_1d():
    """Test that sigma clipper respects weights, handles NaNs, and masks out outlier distribution."""
    rng = np.random.default_rng(42)

    n_inlier = 1000
    n_outlier = 100
    n_nan = 20  # NaN values placed at random
    n_skew = 500  # pixels at a skewing value with very low weight

    # inliers drawn from N(0, 1)
    inlier_data = rng.normal(0, 1, size=n_inlier)
    inlier_weights = np.ones(n_inlier)

    # outliers drawn from N(±10, 0.5)
    outlier_data = np.concatenate(
        [
            rng.normal(-10, 0.5, size=n_outlier // 2),
            rng.normal(10, 0.5, size=n_outlier // 2),
        ]
    )
    outlier_weights = np.ones(n_outlier)

    # Skewing pixels at 3.0: with low weight the mean stays near 0
    skew_center = 3.0
    skew_data = rng.normal(skew_center, 3.0, size=n_skew)
    skew_weights = np.full(n_skew, 0.001)

    N = n_inlier + n_outlier + n_skew

    # Concatenate groups; index ranges are known before NaN injection
    data = np.concatenate([inlier_data, outlier_data, skew_data])
    weights = np.concatenate([inlier_weights, outlier_weights, skew_weights])

    inlier_idx = np.arange(n_inlier)
    outlier_idx = np.arange(n_inlier, n_inlier + n_outlier)

    # Randomly distribute NaNs across all distributions
    nan_idx = rng.choice(N, size=n_nan, replace=False)
    data[nan_idx] = np.nan

    clipper = WeightedSigmaClip(sigma=3.0)
    mask = ~clipper(data, weights)

    # NaN pixels must be masked
    assert np.all(mask[nan_idx]), "NaN pixels should be masked"

    # The vast majority of non-NaN inliers should survive (not be masked)
    non_nan_inlier_idx = np.setdiff1d(inlier_idx, nan_idx)
    inlier_pass_rate = np.sum(~mask[non_nan_inlier_idx]) / len(non_nan_inlier_idx)
    assert inlier_pass_rate > 0.95, "Most inliers should remain unmasked"

    # The vast majority of non-NaN outliers should be clipped
    non_nan_outlier_idx = np.setdiff1d(outlier_idx, nan_idx)
    outlier_clip_rate = np.sum(mask[non_nan_outlier_idx]) / len(non_nan_outlier_idx)
    assert outlier_clip_rate > 0.90, "Most outliers should be masked"


def test_sigma_clip_along_axis():
    """Test more realistic 2-D case where sigma clip is in cross-dispersion direction."""

    pass
