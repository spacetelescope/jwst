import numpy as np
import shapely.geometry
from stcal.alignment import sregion_to_footprint

__all__ = ["combine_sregions"]


def _combine_sregions(sregion_list, wcs, intersect=False):
    """
    Combine s_regions from input models to compute the s_region for the resampled data.

    Parameters
    ----------
    sregion_list : list of str
        List of s_regions from input models.
    wcs : gwcs.WCS
        WCS object for the resampled data.
        Needs to have "world" and "detector" frames, and the corresponding transforms
        must take in (RA, Dec) as their first two inputs
        and return (x, y) as their first two outputs.
    intersect : bool, optional
        If True, intersect the combined footprint with the WCS footprint.
        Default is False.

    Returns
    -------
    sregion : str
        The s_region for the resample data.
    """
    footprints = [sregion_to_footprint(sregion) for sregion in sregion_list]
    footprints = np.array(footprints)

    # convert to pixel coordinates to feed into polygon combine method
    # need to hack this a bit to ignore additional inputs and outputs
    # which are typically wavelength and order
    footprints_flat = footprints.reshape(-1, 2)
    world2det = wcs.get_transform("world", "detector")
    args = [1.0] * world2det.n_inputs
    args[0:2] = [footprints_flat[:, 0], footprints_flat[:, 1]]
    xy = world2det(*args)
    x, y = xy[0], xy[1]
    footprints_pixels = np.vstack([x, y]).T.reshape(footprints.shape)

    combined_polygons = _combine_footprints(footprints_pixels)

    # convert back to world coordinates
    det2world = wcs.get_transform("detector", "world")
    args = [1.0] * det2world.n_inputs
    combined_polygons_world = []
    for polygon in combined_polygons:
        args[0:2] = [polygon[:, 0], polygon[:, 1]]
        radec = det2world(*args)
        ra, dec = radec[0], radec[1]  # hack to ignore additional outputs
        combined_polygons_world.append(np.vstack([ra, dec]).T)

    # make s_region string
    sregion = _polygons_to_sregion(combined_polygons_world)
    if intersect:
        # TODO: add intersection with WCS footprint here
        # Relevant for resample when called on a custom WCS that may not cover
        # the full combined footprint
        pass
    return sregion


def combine_sregions(sregion_list, det2world, intersect=False):
    """
    Combine s_regions from input models to compute the s_region for the resampled data.

    Parameters
    ----------
    sregion_list : list of str
        List of s_regions from input models.
    det2world : astropy.modeling.CompoundModel
        WCS detector-to-world transform for the resampled data.
        Must take in exactly two inputs (x, y) and return exactly two outputs (RA, Dec).
        Must have a valid inverse transform.
    intersect : bool, optional
        If True, intersect the combined footprint with the WCS footprint.
        Default is False.

    Returns
    -------
    sregion : str
        The s_region for the resample data.
    """
    footprints = [sregion_to_footprint(sregion) for sregion in sregion_list]
    footprints = np.array(footprints)

    # convert to pixel coordinates to feed into polygon combine method
    footprints_flat = footprints.reshape(-1, 2)
    world2det = det2world.inverse
    x, y = world2det(footprints_flat[:, 0], footprints_flat[:, 1])
    footprints_pixels = np.vstack([x, y]).T.reshape(footprints.shape)

    combined_polygons = _combine_footprints(footprints_pixels)

    # convert back to world coordinates
    combined_polygons_world = []
    for polygon in combined_polygons:
        ra, dec = det2world(polygon[:, 0], polygon[:, 1])
        combined_polygons_world.append(np.vstack([ra, dec]).T)

    # make s_region string
    sregion = _polygons_to_sregion(combined_polygons_world)
    if intersect:
        # TODO: add intersection with WCS footprint here
        # Relevant for resample when called on a custom WCS that may not cover
        # the full combined footprint
        pass
    return sregion


def _polygons_to_sregion(polygons):
    """
    Create a FITS S_REGION from a list of polygons.

    Parameters
    ----------
    polygons : list[np.ndarray]
        List of polygons. Each polygon should have shape (V, 2), where V is the number of vertices.

    Returns
    -------
    str
        FITS S_REGION string.
    """
    s_region = ""
    for poly in polygons:
        poly = " ".join(f"{x:.9f} {y:.9f}" for x, y in poly)
        s_region += f"POLYGON ICRS  {poly}  "
    return s_region.strip()


def _combine_footprints(footprints):
    """
    Combine a list of footprints into one or more combined footprints.

    Parameters
    ----------
    footprints : list of np.ndarray
        List of footprints, where each footprint is a 2D array of shape (N, 2).

    Returns
    -------
    list of np.ndarray
        List of combined footprints, where each footprint is a 2D array of shape (M, 2).
    """
    footprints_shapely = [shapely.geometry.Polygon(footprint) for footprint in footprints]
    combined_footprints = shapely.unary_union(footprints_shapely)
    if isinstance(combined_footprints, shapely.geometry.Polygon):
        combined_footprints = [combined_footprints]
    elif isinstance(combined_footprints, shapely.geometry.MultiPolygon):
        combined_footprints = combined_footprints.geoms
    combined_polys = []
    for poly in combined_footprints:
        x, y = poly.exterior.coords.xy
        combined_poly = np.vstack([x, y]).T
        combined_poly = combined_poly[:-1]  # remove duplicate last point
        combined_poly = _simplify_by_angle(combined_poly)
        combined_polys.append(combined_poly)
    return combined_polys


def _simplify_by_angle(coords, point_thresh=1e-6, angle_thresh=1e-6):
    """
    Simplify a polygon by removing points that are collinear with their neighbors.

    Parameters
    ----------
    coords : np.ndarray
        2D array of shape (N, 2) representing the polygon vertices.

    Returns
    -------
    np.ndarray
        2D array of shape (M, 2) representing the simplified polygon vertices.
    """
    n = len(coords)
    if n < 3:
        return coords

    # Indices for previous, current, next points (wrap around)
    idx_prev = np.arange(-1, n - 1)
    idx_curr = np.arange(n)
    idx_next = np.arange(1, n + 1) % n

    p0 = coords[idx_prev]
    p1 = coords[idx_curr]
    p2 = coords[idx_next]

    # Vectors
    v1 = p1 - p0
    v2 = p2 - p1

    # Check closeness
    close1 = (np.abs(v1[:, 0]) < point_thresh) & (np.abs(v1[:, 1]) < point_thresh)
    close2 = (np.abs(v2[:, 0]) < point_thresh) & (np.abs(v2[:, 1]) < point_thresh)
    close = close1 | close2

    # Slopes
    m1 = v1[:, 1] / (v1[:, 0] + 1e-12)
    m2 = v2[:, 1] / (v2[:, 0] + 1e-12)
    delta_theta = np.arctan((m2 - m1) / (1 + m1 * m2))

    collinear = close | (np.abs(delta_theta) < angle_thresh)

    # Keep only non-collinear points
    return coords[~collinear]
