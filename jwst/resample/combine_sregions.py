import numpy as np
import shapely.geometry
from stcal.alignment import sregion_to_footprint

__all__ = ["combine_sregions"]


def combine_sregions(sregion_list, det2world, intersect_footprint=None):
    """
    Combine s_regions from input models to compute the s_region for the resampled data.

    Parameters
    ----------
    sregion_list : list of str
        List of s_regions from input models.
    det2world : `~astropy.modeling.Model`
        WCS detector-to-world transform for the resampled data.
        Must take in exactly two inputs (x, y) and return exactly two outputs (RA, Dec).
        Must have a valid inverse transform.
    intersect_footprint : np.ndarray, optional
        Footprint of the output WCS in world coordinates, shape (N, 2).
        If provided, the combined footprint from the input s_region list
        will be intersected with this footprint.

    Returns
    -------
    str
        The combined s_region.
    """
    footprints = np.array([sregion_to_footprint(sregion) for sregion in sregion_list])

    # convert from world to pixel coordinates
    footprints_flat = footprints.reshape(-1, 2)
    world2det = det2world.inverse
    x, y = world2det(footprints_flat[:, 0], footprints_flat[:, 1])
    footprints_pixels = np.vstack([x, y]).T.reshape(footprints.shape)

    # combine footprints with Shapely
    combined_polygons = _combine_footprints(footprints_pixels)

    # intersect with output WCS footprint
    if intersect_footprint is not None:
        x, y = world2det(intersect_footprint[:, 0], intersect_footprint[:, 1])
        intersect_footprint_pixels = np.vstack([x, y]).T
        final_polygons = _intersect_with_bbox(combined_polygons, intersect_footprint_pixels)
        if not final_polygons:
            raise ValueError("No overlap between input s_regions and intersection footprint")
    else:
        final_polygons = combined_polygons

    # convert back from pixel to world coordinates
    combined_polygons_world = []
    for polygon in final_polygons:
        ra, dec = det2world(polygon[:, 0], polygon[:, 1])
        combined_polygons_world.append(np.vstack([ra, dec]).T)

    # turn lists of indices into a single S_REGION string
    sregion = _polygons_to_sregion(combined_polygons_world)

    return sregion


def _polygons_to_sregion(polygons):
    """
    Create an S_REGION from a list of polygons.

    Parameters
    ----------
    polygons : list[np.ndarray]
        List of polygons. Each polygon should have shape (V, 2), where V is the number of vertices.
        V can be different for each polygon.

    Returns
    -------
    str
        S_REGION string.
    """
    s_region = ""
    for poly in polygons:
        poly = " ".join(f"{x:.9f} {y:.9f}" for x, y in poly)
        s_region += f"POLYGON ICRS  {poly}  "
    return s_region.strip()


def _combine_footprints(footprints):
    """
    Combine a list of footprints into one or more combined footprints using a Shapely union.

    Parameters
    ----------
    footprints : list of np.ndarray
        List of footprints, where each footprint is a 2-D array of shape (N, 2).

    Returns
    -------
    list of np.ndarray
        List of combined footprints, where each footprint is a 2-D array of shape (M, 2).
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


def _intersect_with_bbox(polygons, bbox):
    """
    Intersect a list of polygons with a bounding box.

    Parameters
    ----------
    polygons : list[np.ndarray]
        List of polygons. Each polygon should have shape (V, 2), where V is the number of vertices.
    bbox : np.ndarray
        2D array of shape (N, 2) representing the bounding box vertices.

    Returns
    -------
    np.ndarray
        2D array of shape (M, 2) representing the intersected polygon vertices.
    """
    intersect_polygon = shapely.geometry.Polygon(bbox)
    final_polygons = []
    for polygon in polygons:
        polygon = shapely.geometry.Polygon(polygon)
        intersection = shapely.intersection(polygon, intersect_polygon)
        if not intersection.is_empty:
            x, y = intersection.exterior.coords.xy
            poly_out = np.vstack([x, y]).T
            poly_out = poly_out[:-1]  # remove duplicate last point
            poly_out = _simplify_by_angle(poly_out)
            final_polygons.append(poly_out)
    return final_polygons


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
    # Indices for previous, current, next points (wrap around)
    n = len(coords)
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
