"""
Polygon filling algorithm.

"""
# Original author: Nadezhda Dencheva
#
# modifications by Mihai Cara: removed functionality not needed for the
# skymatch algorithm and modified the code to be able to work with polygons
# that have vertices with negative coordinates. Polygon vertices are now
# *internally* (to region.py) rounded to integers so that Polygon will not
# crash when input vertices are floats. Fixed a bug in _construct_ordered_GET
# that was causing varying polygon filling for different ordering of the
# vertices. Finally, modified the algorithm to fill the right-most pixels
# as well as top-most row of the polygon.
#
# NOTE: Algorithm description can be found, e.g., here:
#    http://www.cs.rit.edu/~icss571/filling/how_to.html
#    http://www.cs.uic.edu/~jbell/CourseNotes/ComputerGraphics/PolygonFilling.html
#
from collections import OrderedDict
import numpy as np

__all__ = ['Region', 'Edge', 'Polygon']
__taskname__ = 'region'
__author__ = 'Nadezhda Dencheva, Mihai Cara'


class ValidationError(Exception):
    def __init__(self, message):
        self._message = message

    def __str__(self):
        return self._message


class Region():
    """
    Base class for regions.

    Parameters
    -------------
    rid : int or string
        region ID
    coordinate_system : astropy.wcs.CoordinateSystem instance or a string
        in the context of WCS this would be an instance of wcs.CoordinateSysem
    """

    def __init__(self, rid, coordinate_system):
        self._coordinate_system = coordinate_system
        self._rid = rid

    def __contains__(self, x, y):
        """
        Determines if a pixel is within a region.

        Parameters
        ----------
        x, y : float
            x , y values of a pixel

        Returns
        -------
        True or False

        Subclasses must define this method.
        """
        raise NotImplementedError("__contains__")

    def scan(self, mask):
        """
        Sets mask values to region id for all pixels within the region.
        Subclasses must define this method.

        Parameters
        ----------
        mask : ndarray
            a byte array with the shape of the observation to be used as a mask

        Returns
        -------
        mask : array where the value of the elements is the region ID or 0 (for
            pixels which are not included in any region).
        """
        raise NotImplementedError("scan")


class Polygon(Region):
    """
    Represents a 2D polygon region with multiple vertices

    Parameters
    ----------
    rid : string
         polygon id
    vertices : list of (x,y) tuples or lists
         The list is ordered in such a way that when traversed in a
         counterclockwise direction, the enclosed area is the polygon.
         The last vertex must coincide with the first vertex, minimum
         4 vertices are needed to define a triangle
    coord_system : string
        coordinate system

    """

    def __init__(self, rid, vertices, coord_system="Cartesian"):
        assert len(vertices) >= 4, ("Expected vertices to be "
                                    "a list of minimum 4 tuples (x,y)")
        super(Polygon, self).__init__(rid, coord_system)

        # self._shiftx & self._shifty are introduced to shift the bottom-left
        # corner of the polygon's bounding box to (0,0) as a (hopefully
        # temporary) workaround to a limitation of the original code that the
        # polygon must be completely contained in the image. It seems that the
        # code works fine if we make sure that the bottom-left corner of the
        # polygon's bounding box has non-negative coordinates.
        self._shiftx = 0
        self._shifty = 0
        for vertex in vertices:
            x, y = vertex
            if x < self._shiftx:
                self._shiftx = x
            if y < self._shifty:
                self._shifty = y
        v = [(i - self._shiftx, j - self._shifty) for i, j in vertices]

        # convert to integer coordinates:
        self._vertices = np.asarray(list(map(_round_vertex, v)))
        self._shiftx = int(round(self._shiftx))
        self._shifty = int(round(self._shifty))

        self._bbox = self._get_bounding_box()
        self._scan_line_range = \
            list(range(self._bbox[1], self._bbox[3] + self._bbox[1] + 1))
        # constructs a Global Edge Table (GET) in bbox coordinates
        self._GET = self._construct_ordered_GET()

    def _get_bounding_box(self):
        x = self._vertices[:, 0].min()
        y = self._vertices[:, 1].min()
        w = self._vertices[:, 0].max() - x
        h = self._vertices[:, 1].max() - y
        return x, y, w, h

    def _construct_ordered_GET(self):
        """
        Construct a Global Edge Table (GET)

        The GET is an OrderedDict. Keys are scan  line numbers,
        ordered from bbox.ymin to bbox.ymax, where bbox is the
        bounding box of the polygon.
        Values are lists of edges for which edge.ymin==scan_line_number.

        Returns
        -------
        GET: OrderedDict
            {scan_line: [edge1, edge2]}
        """
        # edges is a list of Edge objects which define a polygon
        # with these vertices
        edges = self.get_edges()
        GET = OrderedDict.fromkeys(self._scan_line_range)
        ymin = np.asarray([e._ymin for e in edges])
        for i in self._scan_line_range:
            ymin_ind = (ymin == i).nonzero()[0]
            # a hack for incomplete filling .any() fails if 0 is in ymin_ind
            # if ymin_ind.any():
            yminindlen, = ymin_ind.shape
            if yminindlen:
                GET[i] = [edges[ymin_ind[0]]]
                for j in ymin_ind[1:]:
                    GET[i].append(edges[j])
        return GET

    def get_edges(self):
        """
        Create a list of Edge objects from vertices
        """
        edges = []
        for i in range(1, len(self._vertices)):
            name = 'E' + str(i - 1)
            edges.append(Edge(name=name, start=self._vertices[i - 1], stop=self._vertices[i]))
        return edges

    def scan(self, data):
        """
        This is the main function which scans the polygon and creates the mask

        Parameters
        ----------
        data : array
            the mask array
            it has all zeros initially, elements within a region are set to
            the region's ID

        Algorithm:
        - Set the Global Edge Table (GET)
        - Set y to be the smallest y coordinate that has an entry in GET
        - Initialize the Active Edge Table (AET) to be empty
        - For each scan line:
          1. Add edges from GET to AET for which ymin==y
          2. Remove edges from AET fro which ymax==y
          3. Compute the intersection of the current scan line with all edges in the AET
          4. Sort on X of intersection point
          5. Set elements between pairs of X in the AET to the Edge's ID

        """
        # TODO: 1.This algorithm does not mark pixels in the top row and left most column.
        # Pad the initial pixel description on top and left with 1 px to prevent this.
        # 2. Currently it uses intersection of the scan line with edges. If this is
        # too slow it should use the 1/m increment (replace 3 above) (or the increment
        # should be removed from the GET entry).

        # see comments in the __init__ function for the reason of introducing
        # polygon shifts (self._shiftx & self._shifty). Here we need to shift
        # it back.

        (ny, nx) = data.shape

        y = np.min(list(self._GET.keys()))

        AET = []
        scline = self._scan_line_range[-1]

        while y <= scline:

            if y < scline:
                AET = self.update_AET(y, AET)

            if self._bbox[2] <= 0:
                y += 1
                continue

            scan_line = Edge('scan_line', start=[self._bbox[0], y],
                             stop=[self._bbox[0] + self._bbox[2], y])
            x = [int(np.ceil(e.compute_AET_entry(scan_line)[1])) for e in AET if e is not None]
            xnew = np.sort(x)
            ysh = y + self._shifty

            if ysh < 0 or ysh >= ny:
                y += 1
                continue

            for i, j in zip(xnew[::2], xnew[1::2]):
                xstart = max(0, i + self._shiftx)
                xend = min(j + self._shiftx, nx - 1)
                data[ysh][xstart:xend + 1] = self._rid

            y += 1

        return data

    def update_AET(self, y, AET):
        """
        Update the Active Edge Table (AET)

        Add edges from GET to AET for which ymin of the edge is
        equal to the y of the scan line.
        Remove edges from AET for which ymax of the edge is
        equal to y of the scan line.

        """
        edge_cont = self._GET[y]
        if edge_cont is not None:
            for edge in edge_cont:
                if edge._start[1] != edge._stop[1] and edge._ymin == y:
                    AET.append(edge)
        for edge in AET[::-1]:
            if edge is not None:
                if edge._ymax == y:
                    AET.remove(edge)
        return AET

    def __contains__(self, px):
        """even-odd algorithm or smth else better should be used"""
        # minx = self._vertices[:,0].min()
        # maxx = self._vertices[:,0].max()
        # miny = self._vertices[:,1].min()
        # maxy = self._vertices[:,1].max()
        return px[0] >= self._bbox[0] and px[0] <= self._bbox[0] + self._bbox[2] and \
            px[1] >= self._bbox[1] and px[1] <= self._bbox[1] + self._bbox[3]


class Edge():
    """
    Edge representation

    An edge has a "start" and "stop" (x,y) vertices and an entry in the
    GET table of a polygon. The GET entry is a list of these values:

    [ymax, x_at_ymin, delta_x/delta_y]

    """

    def __init__(self, name=None, start=None, stop=None, next=None):
        self._start = None
        if start is not None:
            self._start = np.asarray(start)
        self._name = name
        self._stop = stop
        if stop is not None:
            self._stop = np.asarray(stop)
        self._next = next

        if self._stop is not None and self._start is not None:
            if self._start[1] < self._stop[1]:
                self._ymin = self._start[1]
                self._yminx = self._start[0]
            else:
                self._ymin = self._stop[1]
                self._yminx = self._stop[0]
            self._ymax = max(self._start[1], self._stop[1])
            self._xmin = min(self._start[0], self._stop[0])
            self._xmax = max(self._start[0], self._stop[1])
        else:
            self._ymin = None
            self._yminx = None
            self._ymax = None
            self._xmin = None
            self._xmax = None
        self.GET_entry = self.compute_GET_entry()

    @property
    def ymin(self):
        return self._ymin

    @property
    def start(self):
        return self._start

    @property
    def stop(self):
        return self._stop

    @property
    def ymax(self):
        return self._ymax

    def compute_GET_entry(self):
        """
        Compute the entry in the Global Edge Table

        [ymax, x@ymin, 1/m]

        """
        if self._start is None:
            entry = None
        else:
            earr = np.asarray([self._start, self._stop])
            if np.diff(earr[:, 1]).item() == 0:
                return None
            else:
                entry = [self._ymax, self._yminx, (np.diff(earr[:, 0]) / np.diff(earr[:, 1])).item(), None]
        return entry

    def compute_AET_entry(self, edge):
        """
        Compute the entry for an edge in the current Active Edge Table

        [ymax, x_intersect, 1/m]
        note: currently 1/m is not used
        """
        x = self.intersection(edge)[0]
        return [self._ymax, x, self.GET_entry[2]]

    def __repr__(self):
        fmt = ""
        if self._name is not None:
            fmt += self._name
            next = self.next
            while next is not None:
                fmt += "-->"
                fmt += next._name
                next = next.next
        return fmt

    @property
    def next(self):
        return self._next

    @next.setter
    def next(self, edge):
        if self._name is None:
            self._name = edge._name
            self._stop = edge._stop
            self._start = edge._start
            self._next = edge.next
        else:
            self._next = edge

    def intersection(self, edge):
        u = self._stop - self._start
        v = edge._stop - edge._start
        w = self._start - edge._start
        D = np.cross(u, v)

        if np.allclose(np.cross(u, v), 0, rtol=0,
                       atol=1e2 * np.finfo(float).eps):
            return np.array(self._start)

        return np.cross(v, w) / D * u + self._start

    def is_parallel(self, edge):
        u = self._stop - self._start
        v = edge._stop - edge._start
        return np.allclose(np.cross(u, v), 0, rtol=0,
                           atol=1e2 * np.finfo(float).eps)


def _round_vertex(v):
    x, y = v
    return int(round(x)), int(round(y))
