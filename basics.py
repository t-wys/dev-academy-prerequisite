'''General geometrical objects: point, vector, segment, polygon,
polyhedron and plane curve. These objects are base for special cases like
Rectangle or right pyramid.

Note
----
Instances of classes in this module are not intended to be updated after
instantiation. Reason for this limitation is the scope of the task,
which is limited to calcuations of surface area or volume, but no manipulations
on these objects.
'''
import numbers
import warnings
from math import sqrt, acos
try:
    # Supported in Python 3.5+
    from math import isclose
except ImportError:
    def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
        '''Check whether two float numbers are close.
        Taken from PEP 485.

        Parameters
        ----------
        a, b : float
            Numbers to be compared
        rel_tol, abs_tol : float, optional
            Relative and absolute tolerance

        Returns
        -------
        bool
            True if numbers are close
        '''
        return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


class Point(object):
    '''Basic 0D object, building block of everything else.'''
    def __init__(self, x, y, z=0):
        '''Creates new point from at least two numbers.

        Parameters
        ----------
        x, y : number
            At least X and Y coordinates are mandatory
        z : number, optional
            Z coordinate is optional, by default = 0

        Raises
        ------
        TypeError
            If any of coordinates is not a number
        '''
        if not all(isinstance(coord, numbers.Number) for coord in (x, y, z)):
            raise TypeError('X, Y, or Z is not a number')
        self._coords = (x, y, z)

    def __repr__(self):
        return f'Point ({self.x:.2f}, {self.y:.2f}, {self.z:.2f})'

    def __getitem__(self, idx):
        return self._coords[idx]

    def __eq__(self, point):
        return all(c0 == c1 for c0, c1 in zip(self, point))

    @classmethod
    def from_str(cls, string):
        '''Allow to create new point from a string which
        follow specific formatting

        Parameters
        ----------
        string : str
            Expected format: X,Y or X,Y,Z where X, Y and Z
            are digits with (optional) decimal mark. Example:
            3,4.56 or 1.1,2.2,3
        '''
        return cls(*(float(coord) for coord in string.split(',')))

    @property
    def x(self):
        '''Access to X coordinate'''
        return self._coords[0]

    @property
    def y(self):
        '''Access to Y coordinate'''
        return self._coords[1]

    @property
    def z(self):
        '''Access to Z coordinate'''
        return self._coords[2]

    def as_tuple(self):
        '''Access to tuple with (x, y, z) coordinates'''
        return self._coords

    def xy_distance_to(self, point):
        '''Euclidean distance between two XYZ points projected onto XY plane.

        Parameters
        ----------
        point : Point or iterable
            The second point for calculation

        Returns
        -------
        float
        '''
        return sqrt((self.x - point.x)**2 + (self.y - point.y)**2)

    def distance_to(self, point):
        '''Euclidean distance between two points given by X, Y, Z coordinates.

        Parameters
        ----------
        point : Point or iterable
            The second point for calculation

        Returns
        -------
        float
        '''
        return sqrt(sum([(c0 - c1)**2 for c0, c1 in zip(self, point)]))


class Vector(object):
    '''Representation of 3D vector

    Attributes
    ----------
    components : tuple
        Collection of X, Y, Z components
    x, y, z : number
        Components along X, Y and Z axes
    length : float
        sqrt(sum(i**2)) for i in (x, y, z)

    Methods
    -------
    from_points
        Creates new vector from start point (anchor) and end point.
    dot
        Calculates scalar product.
    cross
        Calculates cross product.
    '''
    def __init__(self, x, y, z=0):
        '''Creates new vector and calculates its basic properties

        Parameters
        ----------
        x, y : number
            X and Y components are mandatory
        z : number
            Z component is optional
        '''
        self._components = (x, y, z)
        self._length = sqrt(sum(c**2 for c in self._components))

    @classmethod
    def from_points(cls, start_point, end_point):
        '''Define vector not by its components but derive components
        as a difference between coordinated of end and start points.

        Parameters
        ----------
        start_point, end_point : Point
            starting point p and ending point p+r where r is a vector
            that will be returned by this method

        Returns
        -------
        Vector
        '''
        return cls(*(e - s for s, e in zip(start_point, end_point)))

    def __repr__(self):
        return f'Vector [{self.x}, {self.y}, {self.z}]'

    def __getitem__(self, idx):
        return self._components[idx]

    @property
    def x(self):
        '''Access to X component'''
        return self._components[0]

    @property
    def y(self):
        '''Access to Y component'''
        return self._components[1]

    @property
    def z(self):
        '''Access to Z component'''
        return self._components[2]

    @property
    def components(self):
        '''Access to tuple of components'''
        return self._components

    @property
    def length(self):
        '''Length of a vector'''
        return self._length

    def dot(self, vector):
        '''Dot product (scalar product):
        v1 (dot) v2 = sum(v1{i}*v2{i}) for i in [x, y, z]

        Parameters
        ----------
        vector : Vector
            another vector

        Returns
        -------
        number
            scalar product of two vectors
        '''
        return sum(c0*c1 for c0, c1 in zip(self, vector))

    def cross(self, vector):
        '''Cross product

        Parameters
        ----------
        vector : Vector
            another vector

        Returns
        -------
        Vector
            vector product of two vectors
        '''
        return Vector(
            self.y*vector.z - self.z*vector.y,
            self.z*vector.x - self.x*vector.z,
            self.x*vector.y - self.y*vector.x,
        )

    def __add__(self, vector):
        return Vector(*(c0 + c1 for c0, c1 in zip(self, vector)))

    def __neg__(self):
        return Vector(*(-c for c in self))

    def __sub__(self, vector):
        return self + (-vector)

    def __mul__(self, scalar):
        return Vector(*(scalar*c for c in self))

    def __rmul__(self, scalar):
        return self * scalar


class Segment(object):
    '''Representation of a line segment.

    Attributes
    ----------
    start_point, end_point : Point
        Segment spans between these two point 
    length : float
        Length ot segment
    xy_projection_length : float
        Length of segment's projection onto XY plane

    Methods
    -------
    intersects_with
        Performs vector-based intersection check with other segment
    '''
    def __init__(self, p0, p1):
        '''Creates new segment given two points. Also calculates
        basic properties.

        Parameters
        ----------
        p0, p1 : Point
            TODO
        '''
        self._p0 = p0
        self._p1 = p1
        self._length = self._p0.distance_to(self._p1)
        self._xy_length = self._p0.xy_distance_to(self._p1)

    def __repr__(self):
        return f'Segment between {self._p0} and {self._p1}'

    @property
    def length(self):
        '''Segment length

        Returns
        -------
        length of the segment : float
        '''
        return self._length

    @property
    def xy_projection_length(self):
        '''Length of the projection onto XY plane.

        Returns
        -------
        length of projection : float
        '''
        return self._xy_length

    @property
    def start_point(self):
        '''Access to starting point'''
        return self._p0

    @property
    def end_point(self):
        '''Access to end point'''
        return self._p1

    def intersects_with(self, segment):
        '''Checks whether the segment intersects with another one.
        This function operates only on XY plane.

        This follows solution posted on Stack Overflow by Gareth Rees:
        https://stackoverflow.com/a/565282

        Parameters
        ----------
        segment : Segment
            segment to compare againts

        Returns
        -------
        bool
        '''
        origin_1 = Vector.from_points(Point(0, 0), self._p0)
        origin_2 = Vector.from_points(Point(0, 0), segment.start_point)
        end_1 = Vector.from_points(self._p0, self._p1)
        end_2 = Vector.from_points(segment.start_point, segment.end_point)

        u_11 = (origin_2 - origin_1).cross(end_1)
        u_12 = (origin_2 - origin_1).cross(end_2)
        u_2 = end_2.cross(end_1)

        # Colinear, may overlap
        if u_11.z == 0 and u_2.z == 0:
            t_0 = ((origin_2 - origin_1).dot(end_1)) / end_1.dot(end_1)
            t_1 = t_0 + end_2.dot(end_1) / end_1.dot(end_1)
            if 0 < t_0 < 1 or 0 < t_1 < 1:
                return True
        # Parallel and non-intersecting
        elif u_11.z != 0 and u_2.z == 0:
            return False
        # Intersect
        elif u_2.z != 0:
            # For some reason I need "-" to make this work as expected.
            param_1 = -u_11.z / u_2.z
            param_2 = -u_12.z / u_2.z
            # Originally it is 0 <= param_1 <= 1 and 0 <= param_2 <= 1,
            # but single common end-start point is not considered here as
            # an intersection.
            if 0 < param_1 < 1 and 0 < param_2 < 1:
                return True
        return False


class PlaneCurve(object):
    '''Base class for plane curves (line, circle, ellipse, ...).

    Task specifies only one curve (ellipse) and there is no general
    formula for perimeter or surface area (in some cases it does not exist!).
    Thus, this class is a stub and base class for `Ellipse`.
    '''
    def __init__(self):
        pass

    def __repr__(self):
        return 'Plane curve'


class Polygon(object):
    '''A non-intersecting polygon (simple polygon).

    Attributes
    ----------
    perimeter : number
        Sum of lengths of edges
    surface_area : number
        Self explanatory
    is_regular : bool
        Indicates whether all edges and angles are equal
    vertices : list of `Point`s
        List of points on which span the polygon
    edges : list of `Segment`s
        List of line segments that span on pairs of consecutive vertices
    angles : list of floats
        List of angles between to neighbouring segments. Expressed in radians.
    '''
    def __init__(self, vertices):
        '''Creates new polygon

        Parameters
        ----------
        vertices : list of `Point`s

        Raises
        ------
        ValueError
            If number of vertices if smaller than 3
        '''
        if len(vertices) < 3:
            raise ValueError('Passed %s vertices. Expecting at least 3' % len(vertices))
        self._vertices = vertices.copy()
        if any(vertex.z != 0 for vertex in self._vertices):
            warnings.warn(
                'One or more vertices has non-zero Z coordinate. '\
                'Note that Z values will be ignored in calculations.',
                Warning
            )
        self._edges = [
            Segment(prev_vertex, next_vertex)
            for prev_vertex, next_vertex in
            zip(self._vertices, (*self._vertices[1:], self._vertices[0]))
        ]
        self._validate_polygon()
        self._angles = []
        for prev_vertex, cur_vertex, next_vertex in zip(
                (self._vertices[-1], *self._vertices[:-1]),
                self._vertices,
                (*self._vertices[1:], self._vertices[0])
        ):
            prev_vector = Vector.from_points(prev_vertex, cur_vertex)
            next_vector = Vector.from_points(cur_vertex, next_vertex)
            self._angles.append(
                acos(prev_vector.dot(next_vector) / (prev_vector.length * next_vector.length))
            )

        self._perimeter = self._calculate_perimeter()
        self._surface_area = self._calculate_surface_area()
        self._regular = all(
            isclose(edge.length, self._edges[0].length) for edge in self._edges
        ) and all(
            isclose(angle, self._angles[0]) for angle in self._angles
        )

    def __repr__(self):
        return f'Polygon {len(self._vertices)}-vertices'

    def _validate_polygon(self):
        '''Checks whether segments defined by vertices do not intersect.

        Raises
        ------
        Exception
            If an intersection is detected
        '''
        for i, edge_1 in enumerate(self._edges[:-1]):
            for edge_2 in self._edges[(i+1):]:
                if edge_1.intersects_with(edge_2):
                    raise Exception(
                        'Found intersection between two edges. '\
                        'Make sure that vertices are given in right order '\
                        '(either clockwise or counter-clockwise).'
                    )
        return True

    def _calculate_perimeter(self):
        '''Perimeter -- sum of lengths of segments.

        Returns
        -------
        float
        '''
        return sum([edge.xy_projection_length for edge in self._edges])

    def _calculate_surface_area(self):
        '''Calculates perimeter based on general formula (vertices coordinates).
        The formula is valid for simple polygons.

        Formula gives signed area (sign depends on how the list of points
        is traversed: clockwise or counter-clockwise), hence abs()

        Returns
        -------
        float
        '''
        return abs(0.5*sum([
            prev_vertex.x*next_vertex.y - next_vertex.x*prev_vertex.y
            for prev_vertex, next_vertex in
            zip(self._vertices, (*self._vertices[1:], self._vertices[0]))
        ]))

    @property
    def perimeter(self):
        '''Perimeter of the polygon

        Returns
        -------
        number
        '''
        return self._perimeter

    @property
    def surface_area(self):
        '''Surface area of the polygon

        Returns
        -------
        number
        '''
        return self._surface_area

    @property
    def is_regular(self):
        '''Indicates whether the polygon has all egdes and angles equal.

        Returns
        -------
        bool
        '''
        return self._regular

    @property
    def vertices(self):
        '''List of vertices.

        Returns
        -------
        list of `Point`s
        '''
        return self._vertices

    @property
    def edges(self):
        '''List of egdes.

        Returns
        -------
        list of `Segment`s
        '''
        return self._edges

    @property
    def angles(self):
        '''List of angles (in radians)

        Returns
        -------
        list of floats
        '''
        return self._angles


class Polyhedron(object):
    '''3D object composed of two bases and faces.
    This is not the most general case.

    Attributes
    ----------
    surface_area : number
        Self-explanatory
    volume : number
        Self-explanatory
    '''
    def __init__(self, base, height):
        '''Creates a polyhedron with given base and height

        Parameters
        ----------
        base : Polygon
            Instance of `Polygon` class
        height : number
            Distance between two bases

        Raises
        ------
        TypeError
            If base is not an instance of `Polygon` class
        ValueError
            If height is not a number or is not greater than zero
        '''
        if not isinstance(base, Polygon):
            raise TypeError('Base must be an instance of `Polygon`')
        if not isinstance(height, numbers.Number):
            raise TypeError('Height must be a numeric type')
        if height <= 0:
            raise ValueError('Height must be a positive number')
        self._base = base
        self._height = height

        self._surface_area = self._calculate_surface_area()
        self._volume = self._calculate_volume()

    def __repr__(self):
        return f'Polyhedron with {self._base} base'

    def _calculate_surface_area(self):
        '''Surface area, just:
        surf.area = 2 * area{base} + perimeter{base} * height
        '''
        return 2*self._base.surface_area + self._base.perimeter*self._height

    def _calculate_volume(self):
        '''Volume, just: volume = area{base} * height'''
        return self._base.surface_area * self._height

    @property
    def surface_area(self):
        '''Access surface area of polyhedron'''
        return self._surface_area

    @property
    def volume(self):
        '''Access volume of polyhedron'''
        return self._volume
