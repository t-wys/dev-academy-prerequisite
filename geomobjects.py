'''Geometrical objects with some special properties: rectangle, square,
triangle, ellipse, cube and right pyramid.

Note
----
Instances of classes in this module are not intended to be updated after
instantiation. Reason for this limitation is the scope of the task,
which is limited to calcuations of surface area or volume, but no manipulations
on these objects.
'''
import numbers
import warnings
from math import sqrt, sin, cos, pi
from basics import Polygon, Polyhedron, PlaneCurve, Point, Segment

_SCIPY_SUPPORT = True

try:
    from scipy.special import ellipe
except ImportError:
    warnings.warn('Could not import scipy. Perimeter of ellipse will be approximated.')
    _SCIPY_SUPPORT = False


class RegularPolyhedron(Polyhedron):
    '''Polyhedron with regular polygons as bases.

    See docs for `Polyhedron` for a list of methods and attributes.
    '''
    def __init__(self, base, height):
        if not isinstance(base, RegularPolygon) or not base.is_regular:
            raise ValueError('Base of regular polyhedron must be an instance of `RegularPolygon`')
        super().__init__(base, height)

    @classmethod
    def from_base_parameters(cls, num_vertices, circumcircle_radius, height):
        '''TODO'''
        return cls(RegularPolygon(num_vertices, circumcircle_radius), height)


class RightPyramid(RegularPolyhedron):
    '''Polyhedron with single base, which is a regular polygon.
    Side faces are triangles that meet at the pyramid's apex.

    See docs for `Polyhedron` for a list of methods and attributes.
    '''
    def _calculate_surface_area(self):
        '''TODO'''
        return (
            self._base.surface_area
            + 0.5 * self._base.perimeter * sqrt(
                self._height**2 + self._base.inscribed_radius**2
            ))

    def _calculate_volume(self):
        '''TODO'''
        return super()._calculate_volume() / 3


class Rectangle(Polygon):
    '''A polygon with four sides and four right angles.

    See docs for `Polygon` for a list of methods and attributes.
    '''
    def __init__(self, a, b):
        '''Create rectangle with sides of length `a` and `b`

        Parameters
        ----------
        a, b : number
            Length of the edges of the rectangle
        '''
        if a <= 0 or b <= 0:
            raise ValueError('Both edges have to be positive numbers')
        self._a = a
        self._b = b
        vertices = [
            Point(0, 0), Point(self._a, 0),
            Point(self._a, self._b), Point(0, self._b)
        ]
        super().__init__(vertices)

    def __repr__(self):
        return f'Rectangle {self._a} by {self._b}'


class Square(Rectangle):
    '''A polygon with four equal edges and four right angles.

    See docs for `Polygon` for a list of methods and attributes.
    '''
    def __init__(self, a):
        '''Create square with side length equal to `a`

        Parameters
        ----------
        a : number
            Length of the square's edge
        '''
        super().__init__(a, a)

    def __repr__(self):
        return f'Square with edge length {self._a}'


class RegularPolygon(Polygon):
    '''Polygon with all edges and angles equal

    See docs for `Polygon` for a list of methods and attributes.
    '''
    def __init__(self, num_vertices, circumcircle_radius):
        '''Creates new regular polygon.

        Parameters
        ----------
        num_vertices : int
            Number of vertices in the resulting polygon
        circumcircle_radius : number
            Radius of the circle on which circumference vertices lie

        Raises
        ------
        TypeError
            If num_vertices is not an integer
        ValueError
            If circumcircle_radius is <= 0
        '''
        if not isinstance(num_vertices, int):
            raise TypeError('Number of vertices needs to be of type "int"')
        if circumcircle_radius <= 0:
            raise ValueError('Circumcircle radius have to be greater than zero')
        self._circumcircle_radius = circumcircle_radius
        self._inscribed_radius = self._circumcircle_radius * cos(pi / num_vertices)
        angle = 2*pi / num_vertices
        vertices = [
            Point(
                self._circumcircle_radius * sin(i*angle),
                self._circumcircle_radius * cos(i*angle)
            )
            for i in range(num_vertices)]
        super().__init__(vertices)

    def __repr__(self):
        return f'Regular polygon {len(self._vertices)}-vertices'

    @property
    def circumcircle_radius(self):
        '''Radius of circle which contains '''
        return self._circumcircle_radius

    @property
    def inscribed_radius(self):
        '''Radius of circle fully contained within the polygon'''
        return self._inscribed_radius


class Triangle(Polygon):
    '''Three-segment polygon
    
    See docs for `Polygon` for a list of methods and attributes.
    '''
    def __init__(self, p0, p1, p2):
        '''Creates new triangle

        Parameters
        ----------
        p0, p1, p2 : Point
            Instances of basics.Point

        Raises
        ------
        ValueError
            If points passed to constructor do not satisfy triangle inequality
        '''
        edges = [
            Segment(p0, p1),
            Segment(p1, p2),
            Segment(p2, p0),
        ]
        edge_lengths = sorted(map(lambda edge: edge.xy_projection_length, edges), reverse=True)
        if not edge_lengths[0] < edge_lengths[1] + edge_lengths[2]:
            raise ValueError('Provided points do not satisfy triangle inequality')
        super().__init__([p0, p1, p2])

    def __repr__(self):
        return f'Triangle {", ".join(map(str, self.vertices))}'


class Cube(Polyhedron):
    '''Just a cube'''
    def __init__(self, a):
        '''Create cube with edge length equal to `a`

        Parameters
        ----------
        a : number
            edge length
        '''
        super().__init__(Square(a), a)

    def __repr__(self):
        return f'Cube with edge length {self._height}'


class Ellipse(PlaneCurve):
    '''Ellipse, 2D plane curve with equation given by:

    .. math::\frac{x^2}{a^2} + \frac{y^2}{b^2} = 1
    '''
    def __init__(self, major, minor):
        '''Create new ellipe with given major and minor semi-axes.

        Parameters
        ----------
        a : number
            Major semi-axis
        b : number
            Minor semi-axis

        Raises
        ------
        TypeError
            If at least on parameter is not a numbers
        '''
        if not all(isinstance(axis, numbers.Number) for axis in (major, minor)):
            raise TypeError('Major or minor semi-axis is not a number')
        super().__init__()
        self._major = major
        self._minor = minor
        self._eccentricity = sqrt(1 - (self._minor/self._major)**2)

        self._surface_area = self._calculate_surface_area()
        perimeter_func = (
            self._calculate_perimeter
            if _SCIPY_SUPPORT else
            self._approximate_perimeter
        )
        self._perimeter = perimeter_func()

    def __repr__(self):
        return f'Ellipse major {self._major}, minor {self._minor}'

    def _calculate_surface_area(self):
        '''Uses standard formula: pi*(major)*(minor)'''
        return pi*self._major*self._minor

    def _approximate_perimeter(self):
        '''Used only is SciPy is not available.

        The formula taken form
        https://www.johndcook.com/blog/2013/05/05/ramanujan-circumference-ellipse/
        '''
        t = ((self._major - self._minor)/(self._major + self._minor))**2
        return pi*(self._major + self._minor)*(1 + 3*t/(10 + sqrt(4 - 3*t)))

    def _calculate_perimeter(self):
        '''Uses elliptic integral of second kind provided by SciPy'''
        return 4*self._major*ellipe(self._eccentricity**2)

    @property
    def eccentricity(self):
        '''Given by sqrt(1 - (minor/major)**2)'''
        return self._eccentricity

    @property
    def semi_major_axis(self):
        '''Access value for semi-major axis'''
        return self._major

    @property
    def semi_minor_axis(self):
        '''Access value for semi-minor axis'''
        return self._minor

    @property
    def perimeter(self):
        '''Access perimeter (or circumference) of the ellipse'''
        return self._perimeter

    @property
    def surface_area(self):
        '''Access surface area of the ellipse'''
        return self._surface_area
