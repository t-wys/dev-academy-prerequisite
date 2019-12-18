'''Entry point for the script'''
import argparse
from basics import Point
from geomobjects import (Rectangle, Square, Triangle, Ellipse, Cube, RightPyramid)


def main():
    '''Run this function'''
    parser = argparse.ArgumentParser(
        description='Calculates basic parameters of geometric objects (perimeter/surface area/volume).'
    )

    subparsers = parser.add_subparsers(
        dest='geom_object',
        help='For more options type chosen object name followed by -h/--help.'
    )

    parser_tri = subparsers.add_parser('triangle')
    parser_tri.add_argument(
        'point', nargs=3, metavar='POINT', type=str,
        help='Point 2D coordinates in form: {X},{Y} where X, Y are floats, e.g.: 1.2,5 for point x=1.2 and y=5.0'
    )

    parser_rect = subparsers.add_parser('rectangle')
    parser_rect.add_argument(
        'edge', nargs=2, metavar='EDGE', type=float,
        help='Lenght of rectangle\'s edge, floating point number'
    )

    parser_sqr = subparsers.add_parser('square')
    parser_sqr.add_argument(
        'edge', nargs=1, metavar='EDGE', type=float,
        help='Length of square\'s edge, floating point number'
    )

    parser_ell = subparsers.add_parser('ellipse')
    parser_ell.add_argument(
        'major', metavar='MAJOR', type=float,
        help='Length of semi-major axis, floating point number'
    )
    parser_ell.add_argument(
        'minor', metavar='MINOR', type=float,
        help='Length of semi-minor axis. floating point number'
    )

    parser_cube = subparsers.add_parser('cube')
    parser_cube.add_argument(
        'edge', metavar='EDGE', type=float,
        help='Length of cube\'s edge, floating point number'
    )

    parser_pyr = subparsers.add_parser('right-pyramid')
    parser_pyr.add_argument(
        'num_vertices', metavar='NUM-VERTICES', type=int,
        help='Number of vertices in base, integer'
    )
    parser_pyr.add_argument(
        'radius', metavar='RADIUS', type=float,
        help='Radius of circumcircle, floating point number'
    )
    parser_pyr.add_argument(
        'height', metavar='HEIGHT', type=float,
        help='Height of pyramid, floating point number'
    )

    args = parser.parse_args()

    if args.geom_object is None:
        # required=True for add_subparsers added in Python 3.7
        # so let's just print help if there are no arguments
        parser.print_help()
        return

    response = 'Requested object: {geom_obj}\n'
    response_2d = 'perimeter: {perimeter:.2g}, surface area: {surf_area:.2g}'
    response_3d = 'surface area: {surf_area:.2g}, volume: {volume:.2g}'

    if args.geom_object == 'square':
        square = Square(args.edge[0])
        response = response.format(geom_obj=square)
        response += response_2d.format(
            perimeter=square.perimeter, surf_area=square.surface_area
        )

    if args.geom_object == 'rectangle':
        rectangle = Rectangle(*args.edge)
        response = response.format(geom_obj=rectangle)
        response += response_2d.format(
            perimeter=rectangle.perimeter, surf_area=rectangle.surface_area
        )

    if args.geom_object == 'triangle':
        triangle = Triangle(*[Point.from_str(pnt) for pnt in args.point])
        response = response.format(geom_obj=triangle)
        response += response_2d.format(
            perimeter=triangle.perimeter, surf_area=triangle.surface_area
        )

    if args.geom_object == 'ellipse':
        ellipse = Ellipse(args.major, args.minor)
        response = response.format(geom_obj=ellipse)
        response += response_2d.format(
            perimeter=ellipse.perimeter, surf_area=ellipse.surface_area
        )

    if args.geom_object == 'cube':
        cube = Cube(args.edge)
        response = response.format(geom_obj=cube)
        response += response_3d.format(
            surf_area=cube.surface_area, volume=cube.volume
        )

    if args.geom_object == 'right-pyramid':
        pyramid = RightPyramid.from_base_parameters(
            args.num_vertices, args.radius,
            args.height
        )
        response = response.format(geom_obj=pyramid)
        response += response_3d.format(
            surf_area=pyramid.surface_area, volume=pyramid.volume
        )

    print(response)


if __name__ == "__main__":
    main()
