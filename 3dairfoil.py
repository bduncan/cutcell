#!/usr/bin/python
# vim: set fileencoding=utf-8:

"""A program to extend the dimensions of an OFF file.

XFoil can be used to create a 2D section of an airfoil, for example NACA
airfoils. This program intends to take the 2D section as a set of points and
extend it into the third dimension. It also creates facets which join the
points, so as to create a complete OFF description of the airfoil.

The input should be any number of x-y points. The output will include an OFF
header, followed by one copy of the points at z=0 and a copy at z=1. The facets
consist of a section of rectangular elements which join pairs of points on each
of the copies to the opposite set of points. A second set of facets joins pairs
of points on a single section to a point at the centroid of all points in that
set. These are triangles.

"""

import sys
from contextlib import contextmanager


@contextmanager
def opendash(args, i=1, mode='r'):
    """Open a file.

    Open the file given at the index of a list, or stdin if that index doesn't
    exist or is the string '-'.

    """

    if i >= len(args) or args[i] == '-':
        yield sys.stdin
    else:
        yield open(args[i], mode)


def main():
    # Read the x,y points from the named file or standard input.
    points = []
    with opendash(sys.argv) as f:
        for line in f:
            points += [[float(x) for x in line.split()], ]

    np = len(points)
    with sys.stdout as f:

        # Output the header
        print >> f, "OFF"
        print >> f, 2 * np + 2, 3 * np, 0

        # Output the 2D airfoil at z=0 and z=1
        for point in points:
            print >> f, point[0], point[1], 0
        for point in points:
            print >> f, point[0], point[1], 1

        # Put a point at the centre of each face in order to join the edges as
        # triangles
        print >> f, sum(p[0] for p in points) / np, 0, 0
        print >> f, sum(p[0] for p in points) / np, 0, 1

        # Join pairs of points on each face into cuboids.
        for i in range(np):
            print >> f, 4, i, (i + 1) % np, np + (i + 1) % np, np + i

        # Join pairs of points on each face with the centroid to make
        # triangles.
        for i in range(np):
            print >> f, 3, (i + 1) % np, i, 2 * np
        for i in range(np):
            print >> f, 3, np + i, np + ((i + 1) % np), 2 * np + 1

if __name__ == '__main__':
    main()
