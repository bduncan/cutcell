/* Simple Nef_polyhedron_3
 *
 * This program creates a simple Nef_polyhedron using the intersection method
 * shown in examples/Nef_3/point_set_operations.cpp. The shape generated is a
 * cube which is then output in Nef format.
 *
 * This file was copied from CGAL-3.8/examples/Nef_3/point_set_operations.cpp,
 * which has no copyright information. Quoting CGAL-3.8/LICENSE:
 *
 * All other files that do not have an explicit copyright notice (e.g., all
 * examples and some demos) are licensed under a very permissive license. The
 * exact license text can be found in the file LICENSE.FREE_USE.
 *
 * Therefore this file is licensed under the CGAL open source license, copied below:
 *
 * Copyright (c) 1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007
 * Utrecht University (The Netherlands),
 * ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
 * INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
 * (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
 * and Tel-Aviv University (Israel).  All rights reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 */
#include <CGAL/Gmpz.h>
#include <CGAL/Extended_homogeneous.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <iostream>

typedef CGAL::Extended_homogeneous<CGAL::Gmpz>  Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel>  Nef_polyhedron;
typedef Kernel::Plane_3 Plane_3;

int main() {

    Nef_polyhedron N1(Plane_3( 1, 0, 0,-1));
    Nef_polyhedron N2(Plane_3(-1, 0, 0, 0));
    Nef_polyhedron N3(Plane_3( 0, 1, 0,-1));
    Nef_polyhedron N4(Plane_3( 0,-1, 0, 0));
    Nef_polyhedron N5(Plane_3( 0, 0, 1,-1));
    Nef_polyhedron N6(Plane_3( 0, 0,-1, 0));
    Nef_polyhedron Cube = N1 * N2 * N3 * N4 * N5 * N6;

    // Output the Nef_polyhedron in the Nef file format.
    std::cout << Cube;
    return 0;
}
