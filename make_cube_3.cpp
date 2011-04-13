/* Simple Nef_olyhedron
 *
 * This program creates a simple Nef_polyhedron using the intersection method
 * shown in examples/Nef_3/point_set_operations.cpp. The shape generated is a
 * cube which is then output in Nef format.
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
