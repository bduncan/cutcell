#define CGAL_IA_NO_X86_OVER_UNDER_FLOW_PROTECT
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/IO/Polyhedron_VRML_2_ostream.h>
#include <CGAL/centroid.h>
#include <list>
#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef Polyhedron::Vertex_iterator Vertex_iterator;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Direction_3 Direction;
typedef Kernel::Aff_transformation_3 Aff_transformation;
typedef Kernel::Plane_3 Plane_3;
typedef Kernel::Vector_3 Vector_3;
typedef Kernel::Point_3 Point_3;

static char cube[] = "Selective Nef Complex\n\
standard\n\
vertices   8\n\
halfedges  24\n\
facets     12\n\
volumes    2\n\
shalfedges 48\n\
shalfloops 0\n\
sfaces     16\n\
0 { 0 2, 0 5, 0 1, -2 | -1 -1 -1 1 } 1\n\
1 { 3 5, 6 11, 2 3, -2 | -1 -1 1 1 } 1\n\
2 { 6 8, 12 17, 4 5, -2 | -1 1 -1 1 } 1\n\
3 { 9 11, 18 23, 6 7, -2 | -1 1 1 1 } 1\n\
4 { 12 14, 24 29, 8 9, -2 | 1 -1 -1 1 } 1\n\
5 { 15 17, 30 35, 10 11, -2 | 1 -1 1 1 } 1\n\
6 { 18 20, 36 41, 12 13, -2 | 1 1 -1 1 } 1\n\
7 { 21 23, 42 47, 14 15, -2 | 1 1 1 1 } 1\n\
0 { 3, 0, 0 0 | 0 0 1 1 } 1\n\
1 { 6, 0, 0 1 | 0 1 0 1 } 1\n\
2 { 12, 0, 0 3 | 1 0 0 1 } 1\n\
3 { 0, 1, 0 6 | 0 0 -1 1 } 1\n\
4 { 9, 1, 0 7 | 0 1 0 1 } 1\n\
5 { 15, 1, 0 9 | 1 0 0 1 } 1\n\
6 { 1, 2, 0 12 | 0 -1 0 1 } 1\n\
7 { 10, 2, 0 13 | 0 0 1 1 } 1\n\
8 { 18, 2, 0 15 | 1 0 0 1 } 1\n\
9 { 4, 3, 0 18 | 0 -1 0 1 } 1\n\
10 { 7, 3, 0 19 | 0 0 -1 1 } 1\n\
11 { 21, 3, 0 21 | 1 0 0 1 } 1\n\
12 { 2, 4, 0 24 | -1 0 0 1 } 1\n\
13 { 16, 4, 0 25 | 0 0 1 1 } 1\n\
14 { 19, 4, 0 27 | 0 1 0 1 } 1\n\
15 { 5, 5, 0 30 | -1 0 0 1 } 1\n\
16 { 13, 5, 0 31 | 0 0 -1 1 } 1\n\
17 { 22, 5, 0 33 | 0 1 0 1 } 1\n\
18 { 8, 6, 0 36 | -1 0 0 1 } 1\n\
19 { 14, 6, 0 37 | 0 -1 0 1 } 1\n\
20 { 23, 6, 0 39 | 0 0 1 1 } 1\n\
21 { 11, 7, 0 42 | -1 0 0 1 } 1\n\
22 { 17, 7, 0 43 | 0 -1 0 1 } 1\n\
23 { 20, 7, 0 45 | 0 0 -1 1 } 1\n\
0 { 11, 1 , , 0 | -1 0 0 -1 } 1\n\
1 { 10, 2 , , 0 | 0 -1 0 -1 } 1\n\
2 { 9, 5 , , 0 | 0 0 -1 -1 } 1\n\
3 { 8, 10 , , 0 | 0 0 1 -1 } 1\n\
4 { 7, 17 , , 0 | 0 1 0 -1 } 1\n\
5 { 6, 28 , , 0 | 1 0 0 -1 } 1\n\
6 { 5, 29 , , 1 | -1 0 0 1 } 1\n\
7 { 4, 16 , , 1 | 0 -1 0 1 } 1\n\
8 { 3, 11 , , 1 | 0 0 -1 1 } 1\n\
9 { 2, 4 , , 1 | 0 0 1 1 } 1\n\
10 { 1, 3 , , 1 | 0 1 0 1 } 1\n\
11 { 0, 0 , , 1 | 1 0 0 1 } 1\n\
0 { 0 } 0\n\
1 { 1 } 1\n\
0 { 1, 3, 4, 0, 0, 7, 12, 11 | -1 0 0 0 } 1\n\
1 { 0, 5, 2, 1, 1, 13, 6, 0 | 1 0 0 0 } 1\n\
2 { 3, 1, 5, 0, 1, 9, 24, 1 | 0 1 0 0 } 1\n\
3 { 2, 4, 0, 2, 0, 25, 8, 10 | 0 -1 0 0 } 1\n\
4 { 5, 0, 3, 1, 0, 15, 26, 9 | 0 0 -1 0 } 1\n\
5 { 4, 2, 1, 2, 1, 27, 14, 2 | 0 0 1 0 } 1\n\
6 { 7, 9, 10, 3, 3, 1, 18, 0 | 1 0 0 0 } 1\n\
7 { 6, 11, 8, 4, 2, 19, 0, 11 | -1 0 0 0 } 1\n\
8 { 9, 7, 11, 3, 2, 3, 30, 10 | 0 -1 0 0 } 1\n\
9 { 8, 10, 6, 5, 3, 31, 2, 1 | 0 1 0 0 } 1\n\
10 { 11, 6, 9, 4, 3, 21, 32, 3 | 0 0 -1 0 } 1\n\
11 { 10, 8, 7, 5, 2, 33, 20, 8 | 0 0 1 0 } 1\n\
12 { 13, 15, 16, 6, 4, 0, 19, 11 | -1 0 0 0 } 1\n\
13 { 12, 17, 14, 7, 5, 18, 1, 0 | 1 0 0 0 } 1\n\
14 { 15, 13, 17, 6, 5, 5, 36, 2 | 0 0 1 0 } 1\n\
15 { 14, 16, 12, 8, 4, 37, 4, 9 | 0 0 -1 0 } 1\n\
16 { 17, 12, 15, 7, 4, 23, 38, 7 | 0 1 0 0 } 1\n\
17 { 16, 14, 13, 8, 5, 39, 22, 4 | 0 -1 0 0 } 1\n\
18 { 19, 21, 22, 9, 7, 6, 13, 0 | 1 0 0 0 } 1\n\
19 { 18, 23, 20, 10, 6, 12, 7, 11 | -1 0 0 0 } 1\n\
20 { 21, 19, 23, 9, 6, 11, 42, 8 | 0 0 1 0 } 1\n\
21 { 20, 22, 18, 11, 7, 43, 10, 3 | 0 0 -1 0 } 1\n\
22 { 23, 18, 21, 10, 7, 17, 44, 4 | 0 -1 0 0 } 1\n\
23 { 22, 20, 19, 11, 6, 45, 16, 7 | 0 1 0 0 } 1\n\
24 { 25, 27, 28, 12, 8, 2, 31, 1 | 0 1 0 0 } 1\n\
25 { 24, 29, 26, 13, 9, 30, 3, 10 | 0 -1 0 0 } 1\n\
26 { 27, 25, 29, 12, 9, 4, 37, 9 | 0 0 -1 0 } 1\n\
27 { 26, 28, 24, 14, 8, 36, 5, 2 | 0 0 1 0 } 1\n\
28 { 29, 24, 27, 13, 8, 35, 40, 5 | -1 0 0 0 } 1\n\
29 { 28, 26, 25, 14, 9, 41, 34, 6 | 1 0 0 0 } 1\n\
30 { 31, 33, 34, 15, 11, 8, 25, 10 | 0 -1 0 0 } 1\n\
31 { 30, 35, 32, 16, 10, 24, 9, 1 | 0 1 0 0 } 1\n\
32 { 33, 31, 35, 15, 10, 10, 43, 3 | 0 0 -1 0 } 1\n\
33 { 32, 34, 30, 17, 11, 42, 11, 8 | 0 0 1 0 } 1\n\
34 { 35, 30, 33, 16, 11, 29, 46, 6 | 1 0 0 0 } 1\n\
35 { 34, 32, 31, 17, 10, 47, 28, 5 | -1 0 0 0 } 1\n\
36 { 37, 39, 40, 18, 12, 14, 27, 2 | 0 0 1 0 } 1\n\
37 { 36, 41, 38, 19, 13, 26, 15, 9 | 0 0 -1 0 } 1\n\
38 { 39, 37, 41, 18, 13, 16, 45, 7 | 0 1 0 0 } 1\n\
39 { 38, 40, 36, 20, 12, 44, 17, 4 | 0 -1 0 0 } 1\n\
40 { 41, 36, 39, 19, 12, 28, 47, 5 | -1 0 0 0 } 1\n\
41 { 40, 38, 37, 20, 13, 46, 29, 6 | 1 0 0 0 } 1\n\
42 { 43, 45, 46, 21, 15, 20, 33, 8 | 0 0 1 0 } 1\n\
43 { 42, 47, 44, 22, 14, 32, 21, 3 | 0 0 -1 0 } 1\n\
44 { 45, 43, 47, 21, 14, 22, 39, 4 | 0 -1 0 0 } 1\n\
45 { 44, 46, 42, 23, 15, 38, 23, 7 | 0 1 0 0 } 1\n\
46 { 47, 42, 45, 22, 15, 34, 41, 6 | 1 0 0 0 } 1\n\
47 { 46, 44, 43, 23, 14, 40, 35, 5 | -1 0 0 0 } 1\n\
0 { 0, 0 , , , 0 } 0\n\
1 { 0, 2 , , , 1 } 1\n\
2 { 1, 8 , , , 0 } 0\n\
3 { 1, 6 , , , 1 } 1\n\
4 { 2, 12 , , , 0 } 0\n\
5 { 2, 14 , , , 1 } 1\n\
6 { 3, 20 , , , 0 } 0\n\
7 { 3, 18 , , , 1 } 1\n\
8 { 4, 24 , , , 1 } 1\n\
9 { 4, 26 , , , 0 } 0\n\
10 { 5, 32 , , , 1 } 1\n\
11 { 5, 30 , , , 0 } 0\n\
12 { 6, 36 , , , 1 } 1\n\
13 { 6, 38 , , , 0 } 0\n\
14 { 7, 44 , , , 1 } 1\n\
15 { 7, 42 , , , 0 } 0\n\
/* end Selective Nef complex */\n";

int main() {
    Nef_polyhedron UnitCube;//(Plane_3(0.5,0.5,0.5,1));
    Direction d(Vector(1, 1, 1));
    std::cerr << "Making cube" << std::endl;
    std::istringstream in(cube); // Load the definition of a unit cube from string
    in >> UnitCube;
    Plane_3 p(0.5, 0.5, 0.5, 0.5);

    Polyhedron P;
    UnitCube.convert_to_polyhedron(P);
    for (Polyhedron::Facet_const_iterator fi = P.facets_begin(); fi != P.facets_end(); ++fi) {
        std::cerr << "Facet... " << &fi << std::endl;
        std::cerr << "fi->plane(): " << fi->plane() << std::endl;
        std::cerr << "p:           " << p << std::endl;
        std::cerr << "fi->plane().orthogonal_direction(): " << fi->plane().orthogonal_direction() << std::endl;
        Direction dead = fi->plane().orthogonal_direction();
        std::cerr << "fi->plane().orthogonal_direction() dx dy dz: " << dead.dx() << ", " << dead.dy() << ", " << dead.dz() << std::endl;
        std::cerr << "p.orthogonal_direction():           " << p.orthogonal_direction()<< std::endl;
        std::cerr << "d: " << d << std::endl;
        std::cerr << "p.orthogonal_direction() == d:  " << (p.orthogonal_direction() == d) << std::endl;
        std::cerr << "fi->plane().orthogonal_direction() == d: " << (fi->plane().orthogonal_direction() == d) << std::endl;
    }

    return 0;
}
