#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef Polyhedron::Vertex_iterator Vertex_iterator;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Aff_transformation_3 Aff_transformation;

enum Cell { Solid, Fluid, Cut };

int main() {
    Nef_polyhedron UnitCube;
    const unsigned NX = 5, NY = 5, NZ = 5;

    std::cout << "Making cube" << std::endl;
    std::ifstream in("cube.nef3"); // Load the definition of a unit cube from file
    in >> UnitCube;

    Nef_polyhedron N[NX][NY][NZ];
    Nef_polyhedron N1(UnitCube);
    Cell cell[NX][NY][NZ];
    std::cout << "Making grid" << std::endl;
    for (int x = 0; x < NX; x++)
        for (int y = 0; y < NY; y++)
            for (int z = 0; z < NZ; z++) {
                N[x][y][z] = UnitCube;
                Aff_transformation Aff(CGAL::TRANSLATION, Vector(x,y,z));
                N[x][y][z].transform(Aff);
            }

    Aff_transformation Aff1(CGAL::TRANSLATION, Vector(0.5,0.5,0.5));
    Aff_transformation Aff2(CGAL::SCALING, 1.5);
    std::cout << "Transforming..." << std::endl;
    N1.transform(Aff1);
    N1.transform(Aff2);

    std::cout << "Computing Union. Solid = " << Solid << ", Fluid = " << Fluid << ", Cut = " << Cut << std::endl;
    for (int x = 0; x < NX; x++) {
        for (int y = 0; y < NY; y++)
            for (int z = 0; z < NZ; z++) {
                Nef_polyhedron I = N[x][y][z] - N1;
                if (I.is_empty())
                    cell[x][y][z] = Solid;
                else if (I == N[x][y][z])
                    cell[x][y][z] = Fluid;
                else
                    cell[x][y][z] = Cut;
                std::cout << cell[x][y][z];
            }
        std::cout << std::endl;
    }

    if (N1.is_simple()) {
        Polyhedron P;
        N1.convert_to_polyhedron(P);
        //std::cout << N;
        CGAL::set_ascii_mode( std::cout);
        for (Vertex_iterator v = P.vertices_begin(); v != P.vertices_end(); ++v)
            std::cout << v->point() << std::endl;
        return 0;
    } else {
        std::cerr << "N is not a 2-manifold" << std::endl;
        return 1;
    }
}
