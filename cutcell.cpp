#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <iostream>

typedef CGAL::Filtered_kernel<
        CGAL::Simple_cartesian<
        CGAL::Quotient<
        CGAL::MP_Float> > > Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef Polyhedron::Vertex_iterator Vertex_iterator;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Aff_transformation_3 Aff_transformation;

template <class Poly>
typename Poly::Halfedge_handle make_cube_3( Poly& P) {
    // appends a cube of size [0,1]^3 to the polyhedron P.
    CGAL_precondition( P.is_valid());
    typedef typename Poly::Point_3         Point;
    typedef typename Poly::Halfedge_handle Halfedge_handle;
    Halfedge_handle h = P.make_tetrahedron( Point( 1, 0, 0),
                                            Point( 0, 0, 1),
                                            Point( 0, 0, 0),
                                            Point( 0, 1, 0));
    Halfedge_handle g = h->next()->opposite()->next();             // Fig. (a)
    P.split_edge( h->next());
    P.split_edge( g->next());
    P.split_edge( g);                                              // Fig. (b)
    h->next()->vertex()->point()     = Point( 1, 0, 1);
    g->next()->vertex()->point()     = Point( 0, 1, 1);
    g->opposite()->vertex()->point() = Point( 1, 1, 0);            // Fig. (c)
    Halfedge_handle f = P.split_facet( g->next(),
                                       g->next()->next()->next()); // Fig. (d)
    Halfedge_handle e = P.split_edge( f);
    e->vertex()->point() = Point( 1, 1, 1);                        // Fig. (e)
    P.split_facet( e, f->next()->next());                          // Fig. (f)
    CGAL_postcondition( P.is_valid());
    return h;
}

enum Cell { Solid, Fluid, Cut };

int main() {
    Polyhedron P;

    std::cout << "Making cube" << std::endl;
    make_cube_3(P);

    Nef_polyhedron N[10][10][10];
    Nef_polyhedron N1;
    Cell cell[10][10][10];
    std::cout << "Making grid" << std::endl;
    for (int x = 0; x < 10; x++)
        for (int y = 0; y < 10; y++)
            for (int z = 0; z < 10; z++) {
                N[x][y][z] = P;
                Aff_transformation Aff(CGAL::TRANSLATION, Vector(x,y,z));
                N[x][y][z].transform(Aff);
            }

    N1.convert_to_polyhedron(P);
    //std::cout << N;
    CGAL::set_ascii_mode( std::cout);
    for (Vertex_iterator v = P.vertices_begin(); v != P.vertices_end(); ++v)
        std::cout << v->point() << std::endl;
    Aff_transformation Aff(2,0,0,0.5,
                           0,2,0,0.5,
                           0,0,2,0.5);
    std::cout << "Transforming..." << std::endl;
    N1.transform(Aff);
    N1.convert_to_polyhedron(P);

    std::cout << "Computing Union" << std::endl;
    for (int x = 0; x < 10; x++)
        for (int y = 0; y < 10; y++)
            for (int z = 0; z < 10; z++) {
                Nef_polyhedron I = N1 * N[x][y][z];
                if (z == 1) {
                    I.convert_to_polyhedron(P);
                }
                if (I.is_empty())
                    cell[x][y][z] = Fluid;
                else if (I == N[x][y][z])
                    cell[x][y][z] = Solid;
                else
                    cell[x][y][z] = Cut;
                std::cout << cell[x][y][z];
            }

    if (N1.is_simple()) {
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
