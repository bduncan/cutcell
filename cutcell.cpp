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


int main() {
    Polyhedron P;

    make_cube_3(P);

    Nef_polyhedron N(P);
    Aff_transformation Aff(CGAL::TRANSLATION, Vector(0.5, 0.5, 0.5));
    N.transform(Aff);
    N.transform(Aff);
    N.transform(Aff);
    if (N.is_simple()) {
        N.convert_to_polyhedron(P);
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
