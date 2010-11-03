#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Extended_homogeneous.h>
#include <CGAL/Quotient.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpzf.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_traits_with_normals_3.h>
#include <CGAL/Polyhedron_traits_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/IO/Polyhedron_VRML_2_ostream.h>
#include <CGAL/centroid.h>
#include <list>
#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
//typedef CGAL::Homogeneous<double> Kernel;
//typedef CGAL::Extended_homogeneous<CGAL::Gmpz> Kernel;
typedef CGAL::Polyhedron_traits_with_normals_3<Kernel> Traits;
//typedef CGAL::Polyhedron_traits_3<Kernel> Traits;
typedef CGAL::Polyhedron_3<Traits> Polyhedron;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef Kernel::Direction_3 Direction;
typedef Kernel::Point_3 Point;

struct Normal_vector {
    template <class Facet>
    typename Facet::Plane_3 operator()( Facet& f) {
        typename Facet::Halfedge_handle h = f.halfedge();
        // Facet::Plane_3 is the normal vector type. We assume the
        // CGAL Kernel here and use its global functions.
        return typename Facet::Plane_3(CGAL::cross_product(
            h->next()->vertex()->point() - h->vertex()->point(),
            h->next()->next()->vertex()->point() - h->next()->vertex()->point()));
    }
};

int main() {
    Point p( 1.0, 0.0, 0.0);
    Point q( 0.0, 1.0, 0.0);
    Point r( 0.0, 0.0, 1.0);
    Point s( 0.0, 0.0, 0.0);
    Polyhedron P;

    P.make_tetrahedron( p, q, r, s);
    Nef_polyhedron N(P);

    Polyhedron P1;
    N.convert_to_polyhedron(P1);
    std::transform(P1.facets_begin(), P1.facets_end(), P1.planes_begin(),
                   Normal_vector());
    Polyhedron::Facet_const_iterator fi = P1.facets_begin();
    std::cerr << "Facet... " << &fi << std::endl;
    std::cerr << "fi->plane(): " << fi->plane() << std::endl;
    Direction d1(fi->plane());
    //++fi;
    std::cout << (d1 == Direction(fi->plane())) << std::endl;

    return 0;
}
