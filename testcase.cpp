#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <boost/multi_array.hpp>
#include <iostream>
#include <sstream>
#include <cutcell.hpp>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;

int main(int argc, char **argv) {
    std::istringstream in(cutcell::cube);
    Nef_polyhedron N1;
    in >> N1;
    Nef_polyhedron::Halffacet_const_iterator fi;
    int i = 0;
    CGAL_forall_facets(fi, N1) {
        std::cout << "Halffacet " << i++ << " plane: " << fi->plane().a()
                                         << " "  << fi->plane().b()
                                         << " "  << fi->plane().c()
                                         << " "  << fi->plane().d()
                                         << std::endl;
        Nef_polyhedron::Halffacet_cycle_const_iterator fc;
        fc = fi->facet_cycles_begin();
        int j = 0;
        CGAL_For_all(fc, fi->facet_cycles_end()) {
            Nef_polyhedron::SHalfedge_const_handle se = Nef_polyhedron::SHalfedge_const_handle(fc);
            Nef_polyhedron::SHalfedge_around_facet_const_circulator hc_start(se);
            Nef_polyhedron::SHalfedge_around_facet_const_circulator hc_end(hc_start);
            CGAL_For_all(hc_start, hc_end)
                std::cout << "Point " << j++ << ": " << hc_start->source()->center_vertex()->point() << std::endl;
        }
    }
    return 0;
}
