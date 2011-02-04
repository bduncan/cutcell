#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <boost/multi_array.hpp>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;

class Cell {
    int x;
};

class Grid {
    public:
    typedef boost::multi_array<Nef_polyhedron, 3> V3Nef;
    typedef boost::multi_array<Cell, 3> V3Cell;
    // This segfaults deep in Nef_polyhedron::operator=
    Grid(int X, int Y) {
        N_ = V3Nef(boost::extents[X][Y][5]);
    }
    // This doesn't.
    Grid(int X, int Y, int Z) : N_(boost::extents[X][Y][Z]) {}
    V3Nef N_;
    V3Cell cell_;
};

int main() {
    std::cerr << "Testcase. localN" << std::endl;
    Grid::V3Nef localN = Grid::V3Nef(boost::extents[5][5][5]);
    std::cout << localN[0][0][0];
    std::cerr << "Testcase. initialiser list" << std::endl;
    Grid g(5, 5, 5);
    std::cout << g.N_[0][0][0] << std::endl;
    std::cerr << "Testcase. assignment in constructor" << std::endl;
    Grid f(5, 5);
    std::cout << f.N_[0][0][0] << std::endl;
    return 0;
}
