/* Cartesian cut cell grid generator
 *
 * Copyright 2010 Bruce Duncan, University of Edinburgh
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _CUTCELL_CUTCELL_HPP_
#define _CUTCELL_CUTCELL_HPP_
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_traits_with_normals_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <iostream>
#include <vector>
#include <boost/multi_array.hpp>

namespace cutcell {

// A macro to disallow the copy constructor and operator= functions
// This should be used in the private: declarations for a class
#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
  TypeName(const TypeName&);               \
  void operator=(const TypeName&)

// This is a filtered, Cartesian, quotient kernel.
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

// We need a Polyhedron which can store the plane normals.
typedef CGAL::Polyhedron_traits_with_normals_3<Kernel> Traits;

// The remaining typedefs are just a shorthand...
typedef CGAL::Polyhedron_3<Traits> Polyhedron;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Direction_3 Direction;
typedef Kernel::Aff_transformation_3 Aff_transformation;

// Export some useful CGAL symbols in this namespace.
extern CGAL::Translation TRANSLATION;
extern CGAL::Scaling SCALING;

// The Nef format representation of a unit cube at the origin. Generated from
// make_cube_3.cpp. To regenerate run:
// ./make_cube_3 | awk 'BEGIN {printf "static char cube[] = \""} { printf "%s\\n\\\n", $0 } END { print "\";"}'
static char cube[] = "Selective Nef Complex\n\
standard\n\
vertices   8\n\
halfedges  24\n\
facets     12\n\
volumes    2\n\
shalfedges 48\n\
shalfloops 0\n\
sfaces     16\n\
0 { 0 2, 0 5, 0 1, -2 | 0 1 0 1 } 1\n\
1 { 3 5, 6 11, 2 3, -2 | 0 0 0 1 } 1\n\
2 { 6 8, 12 17, 4 5, -2 | 0 0 1 1 } 1\n\
3 { 9 11, 18 23, 6 7, -2 | 1 0 0 1 } 1\n\
4 { 12 14, 24 29, 8 9, -2 | 1 0 1 1 } 1\n\
5 { 15 17, 30 35, 10 11, -2 | 0 1 1 1 } 1\n\
6 { 18 20, 36 41, 12 13, -2 | 1 1 0 1 } 1\n\
7 { 21 23, 42 47, 14 15, -2 | 1 1 1 1 } 1\n\
0 { 17, 0, 0 0 | 0 0 1 1 } 1\n\
1 { 5, 0, 0 1 | 0 -1 0 1 } 1\n\
2 { 20, 0, 0 3 | 1 0 0 1 } 1\n\
3 { 7, 1, 0 6 | 0 0 1 1 } 1\n\
4 { 10, 1, 0 7 | 1 0 0 1 } 1\n\
5 { 1, 1, 0 9 | 0 1 0 1 } 1\n\
6 { 13, 2, 0 12 | 1 0 0 1 } 1\n\
7 { 3, 2, 0 13 | 0 0 -1 1 } 1\n\
8 { 16, 2, 0 15 | 0 1 0 1 } 1\n\
9 { 18, 3, 0 18 | 0 1 0 1 } 1\n\
10 { 4, 3, 0 19 | -1 0 0 1 } 1\n\
11 { 12, 3, 0 21 | 0 0 1 1 } 1\n\
12 { 11, 4, 0 24 | 0 0 -1 1 } 1\n\
13 { 6, 4, 0 25 | -1 0 0 1 } 1\n\
14 { 23, 4, 0 27 | 0 1 0 1 } 1\n\
15 { 21, 5, 0 30 | 1 0 0 1 } 1\n\
16 { 8, 5, 0 31 | 0 -1 0 1 } 1\n\
17 { 0, 5, 0 33 | 0 0 -1 1 } 1\n\
18 { 9, 6, 0 36 | 0 -1 0 1 } 1\n\
19 { 22, 6, 0 37 | 0 0 1 1 } 1\n\
20 { 2, 6, 0 39 | -1 0 0 1 } 1\n\
21 { 15, 7, 0 42 | -1 0 0 1 } 1\n\
22 { 19, 7, 0 43 | 0 0 -1 1 } 1\n\
23 { 14, 7, 0 45 | 0 -1 0 1 } 1\n\
0 { 1, 6 , , 0 | 0 -1 0 0 } 1\n\
1 { 0, 7 , , 1 | 0 1 0 0 } 1\n\
2 { 3, 8 , , 0 | 0 0 -1 0 } 1\n\
3 { 2, 9 , , 1 | 0 0 1 0 } 1\n\
4 { 5, 4 , , 0 | 0 1 0 -1 } 1\n\
5 { 4, 5 , , 1 | 0 -1 0 1 } 1\n\
6 { 7, 10 , , 0 | -1 0 0 0 } 1\n\
7 { 6, 11 , , 1 | 1 0 0 0 } 1\n\
8 { 9, 16 , , 0 | 0 0 1 -1 } 1\n\
9 { 8, 17 , , 1 | 0 0 -1 1 } 1\n\
10 { 11, 22 , , 0 | 1 0 0 -1 } 1\n\
11 { 10, 23 , , 1 | -1 0 0 1 } 1\n\
0 { 0 } 0\n\
1 { 1 } 1\n\
0 { 1, 4, 2, 0, 1, 32, 10, 6 | 1 0 0 0 } 1\n\
1 { 0, 3, 5, 1, 0, 11, 33, 7 | -1 0 0 0 } 1\n\
2 { 3, 0, 4, 1, 1, 8, 40, 2 | 0 0 1 0 } 1\n\
3 { 2, 5, 1, 2, 0, 41, 9, 3 | 0 0 -1 0 } 1\n\
4 { 5, 2, 0, 2, 1, 38, 34, 4 | 0 -1 0 0 } 1\n\
5 { 4, 1, 3, 0, 0, 35, 39, 5 | 0 1 0 0 } 1\n\
6 { 7, 10, 8, 3, 3, 12, 20, 0 | 0 1 0 0 } 1\n\
7 { 6, 9, 11, 4, 2, 21, 13, 1 | 0 -1 0 0 } 1\n\
8 { 9, 6, 10, 4, 3, 18, 2, 2 | 0 0 1 0 } 1\n\
9 { 8, 11, 7, 5, 2, 3, 19, 3 | 0 0 -1 0 } 1\n\
10 { 11, 8, 6, 5, 3, 0, 14, 6 | 1 0 0 0 } 1\n\
11 { 10, 7, 9, 3, 2, 15, 1, 7 | -1 0 0 0 } 1\n\
12 { 13, 16, 14, 6, 5, 24, 6, 0 | 0 1 0 0 } 1\n\
13 { 12, 15, 17, 7, 4, 7, 25, 1 | 0 -1 0 0 } 1\n\
14 { 15, 12, 16, 7, 5, 10, 32, 6 | 1 0 0 0 } 1\n\
15 { 14, 17, 13, 8, 4, 33, 11, 7 | -1 0 0 0 } 1\n\
16 { 17, 14, 12, 8, 5, 30, 26, 8 | 0 0 -1 0 } 1\n\
17 { 16, 13, 15, 6, 4, 27, 31, 9 | 0 0 1 0 } 1\n\
18 { 19, 22, 20, 9, 7, 40, 8, 2 | 0 0 1 0 } 1\n\
19 { 18, 21, 23, 10, 6, 9, 41, 3 | 0 0 -1 0 } 1\n\
20 { 21, 18, 22, 10, 7, 6, 24, 0 | 0 1 0 0 } 1\n\
21 { 20, 23, 19, 11, 6, 25, 7, 1 | 0 -1 0 0 } 1\n\
22 { 23, 20, 18, 11, 7, 28, 36, 10 | -1 0 0 0 } 1\n\
23 { 22, 19, 21, 9, 6, 37, 29, 11 | 1 0 0 0 } 1\n\
24 { 25, 28, 26, 12, 9, 20, 12, 0 | 0 1 0 0 } 1\n\
25 { 24, 27, 29, 13, 8, 13, 21, 1 | 0 -1 0 0 } 1\n\
26 { 27, 24, 28, 13, 9, 16, 46, 8 | 0 0 -1 0 } 1\n\
27 { 26, 29, 25, 14, 8, 47, 17, 9 | 0 0 1 0 } 1\n\
28 { 29, 26, 24, 14, 9, 44, 22, 10 | -1 0 0 0 } 1\n\
29 { 28, 25, 27, 12, 8, 23, 45, 11 | 1 0 0 0 } 1\n\
30 { 31, 34, 32, 15, 11, 46, 16, 8 | 0 0 -1 0 } 1\n\
31 { 30, 33, 35, 16, 10, 17, 47, 9 | 0 0 1 0 } 1\n\
32 { 33, 30, 34, 16, 11, 14, 0, 6 | 1 0 0 0 } 1\n\
33 { 32, 35, 31, 17, 10, 1, 15, 7 | -1 0 0 0 } 1\n\
34 { 35, 32, 30, 17, 11, 4, 42, 4 | 0 -1 0 0 } 1\n\
35 { 34, 31, 33, 15, 10, 43, 5, 5 | 0 1 0 0 } 1\n\
36 { 37, 40, 38, 18, 13, 22, 44, 10 | -1 0 0 0 } 1\n\
37 { 36, 39, 41, 19, 12, 45, 23, 11 | 1 0 0 0 } 1\n\
38 { 39, 36, 40, 19, 13, 42, 4, 4 | 0 -1 0 0 } 1\n\
39 { 38, 41, 37, 20, 12, 5, 43, 5 | 0 1 0 0 } 1\n\
40 { 41, 38, 36, 20, 13, 2, 18, 2 | 0 0 1 0 } 1\n\
41 { 40, 37, 39, 18, 12, 19, 3, 3 | 0 0 -1 0 } 1\n\
42 { 43, 46, 44, 21, 15, 34, 38, 4 | 0 -1 0 0 } 1\n\
43 { 42, 45, 47, 22, 14, 39, 35, 5 | 0 1 0 0 } 1\n\
44 { 45, 42, 46, 22, 15, 36, 28, 10 | -1 0 0 0 } 1\n\
45 { 44, 47, 43, 23, 14, 29, 37, 11 | 1 0 0 0 } 1\n\
46 { 47, 44, 42, 23, 15, 26, 30, 8 | 0 0 -1 0 } 1\n\
47 { 46, 43, 45, 21, 14, 31, 27, 9 | 0 0 1 0 } 1\n\
0 { 0, 1 , , , 0 } 0\n\
1 { 0, 0 , , , 1 } 1\n\
2 { 1, 7 , , , 0 } 0\n\
3 { 1, 6 , , , 1 } 1\n\
4 { 2, 13 , , , 0 } 0\n\
5 { 2, 12 , , , 1 } 1\n\
6 { 3, 19 , , , 0 } 0\n\
7 { 3, 18 , , , 1 } 1\n\
8 { 4, 25 , , , 0 } 0\n\
9 { 4, 24 , , , 1 } 1\n\
10 { 5, 31 , , , 0 } 0\n\
11 { 5, 30 , , , 1 } 1\n\
12 { 6, 37 , , , 0 } 0\n\
13 { 6, 36 , , , 1 } 1\n\
14 { 7, 43 , , , 0 } 0\n\
15 { 7, 42 , , , 1 } 1\n\
/* end Selective Nef complex */\n\
";

// A description of the properties of one of the faces of a cell.
class Face {
    public:
    Face() {}
    Face(Kernel::Plane_3 normal, Kernel::Point_3 centroid, Kernel::FT area, bool fluid) : normal_(normal), centroid_(centroid), area_(area), fluid_(fluid) {}

    Kernel::Plane_3 const& normal() const { return normal_; }
    Kernel::Point_3 const& centroid() const { return centroid_; }
    Kernel::FT const& area() const { return area_; }
    bool fluid() const { return fluid_; }

    void normal(Kernel::Plane_3 const& normal) { normal_ = normal; }
    void centroid(Kernel::Point_3 const& centroid) { centroid_ = centroid; }
    void area(Kernel::FT const& area) { area_ = area; }
    void fluid(bool fluid) { fluid_ = fluid; }

    private:
    Kernel::Plane_3 normal_;
    Kernel::Point_3 centroid_;
    Kernel::FT area_;
    bool fluid_; // whether our neighbour is a fluid cell.
};

// Whether the cell is a simple solid or fluid cell, or whether it has been cut
// by a solid object.
enum Type { Solid, Fluid, Cut };

// The ijk index of a cell in the Cartesian space.
class Index_3 {
    public:
    Index_3() : i_(0), j_(0), k_(0) {}
    Index_3(int i, int j, int k) : i_(i), j_(j), k_(k) {}
    int i() const { return i_; }
    int j() const { return j_; }
    int k() const { return k_; }

    private:
    int i_, j_, k_;
};

// A description of the cell properties of one of the cut cells.
class Cell {
    public:
    Cell() {}
    Type const& type() const { return type_; }
    Kernel::Point_3 const& centroid() const { return centroid_; }
    Kernel::FT const& volume() const { return volume_; }
    Index_3 const& parent() const { return parent_; }
    std::vector<Face> const& faces() const { return faces_; }

    void type(Type type) { type_ = type; }
    void centroid(Kernel::Point_3 const& centroid) { centroid_ = centroid; }
    void volume(Kernel::FT const& volume) { volume_ = volume; }
    void parent(Index_3 const& parent) { parent_ = parent; }
    void faces(std::vector<Face> const& faces) { faces_ = faces; }
    void addFace(Face const& f) { faces_.push_back(f); }

    private:
    // Cell properties
    Type type_;
    Kernel::Point_3 centroid_;
    Kernel::FT volume_;
    Index_3 parent_; // I, J, K of parent cell.

    // Face properties
    std::vector<Face> faces_;
};

class Grid {
    friend std::ostream& operator<<(std::ostream&, Grid const &);
    friend class GridFormat;

    public:
    typedef boost::multi_array<Nef_polyhedron, 3> V3Nef;
    typedef boost::multi_array<Cell, 3> V3Cell;
    Grid(int, int, int);
    void cut(Nef_polyhedron const&);
    int output_cgns_file(std::string const &) const;
    V3Nef const& grid() const { return N_; }
    V3Cell const& cell() const { return cell_; }

    private:
    std::ostream& output_vrml(std::ostream&) const;
    std::ostream& output_nef(std::ostream&) const;
    std::ostream& output_cgns(std::ostream&) const;
    static int alloc() { return alloc_; }
    typedef V3Nef::index V3NefIndex;
    typedef V3Cell::index V3CellIndex;
    V3Nef N_;
    V3Cell cell_;
    static const int alloc_;
    DISALLOW_COPY_AND_ASSIGN(Grid);
};

class GridFormat {
    // This is a stream modifier for Grid objects. The stream stores the state
    // in a table which is indexed by the Grid::alloc variable above. This
    // index is allocated at compile time. The state is then retrieved in the
    // operator<< of the Grid class, which then calls the appropriate function
    // on the Grid object.
    // To add a new output format, add an integer constant here and add a case
    // to the switch statement in operator<<, below.
    public:
    static const int OUTPUT_VRML = 1;
    static const int OUTPUT_NEF = 2;
    static const int OUTPUT_CGNS = 3;
    explicit GridFormat(int format) : format(format) {}
    friend std::ostream & operator<<(std::ostream& os, GridFormat const& gf) {
        os.iword(gridAlloc()) = gf.format;
        return os;
    }
    private:
    int format;
    // operator<< needs Grid::alloc, but it can't get it as a friend. Pass it on.
    static int gridAlloc() { return Grid::alloc(); }
    DISALLOW_COPY_AND_ASSIGN(GridFormat);
};

} // namespace cutcell
#endif // _CUTCELL_CUTCELL_HPP_
