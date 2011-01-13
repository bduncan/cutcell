/* Cartesian cut cell grid generator
 *
 * Generate a 3D Cartesian array of unit cubes, with the given solid
 * subtracted from them. Export various properties of the cubes.
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

#include <cgnslib.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_traits_with_normals_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/IO/Polyhedron_VRML_2_ostream.h>
#include <CGAL/centroid.h>
#include <CGAL/Triangulation_3.h>
#include <list>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>

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
typedef CGAL::Triangulation_3<Kernel> Triangulation;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Direction_3 Direction;
typedef Kernel::Aff_transformation_3 Aff_transformation;

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

// A function to compute the plane normal of a Facet
// From CGAL-3.7/examples/Polyhedron/polyhedron_prog_normals.cpp
struct Normal_vector {
    template <class Facet>
    typename Facet::Plane_3 operator()(Facet const& f) {
        typename Facet::Halfedge_const_handle h = f.halfedge();
        // Facet::Plane_3 is the normal vector type. We assume the
        // CGAL Kernel here and use its global functions.
        return typename Facet::Plane_3(CGAL::cross_product(
            h->next()->vertex()->point() - h->vertex()->point(),
            h->next()->next()->vertex()->point() - h->next()->vertex()->point()));
    }
};

// A function to compare two Directions, allowing them to be stored in a map.
class compare {
    public:
    bool operator()(Direction const& x, Direction const& y) {
        return x != y;
    }
};

// A description of the properties of one of the faces of a cell.
class Face {
    public:
    Face() {}
    Face(Vector normal, Kernel::Point_3 centroid, Kernel::FT area, bool fluid) : normal_(normal), centroid_(centroid), area_(area), fluid_(fluid) {}

    Vector const& normal() const { return normal_; }
    Kernel::Point_3 const& centroid() const { return centroid_; }
    Kernel::FT const& area() const { return area_; }
    bool fluid() const { return fluid_; }

    void normal(Vector const& normal) { normal_ = normal; }
    void centroid(Kernel::Point_3 const& centroid) { centroid_ = centroid; }
    void area(Kernel::FT const& area) { area_ = area; }
    void fluid(bool fluid) { fluid_ = fluid; }

    private:
    Vector normal_;
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
    typedef std::map<Direction, std::vector<Face>, compare> faceMap;

    Cell() {}
    Type const& type() const { return type_; }
    Kernel::Point_3 const& centroid() const { return centroid_; }
    Kernel::FT const& volume() const { return volume_; }
    Index_3 const& parent() const { return parent_; }
    faceMap const& faces() const { return faces_; }

    void type(Type type) { type_ = type; }
    void centroid(Kernel::Point_3 const& centroid) { centroid_ = centroid; }
    void volume(Kernel::FT const& volume) { volume_ = volume; }
    void parent(Index_3 const& parent) { parent_ = parent; }
    void faces(faceMap const& faces) { faces_ = faces; }
    void addFace(Direction const& d, Face const& f) { faces_[d].push_back(f); }

    private:
    // Cell properties
    Type type_;
    Kernel::Point_3 centroid_;
    Kernel::FT volume_;
    Index_3 parent_; // I, J, K of parent cell.

    // Face properties
    faceMap faces_;
};

class Grid {
    friend std::ostream& operator<<(std::ostream&, Grid &);
    friend class GridFormat;
    public:
    // Shorthands for 3D vector arrays.
    typedef std::vector< std::vector< std::vector< Nef_polyhedron > > > V3Nef;
    typedef std::vector< std::vector< std::vector< Cell > > > V3Cell;
    Grid(int X, int Y, int Z) {
        // Create the Unit Cube from the global string definition.
        Nef_polyhedron UnitCube;
        std::istringstream in(cube);
        in >> UnitCube;
        assert(UnitCube.number_of_vertices() == 8);
        assert(UnitCube.number_of_facets() == 6);
        assert(UnitCube.number_of_edges() == 12);
        assert(UnitCube.number_of_volumes() == 2);

        // Create an array to hold NX*NY*NZ Nef_polyhedron cubes.
        N_ = std::vector< std::vector< std::vector< Nef_polyhedron> > >(X, std::vector< std::vector< Nef_polyhedron> >(Y, std::vector< Nef_polyhedron >(Z, UnitCube)));

        // Create the array to store the cell (and therefore face) properties at
        // each point.
        cell_ = std::vector< std::vector< std::vector< Cell> > >(X, std::vector< std::vector< Cell> >(Y, std::vector< Cell >(Z)));

        // Copy and translate the cube for each point.
        for (int x = 0; x < N_.size(); ++x)
            for (int y = 0; y < N_[x].size(); ++y)
                for (int z = 0; z < N_[x][y].size(); ++z) {
                    Aff_transformation Aff(CGAL::TRANSLATION, Vector(x, y, z));
                    N_[x][y][z].transform(Aff);
                }
    };
    void cut(Nef_polyhedron const& N1) {
        for (int x = 0; x < N_.size(); ++x) {
            for (int y = 0; y < N_[x].size(); ++y)
                for (int z = 0; z < N_[x][y].size(); ++z) {
                    // Compute the intersection of this part of the grid with the
                    // test cube.
                    Nef_polyhedron I = N_[x][y][z] - N1;
                    Polyhedron P;
                    // Convert this new cut Nef_polyhedron I into the Polyhedron P.
                    I.convert_to_polyhedron(P);

                    // Set the type of the new cell.
                    if (I.is_empty())
                        // No points, must be completely inside the solid.
                        cell_[x][y][z].type(Solid);
                    else if (I == N_[x][y][z])
                        // Unchanged, must be completely outside the solid.
                        cell_[x][y][z].type(Fluid);
                    else
                        // Something else, must be a cut cell.
                        cell_[x][y][z].type(Cut);
                    N_[x][y][z] = I;

                    // Set the index pointer to the parent cell.
                    cell_[x][y][z].parent(Index_3(x, y, z));
                    if (cell_[x][y][z].type() == Fluid || cell_[x][y][z].type() == Cut) {
                        std::list<Nef_polyhedron::Point_3> points_3;
                        Nef_polyhedron::Vertex_const_iterator v;

                        // Calculate the centroid of the cell.
                        CGAL_forall_vertices(v, I)
                            points_3.push_back(v->point());
                        assert(points_3.size() >= 6);
                        cell_[x][y][z].centroid(CGAL::centroid(points_3.begin(), points_3.end(), CGAL::Dimension_tag<0>()));

                        // Calculate the volume of the cell, using a Triangulated
                        // volume and summing over tetrahedra.
                        Triangulation T(points_3.begin(), points_3.end());
                        assert(T.is_valid());
                        Kernel::FT volume = 0.0;
                        for (Triangulation::Finite_cells_iterator tcell = T.finite_cells_begin(); tcell != T.finite_cells_end(); ++tcell) {
                            assert(T.is_cell(tcell));
                            volume += T.tetrahedron(tcell).volume();
                        }
                        cell_[x][y][z].volume(volume);
                        assert(cell_[x][y][z].volume() > 0.0); // TODO should also be able to put an upper bound on the volume.

                        // Compute the normal vectors for each facet.
                        // FIXME Use Nef_polyhedron. plane() seems to be provided...?
                        std::transform(P.facets_begin(), P.facets_end(), P.planes_begin(),
                                       Normal_vector());

                        // FIXME convert to iterator over halffacets of Nef_polyhedron I.
                        for (Polyhedron::Facet_const_iterator fi = P.facets_begin(); fi != P.facets_end(); ++fi) {
                            Face newFace;
                            Polyhedron::Halfedge_const_handle h1, h2, h3;

                            assert(fi->is_triangle());
                            // circulate halfedges => vertices
                            h1 = fi->halfedge();
                            h2 = h1->next();
                            h3 = h2->next();
                            assert(h3->next() == h1);

                            // Compute the squared area of this triangular face.
                            newFace.area(CGAL::squared_area(h1->vertex()->point(), h2->vertex()->point(), h3->vertex()->point()));
                            assert(newFace.area() > 0.0); // TODO should also be able to put an upper bound on the area.

                            // Store the plane normal of this face.
                            newFace.normal(fi->plane());

                            // Compute the centroid of this triangular face.
                            points_3.clear();
                            points_3.push_back(h1->vertex()->point());
                            points_3.push_back(h2->vertex()->point());
                            points_3.push_back(h3->vertex()->point());
                            newFace.centroid(CGAL::centroid(points_3.begin(), points_3.end(), CGAL::Dimension_tag<0>()));

                            newFace.fluid(false); // FIXME

                            // Store this new face in the list of faces belonging
                            // to this cell.
                            cell_[x][y][z].addFace(Direction(fi->plane()), newFace);
                        }
                    }
                }
        }
    };

    V3Nef const& grid() const { return N_; }
    V3Cell const& cell() const { return cell_; }

    private:
    // TODO: This function should be const, but
    // Nef_polyhedron::convert_to_polyhedron is not const, so it can't be.
    std::ostream& output_vrml(std::ostream& out) {
        CGAL::VRML_2_ostream vrml_out(out);
        for (int x = 0; x < N_.size(); ++x)
            for (int y = 0; y < N_[x].size(); ++y)
                for (int z = 0; z < N_[x][y].size(); ++z) {
                    Polyhedron P;
                    // Convert this new cut Nef_polyhedron into the Polyhedron P.
                    N_[x][y][z].convert_to_polyhedron(P);
                    // Output the Polyhedron in VRML format.
                    vrml_out << P;
                }
        return out;
    }
    std::ostream& output_nef(std::ostream& out) const {
        Nef_polyhedron big;
        for (int x = 0; x < N_.size(); ++x)
            for (int y = 0; y < N_[x].size(); ++y)
                for (int z = 0; z < N_[x][y].size(); ++z) {
                    // Output the Nef_polyhedron in NEF format.
                    big += N_[x][y][z];
                }
        out << big;
        return out;
    }
    std::ostream& output_cgns(std::ostream& out) const {
        int index_file = 0, index_base = 0, index_zone = 0;
        char const * const NAME = "/tmp/temp.cgns";

        if (cg_open(NAME, CG_MODE_WRITE, &index_file) != CG_OK) {
            std::cerr << cg_get_error() << std::endl;
            return out;
        }
        if (cg_base_write(index_file, "Base", 3, 3, &index_base) != CG_OK) {
            std::cerr << cg_get_error() << std::endl;
            (void)cg_close(index_file);
            std::remove(NAME);
            return out;
        }
        int isize[3][3] = {0}; // TODO Why is isize 3*3?
        int grid_size[3] = { N_.size(), N_[0].size(), N_[0][0].size() };
        // vertex size
        isize[0][0] = grid_size[0] * grid_size[1] * grid_size[2];
        // cell size
        isize[1][0] = (grid_size[0] - 1) * (grid_size[1] - 1) * (grid_size[2] - 1);
        // boundary size
        isize[2][0] = 0;
        if (cg_zone_write(index_file, index_base, "Zone 1", reinterpret_cast<int*>(isize), Unstructured, &index_zone) != CG_OK) {
            std::cerr << cg_get_error() << std::endl;
            (void)cg_close(index_file);
            std::remove(NAME);
        }
        if (cg_close(index_file) != CG_OK) {
            std::cerr << cg_get_error() << std::endl;
            std::remove(NAME);
        }
        std::ifstream in(NAME);
        out << in.rdbuf();
        in.close();
        std::remove(NAME);
        return out;
    }
    static int alloc() { return alloc_; }
    V3Nef N_;
    V3Cell cell_;
    static const int alloc_;
    DISALLOW_COPY_AND_ASSIGN(Grid);
};
const int Grid::alloc_ = std::ios_base::xalloc();

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

// TODO: This function should take a Grid const&, but output_vrml can't be
// const yet (see above), so it can't.
std::ostream& operator<<(std::ostream& out, Grid & g) {
    // Call the appropriate output function on the Grid object depending on
    // the stream state.
    switch (out.iword(g.alloc())) {
        default:
        case GridFormat::OUTPUT_VRML:
            return g.output_vrml(out);
            break;
        case GridFormat::OUTPUT_NEF:
            return g.output_nef(out);
            break;
        case GridFormat::OUTPUT_CGNS:
            return g.output_cgns(out);
    }
}

int main() {
    const int NX = 5, NY = 5, NZ = 5;
    Nef_polyhedron N1;

    std::istringstream in(cube); // Load the definition of a unit cube from string
    in >> N1;                    // into the test cube
    assert(N1.number_of_vertices() == 8);
    assert(N1.number_of_facets() == 6);
    assert(N1.number_of_edges() == 12);
    assert(N1.number_of_volumes() == 2);

    // Put the test cube at an appropriate place
    Aff_transformation Aff1(CGAL::TRANSLATION, Vector(0.5, 0.5, 0.5));
    Aff_transformation Aff2(CGAL::SCALING, 1.5);
    N1.transform(Aff1);
    N1.transform(Aff2);

    // Create the cutting object.
    Grid grid(NX, NY, NZ);

    // Cut the grid to the test cube.
    grid.cut(N1);

    // Output the grid in NEF format.
    std::cout << GridFormat(GridFormat::OUTPUT_CGNS) << grid << std::endl;

    return 0;
}
