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

#include <cutcell.hpp>
#include <cgnslib.h>
#include <CGAL/IO/Polyhedron_VRML_2_ostream.h>
#include <CGAL/centroid.h>
#include <CGAL/Triangulation_3.h>
#include <list>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>

namespace cutcell {

typedef CGAL::Triangulation_3<cutcell::Kernel> Triangulation;

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

Grid::Grid(int X, int Y, int Z) {
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
                Aff_transformation Aff(TRANSLATION, Vector(x, y, z));
                N_[x][y][z].transform(Aff);
            }
};
void Grid::cut(Nef_polyhedron const& N1) {
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

std::ostream& Grid::output_vrml(std::ostream& out) const {
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
std::ostream& Grid::output_nef(std::ostream& out) const {
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
std::ostream& Grid::output_cgns(std::ostream& out) const {
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

const int Grid::alloc_ = std::ios_base::xalloc();

std::ostream& operator<<(std::ostream& out, Grid const & g) {
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

} // namespace cutcell
