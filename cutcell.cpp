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
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/convex_decomposition_3.h>

#include <list>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <string>

namespace cutcell {

typedef CGAL::Triangulation_3<cutcell::Kernel> Triangulation;
typedef CGAL::Delaunay_triangulation_3<cutcell::Kernel> Delaunay_triangulation;

void Grid::cut(Nef_polyhedron const& N1) {
    // Ensure that each multi_array is 3 dimensional
    assert(N_.num_dimensions() == 3 && cell_.num_dimensions() == 3);
    // And that both arrays are the same shape.
    assert(std::equal(N_.shape(), N_.shape() + N_.num_dimensions(), cell_.shape()));
    // Create the Unit Cube from the global string definition.
    Nef_polyhedron UnitCube;
    std::istringstream in(cube);
    in >> UnitCube;
    assert(UnitCube.number_of_vertices() == 8);
    assert(UnitCube.number_of_facets() == 6);
    assert(UnitCube.number_of_edges() == 12);
    assert(UnitCube.number_of_volumes() == 2);

    for (V3NefIndex x = 0; x < N_.shape()[0]; ++x) {
        for (V3NefIndex y = 0; y < N_.shape()[1]; ++y)
            for (V3NefIndex z = 0; z < N_.shape()[2]; ++z) {
                #ifndef NDEBUG
                std::cerr << "Grid cell at " << x << ", " << y << ", " << z << " is a ";
                #endif
                // Compute the intersection of this part of the grid with the
                // test cube.
                Nef_polyhedron I = UnitCube;
                Aff_transformation Aff(TRANSLATION, Vector(static_cast<int>(x), static_cast<int>(y), static_cast<int>(z)));
                I.transform(Aff);
                Nef_polyhedron Nnew = I - N1;

                // Set the type of the new cell.
                if (Nnew.is_empty()) {
                    // No points, must be completely inside the solid.
                    cell_[x][y][z].type(Solid);
                    #ifndef NDEBUG
                    std::cerr << "Solid";
                    #endif
                }
                else if (Nnew == I) {
                    // Unchanged, must be completely outside the solid.
                    cell_[x][y][z].type(Fluid);
                    #ifndef NDEBUG
                    std::cerr << "Fluid";
                    #endif
                }
                else {
                    // Something else, must be a cut cell.
                    cell_[x][y][z].type(Cut);
                    // Assign it to the grid.
                    N_[x][y][z] = Nnew;
                    #ifndef NDEBUG
                    std::cerr << "Cut";
                    #endif
                }

                #ifndef NDEBUG
                std::cerr << " cell." << std::endl;
                #endif
                // Set the index pointer to the parent cell.
                cell_[x][y][z].parent(Index_3(x, y, z));
                if (cell_[x][y][z].type() == Fluid || cell_[x][y][z].type() == Cut) {
                    std::list<Nef_polyhedron::Point_3> points_3;
                    Nef_polyhedron::Vertex_const_iterator v;

                    // Calculate the centroid of the cell.
                    CGAL_forall_vertices(v, Nnew)
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

                    continue; // FIXME This stuff only works if the faces are triangular.
                    // FIXME convert to iterator over halffacets of Nef_polyhedron I.
                    // Nef_polyhedrons have Halffacets
                    // Halffacets have halfedge_cycles
                    // halfedge_cycles are SHalfedges or one SHalfloop
                    Nef_polyhedron::Halffacet_const_iterator fi;
                    CGAL_forall_facets(fi, Nnew) {
                        Face newFace;
                        Nef_polyhedron::SHalfedge_const_handle h[3];
                        size_t i = 0;

                        // circulate 3 halfedges => vertices
                        //for (Nef_polyhedron::Halffacet_cycle_const_iterator hfi = fi->facet_cycles_begin(); hfi != fi->facet_cycles_end(); ++hfi) {
                        Nef_polyhedron::Halffacet_cycle_const_iterator hfc;
                        CGAL_forall_facet_cycles_of(hfc, fi) {
                            assert(hfc.is_shalfedge());
                            Nef_polyhedron::SHalfedge_const_handle se(hfc);
                            Nef_polyhedron::SHalfedge_around_facet_const_circulator hc(se);
                            Nef_polyhedron::SHalfedge_around_facet_const_circulator hc_end(hc);
                            CGAL_For_all(hc, hc_end) {
                                assert(i < sizeof(h)/sizeof(h[0]));
                                h[i] = Nef_polyhedron::SHalfedge_const_handle(hc);
                                std::cerr << "SHalfedge from " << h[i]->source()->point() << " to " << h[i]->target()->point() << std::endl;
                                ++i;
                            }
                        }
                        std::cerr << i << std::endl;
                        if (i != 3)
                            continue;
                        assert(i == 3); // This should have been a triangular facet.

                        // Compute the squared area of this triangular face.
                        newFace.area(CGAL::squared_area(h[0]->source()->point(), h[1]->source()->point(), h[2]->source()->point()));
                        assert(newFace.area() > 0.0); // TODO should also be able to put an upper bound on the area.

                        // Store the plane normal of this face.
                        newFace.normal(fi->plane());

                        // Compute the centroid of this triangular face.
                        points_3.clear();
                        points_3.push_back(h[0]->source()->point());
                        points_3.push_back(h[1]->source()->point());
                        points_3.push_back(h[2]->source()->point());
                        newFace.centroid(CGAL::centroid(points_3.begin(), points_3.end(), CGAL::Dimension_tag<0>()));

                        newFace.fluid(false); // FIXME

                        // Store this new face in the list of faces belonging
                        // to this cell.
                        cell_[x][y][z].addFace(newFace);
                    }
                }
            }
    }
};

std::ostream& Grid::output_vrml(std::ostream& out) const {
    CGAL::VRML_2_ostream vrml_out(out);
    for (V3NefIndex x = 0; x < N_.shape()[0]; ++x)
        for (V3NefIndex y = 0; y < N_.shape()[1]; ++y)
            for (V3NefIndex z = 0; z < N_.shape()[2]; ++z) {
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
    for (V3NefIndex x = 0; x < N_.shape()[0]; ++x)
        for (V3NefIndex y = 0; y < N_.shape()[1]; ++y)
            for (V3NefIndex z = 0; z < N_.shape()[2]; ++z) {
                // Output the Nef_polyhedron in NEF format.
                big += N_[x][y][z];
            }
    out << big;
    return out;
}

template <class Element>
static typename std::vector<Element>::difference_type find_or_insert_index_of(std::vector<Element>& points_3, Element const & point, std::vector<double>& xvec, std::vector<double>& yvec, std::vector<double>& zvec) {
    typename std::vector<Element>::iterator it;
    for (it = points_3.begin(); it != points_3.end(); ++it)
        if (*it == point)
            break;
    if (it == points_3.end()) { // if point not in points_3
        xvec.push_back(CGAL::to_double(point.x()));
        yvec.push_back(CGAL::to_double(point.y()));
        zvec.push_back(CGAL::to_double(point.z()));
        points_3.push_back(point);
        #ifndef NDEBUG
        std::cerr << "Added point " << points_3[points_3.size() - 1] << ": " << xvec[xvec.size() - 1] << " " << yvec[yvec.size() - 1] << " " << zvec[zvec.size() - 1] << std::endl;
        #endif
        // it may have been invalidated by push_back. Make sure we have an iterator at the new point.
        it = points_3.end();
        --it;
    }
    // Push the index of the element in the points_3 set into the connectivity list.
    typename std::vector<Element>::difference_type element_index = std::distance(points_3.begin(), it);
    #ifndef NDEBUG
    std::cerr << "element index: " << element_index << std::endl;
    #endif
    return element_index;
}

int Grid::output_cgns_file(std::string const& name) const {
    int hexa_8_cells = 0, tetra_4_cells = 0;
    // Loop over all the Vertices (Points) in the grid and _uniquely_ add them
    // to the xvec, yvec, zvec vectors of the respective 3D coordinates.
    // The points_3 vector is the key for a naive implementation of a map.
    std::vector<double> xvec, yvec, zvec;
    std::vector<Nef_polyhedron::Point_3> points_3;
    // CGNS requires the indices into the xvec, yvec, zvec vectors to identify
    // the nodes of the polyhedra (in this case only hexahedrons and
    // tetrahedrons are used).
    std::vector<int> hexa_8_elements, tetra_4_elements;
    for (V3NefIndex x = 0; x < N_.shape()[0]; ++x)
        for (V3NefIndex y = 0; y < N_.shape()[1]; ++y)
            for (V3NefIndex z = 0; z < N_.shape()[2]; ++z) {
                if (cell_[x][y][z].type() == Fluid) {
                    // This is a cube.
                    #ifndef NDEBUG
                    std::cerr << "Processing fluid cell at " << x << " " << y << " " << z << std::endl;
                    int nodes = 0;
                    #endif
                    for (unsigned dx = 0; dx < 2; ++dx)
                        for (unsigned dy = 0; dy < 2; ++dy)
                            for (unsigned dz = 0; dz < 2; ++dz) {
                                Nef_polyhedron::Point_3 point(static_cast<double>(x + dx), static_cast<double>(y + dy), static_cast<double>(z + (dy == 0 ? (1 - dz) : dz)));
                                #ifndef NDEBUG
                                std::cerr << "Fluid cell point at " << point << std::endl;
                                ++nodes;
                                #endif
                                hexa_8_elements.push_back(find_or_insert_index_of(points_3, point, xvec, yvec, zvec) + 1); // indices are 1-based.
                            }
                    #ifndef NDEBUG
                    assert(nodes == 8);
                    #endif
                    ++hexa_8_cells;
                }
                else if (cell_[x][y][z].type() == Cut) {
                    #ifndef NDEBUG
                    std::cerr << "Processing cut cell at " << x << " " << y << " " << z << std::endl;
                    #endif
                    Nef_polyhedron N(N_[x][y][z]);
                    CGAL::convex_decomposition_3(N);
                    // the first volume is the outer volume, which is ignored in the decomposition
                    Nef_polyhedron::Volume_const_iterator vit = ++N.volumes_begin();
                    for ( ; vit != N.volumes_end(); ++vit) {
                        if (!vit->mark()) // FIXME: I don't know what this does...
                            continue;
                        Polyhedron P;
                        N.convert_inner_shell_to_polyhedron(vit->shells_begin(), P);
                        Delaunay_triangulation T(P.points_begin(), P.points_end());
                        for (Delaunay_triangulation::Finite_cells_iterator it = T.finite_cells_begin(); it != T.finite_cells_end(); ++it) {
                            for (int i = 0; i < 4; ++i) {
                                Nef_polyhedron::Point_3 point = it->vertex(i)->point();
                                tetra_4_elements.push_back(find_or_insert_index_of(points_3, point, xvec, yvec, zvec) + 1); // indices are 1-based.
                            }
                            ++tetra_4_cells;
                        }
                    }
                }
            }
    int index_file = 0, index_base = 0, index_zone = 0;
    int index_coord_x = 0, index_coord_y = 0, index_coord_z = 0;
    int index_section = 0;
    if (cg_open(name.c_str(), CG_MODE_WRITE, &index_file) != CG_OK) {
        return 1;
    }
    if (cg_base_write(index_file, "Base", 3, 3, &index_base) != CG_OK) {
        std::cerr << cg_get_error() << std::endl;
        (void)cg_close(index_file);
        return 1;
    }
    int isize[9] = {0};
    // vertex size
    isize[0] = points_3.size();
    // cell size
    assert(hexa_8_cells + tetra_4_cells > 0);
    isize[1] = hexa_8_cells + tetra_4_cells;
    // boundary size
    //isize[0][2] = 0;
    if (cg_zone_write(index_file, index_base, "Zone1", isize, Unstructured, &index_zone) != CG_OK) {
        std::cerr << cg_get_error() << std::endl;
        (void)cg_close(index_file);
        return 1;
    }
    #ifndef NDEBUG
    std::cerr << "Writing " << xvec.size() << " coordinates..." << std::endl;
    #endif
    // We told cg_zone_write to expect isize[0] vertex coordinates in each dimension.
    assert(xvec.size() >= isize[0]);
    assert(yvec.size() >= isize[0]);
    assert(zvec.size() >= isize[0]);
    if (cg_coord_write(index_file, index_base, index_zone, RealDouble, "CoordinateX", &xvec[0], &index_coord_x) != CG_OK ||
        cg_coord_write(index_file, index_base, index_zone, RealDouble, "CoordinateY", &yvec[0], &index_coord_y) != CG_OK ||
        cg_coord_write(index_file, index_base, index_zone, RealDouble, "CoordinateZ", &zvec[0], &index_coord_z) != CG_OK) {
        std::cerr << cg_get_error() << std::endl;
        (void)cg_close(index_file);
        return 1;
    }
    if (hexa_8_cells > 0) {
        #ifndef NDEBUG
        std::cerr << "Writing " << hexa_8_cells << " Fluid elements..." << std::endl;
        #endif
        assert(hexa_8_elements.size() >= hexa_8_cells);
        assert(hexa_8_elements.size() / 8 == hexa_8_cells);
        if (cg_section_write(index_file, index_base, index_zone, "FluidGridElements", HEXA_8, 1, hexa_8_cells, /* nbndry = */ 0, &hexa_8_elements[0], &index_section) != CG_OK) {
            std::cerr << cg_get_error() << std::endl;
            (void)cg_close(index_file);
            return 1;
        }
    }
    if (tetra_4_cells > 0) {
        #ifndef NDEBUG
        std::cerr << "Writing " << tetra_4_cells << " Cut cell elements..." << std::endl;
        #endif
        assert(tetra_4_elements.size() >= tetra_4_cells);
        assert(tetra_4_elements.size() / 4 == tetra_4_cells);
        if (cg_section_write(index_file, index_base, index_zone, "CutGridElements", TETRA_4, 1, tetra_4_cells, /* nbndry = */ 0, &tetra_4_elements[0], &index_section) != CG_OK) {
            std::cerr << cg_get_error() << std::endl;
            (void)cg_close(index_file);
            return 1;
        }
    }

    if (cg_close(index_file) != CG_OK) {
        std::cerr << cg_get_error() << std::endl;
        return 1;
    }
    return 0;
}
std::ostream& Grid::output_cgns(std::ostream& out) const {
    const std::string NAME("/tmp/temp.cgns");

    if (output_cgns_file(NAME) != 0) {
        std::remove(NAME.c_str());
        return out;
    }
    std::ifstream in(NAME.c_str());
    out << in.rdbuf();
    in.close();
    std::remove(NAME.c_str());
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
