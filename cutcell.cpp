/* Cartesian cut cell grid generator
 *
 * Generate a 3D Cartesian array of unit cubes, with the given solid
 * subtracted from them. Export various properties of the cubes.
 *
 * Copyright 2010,2011 Bruce Duncan, University of Edinburgh
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
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/convex_decomposition_3.h>
#include <CGAL/bounding_box.h>

#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <string>

namespace cutcell {

typedef CGAL::Delaunay_triangulation_3<cutcell::Kernel> Delaunay_triangulation;

Grid::Grid(int X, int Y, int Z) : N_(boost::extents[X][Y][Z]), cell_(boost::extents[X][Y][Z]), nCells_(), N1_() {
    // Ensure that each multi_array is 3 dimensional
    assert(N_.num_dimensions() == 3 && cell_.num_dimensions() == 3);
    // And that both arrays are the same shape.
    assert(std::equal(N_.shape(), N_.shape() + N_.num_dimensions(), cell_.shape()));
    // Create the Unit Cube from the global string definition.
    std::istringstream in(cube);
    in >> UnitCube_;
    assert(UnitCube_.number_of_vertices() == 8);
    assert(UnitCube_.number_of_facets() == 6);
    assert(UnitCube_.number_of_edges() == 12);
    assert(UnitCube_.number_of_volumes() == 2);
}

void Grid::addSolid(Nef_polyhedron const& N) {
    #ifndef NDEBUG
    std::cerr << "Adding solid with " << N.number_of_vertices() << " vertices." << std::endl;
    std::cerr << "Internal solid has " << N1_.number_of_vertices() << " vertices." << std::endl;
    #endif
    if (N1_.is_empty()) {
        N1_ = N;
    }
    else {
        N1_ += N;
    }
    #ifndef NDEBUG
    std::cerr << "After union, internal solid has " << N1_.number_of_vertices() << " vertices." << std::endl;
    #endif
}

void Grid::cut() {
    // Compute the iso-oriented bounding box of the solid.
    std::vector<Nef_polyhedron::Point_3> points_3;
    points_3.reserve(N1_.number_of_vertices());
    assert(points_3.size() == 0);
    Nef_polyhedron::Vertex_const_iterator v;
    CGAL_forall_vertices(v, N1_)
        points_3.push_back(v->point());
    assert(points_3.size() == N1_.number_of_vertices());
    Kernel::Iso_cuboid_3 c3 = CGAL::bounding_box(points_3.begin(), points_3.end());
    std::vector<int> bb;
    bb.push_back(floor(CGAL::to_double(c3.min_coord(0))));
    bb.push_back(floor(CGAL::to_double(c3.min_coord(1))));
    bb.push_back(floor(CGAL::to_double(c3.min_coord(2))));
    bb.push_back(ceil(CGAL::to_double(c3.max_coord(0))));
    bb.push_back(ceil(CGAL::to_double(c3.max_coord(1))));
    bb.push_back(ceil(CGAL::to_double(c3.max_coord(2))));
    #ifndef NDEBUG
    std::cerr << "Begin cutting. Bounding box of solid object: " << bb[0] << " " << bb[1] << " " << bb[2] << ", "
                                                                 << bb[3] << " " << bb[4] << " " << bb[5] << std::endl;
    #endif
    for (V3NefIndex x = 0; x < N_.shape()[0]; ++x) {
        for (V3NefIndex y = 0; y < N_.shape()[1]; ++y)
            for (V3NefIndex z = 0; z < N_.shape()[2]; ++z) {
                if (x < bb[0] || y < bb[1] || z < bb[2] || x >= bb[3] || y >= bb[4] || z >= bb[5]) {
                    #ifndef NDEBUG
                    std::cerr << "Grid cell at " << x << ", " << y << ", " << z << " is outside the bounding box and therefore Fluid." << std::endl;
                    #endif
                    // This is kind of awkward. This code is repeated below in the Fluid path.
                    cell_[x][y][z].type(Fluid);
                    ++nCells_[Fluid];
                    continue;
                }
                // Compute the intersection of this part of the grid with the
                // test cube.
                Nef_polyhedron I(UnitCube_);
                Aff_transformation Aff(TRANSLATION, Vector(static_cast<int>(x), static_cast<int>(y), static_cast<int>(z)));
                I.transform(Aff);
                Nef_polyhedron Nnew(I - N1_);

                // Set the type of the new cell.
                if (Nnew.is_empty()) {
                    // No points, must be completely inside the solid.
                    cell_[x][y][z].type(Solid);
                    ++nCells_[Solid];
                }
                else if (Nnew == I) {
                    // Unchanged, must be completely outside the solid.
                    cell_[x][y][z].type(Fluid);
                    ++nCells_[Fluid];
                }
                else {
                    // Something else, must be a cut cell.
                    cell_[x][y][z].type(Cut);
                    // Ensure the polyhedron is convex.
                    CGAL::convex_decomposition_3(Nnew);
                    // Assign it to the grid.
                    N_[x][y][z] = Nnew;
                    ++nCells_[Cut];
                }

                #ifndef NDEBUG
                std::cerr << "Grid cell at " << x << ", " << y << ", " << z << " is a " << Typenames[cell_[x][y][z].type()] << " cell." << std::endl;
                #endif
            }
    }
    #ifndef NDEBUG
    std::cerr << "Finished cutting. Grid has:" << std::endl;
    for (unsigned i = 0; i < NUM_TYPES; ++i) {
        std::cerr << nCells_[i] << " " << Typenames[i] << " cells." << std::endl;
    }
    std::cerr << std::endl;
    #endif
}

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
    // Iterating backwards might be faster, due to locality.
    typename std::vector<Element>::reverse_iterator it;
    for (it = points_3.rbegin(); it != points_3.rend(); ++it)
        if (*it == point)
            break;
    if (it == points_3.rend()) { // if point not in points_3
        xvec.push_back(CGAL::to_double(point.x()));
        yvec.push_back(CGAL::to_double(point.y()));
        zvec.push_back(CGAL::to_double(point.z()));
        points_3.push_back(point);
        #ifndef NDEBUG
        std::cerr << "Added point " << points_3[points_3.size() - 1] << ": " << xvec[xvec.size() - 1] << " " << yvec[yvec.size() - 1] << " " << zvec[zvec.size() - 1] << std::endl;
        #endif
        // it may have been invalidated by push_back. Make sure we have an iterator at the new point.
        it = points_3.rbegin();
        assert(*it == point);
    }
    // Push the index of the element in the points_3 set into the connectivity list.
    // Because std::vector has bidirectional random-access iterators, std::distance from rend will be negative and 1-based. Make sure it's positive.
    typename std::vector<Element>::difference_type element_index = -std::distance(--points_3.rend(), it);
    #ifndef NDEBUG
    // This is the alternative, potentially slower, algorithm for finding the location of point in points_3
    // Test that it gets the same answer as the one above.
    // For this we need a forward iterator...
    typename std::vector<Element>::iterator itf;
    for (itf = points_3.begin(); itf != points_3.end(); ++itf)
        if (*itf == point)
            break;
    if (itf == points_3.end()) {
        std::cerr << "Could not find point " << point << " in points_3 vector! Aborting..." << std::endl;
        abort();
    }
    std::cerr << "element index: " << element_index << std::endl;
    assert(element_index == std::distance(points_3.begin(), itf));
    #endif
    return element_index;
}

int Grid::output_cgns_file(std::string const& name) const {
    int hexa_8_cells = 0, tetra_4_cells = 0, quad_4_faces = 0, tri_3_faces = 0;
    // Loop over all the Vertices (Points) in the grid and _uniquely_ add them
    // to the xvec, yvec, zvec vectors of the respective 3D coordinates.
    // The points_3 vector is the key for a naive implementation of a map.

    std::vector<double> xvec, yvec, zvec;
    std::vector<Nef_polyhedron::Point_3> points_3;
    // Try to save time later by allocating the space we need now. For this we need the volume and area of the domain:
    unsigned const cell_area = 2 * ((N_.shape()[0] * N_.shape()[1]) + (N_.shape()[1] * N_.shape()[2]) + (N_.shape()[2] * N_.shape()[0]));
    unsigned const cell_volume = N_.shape()[0] * N_.shape()[1] * N_.shape()[2];
    unsigned const points_area = 2 * (((N_.shape()[0] + 1) * (N_.shape()[1] + 1)) + ((N_.shape()[1] + 1) * (N_.shape()[2] + 1)) + ((N_.shape()[2] + 1) * (N_.shape()[0] + 1)));
    unsigned const points_volume = (N_.shape()[0] + 1) * (N_.shape()[1] + 1) * (N_.shape()[2] + 1);
    // There are approximately X*Y*Z Fluid points plus twice the number of solid points (one on the solid, one on the cut cell boundary)
    unsigned points_capacity = points_volume + 2 * N1_.number_of_vertices();
    #ifndef NDEBUG
    std::cerr << "Grid extent (size + 1) is " << (N_.shape()[0] + 1) << " x " << (N_.shape()[1] + 1) << " x " << (N_.shape()[2] + 1) << " = " << points_volume << " vertices." << std::endl;
    std::cerr << "Solid has " << N1_.number_of_vertices() << " points (expecting twice this number of output vertices)." << std::endl;
    std::cerr << "Pre-allocating " << (3 * sizeof(double) * points_capacity) << " bytes for double points." << std::endl;
    #endif
    xvec.reserve(points_capacity);
    yvec.reserve(points_capacity);
    zvec.reserve(points_capacity);
    #ifndef NDEBUG
    std::cerr << "Pre-allocating " << (sizeof(Nef_polyhedron::Point_3) * points_capacity) << " bytes for Nef_polyhedron::Point_3s." << std::endl;
    #endif
    points_3.reserve(points_capacity);

    // CGNS requires the indices into the xvec, yvec, zvec vectors to identify
    // the nodes of the polyhedra (in this case only hexahedrons and
    // tetrahedrons are used).
    std::vector<int> hexa_8_elements, tetra_4_elements;
    unsigned hexa_8_capacity = 8 * cell_volume, tetra_4_capacity = 4 * N1_.number_of_vertices();
    #ifndef NDEBUG
    std::cerr << "Pre-allocating " << (sizeof(int) * hexa_8_capacity) << " bytes for HEXA_8 elements." << std::endl;
    std::cerr << "Pre-allocating " << (sizeof(int) * tetra_4_capacity) << " bytes for TETRA_4 elements." << std::endl;
    #endif
    hexa_8_elements.reserve(hexa_8_capacity);
    tetra_4_elements.reserve(tetra_4_capacity);

    // For the boundary, we write QUAD_4 and TRI_3 elements in a MIXED section.
    std::vector<int> boundary_elements;
    unsigned boundary_capacity = 5 * cell_area;
    #ifndef NDEBUG
    std::cerr << "Pre-allocating " << (sizeof(int) * boundary_capacity) << " bytes for QUAD_4 MIXED elements." << std::endl;
    #endif
    boundary_elements.reserve(boundary_capacity);

    // We should store the Delaunay_triangulation of the Cut cells to ensure that we get the same answer for both the volume and the boundary.
    typedef boost::multi_array<std::vector<Delaunay_triangulation>, 3> V3Tri;
    V3Tri T(boost::extents[N_.shape()[0]][N_.shape()[1]][N_.shape()[2]]);

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
                    Nef_polyhedron & N = N_[x][y][z];
                    // the first volume is the outer volume, which is ignored in the decomposition
                    Nef_polyhedron::Volume_const_iterator vit = ++N.volumes_begin();
                    for ( ; vit != N.volumes_end(); ++vit) {
                        if (!vit->mark()) // FIXME: I don't know what this does...
                            continue;
                        Polyhedron P;
                        N.convert_inner_shell_to_polyhedron(vit->shells_begin(), P);
                        /* For some reason T[x][y][z].insert(...) produces the wrong result!
                           Using the default constructor and then calling .insert() on a temporary also gives the correct result.
                        */
                        Delaunay_triangulation Tnew(P.points_begin(), P.points_end());
                        for (Delaunay_triangulation::Finite_cells_iterator it = Tnew.finite_cells_begin(); it != Tnew.finite_cells_end(); ++it) {
                            for (int i = 0; i < 4; ++i) {
                                Nef_polyhedron::Point_3 point = it->vertex(i)->point();
                                tetra_4_elements.push_back(find_or_insert_index_of(points_3, point, xvec, yvec, zvec) + 1); // indices are 1-based.
                            }
                            ++tetra_4_cells;
                        }
                        T[x][y][z].push_back(Tnew);
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
        std::cerr << "Writing " << hexa_8_cells << " Fluid elements with " << hexa_8_elements.size() << " nodes..." << std::endl;
        #endif
        assert(hexa_8_elements.size() >= hexa_8_cells);
        assert(hexa_8_elements.size() / 8 == hexa_8_cells);
        if (cg_section_write(index_file, index_base, index_zone, "1_FluidGridElements", HEXA_8, 1, hexa_8_cells, /* nbndry = */ 0, &hexa_8_elements[0], &index_section) != CG_OK) {
            std::cerr << cg_get_error() << std::endl;
            (void)cg_close(index_file);
            return 1;
        }
    }
    if (tetra_4_cells > 0) {
        #ifndef NDEBUG
        std::cerr << "Writing " << tetra_4_cells << " Cut cell elements with " << tetra_4_elements.size() << " nodes..." << std::endl;
        #endif
        assert(tetra_4_elements.size() >= tetra_4_cells);
        assert(tetra_4_elements.size() / 4 == tetra_4_cells);
        if (cg_section_write(index_file, index_base, index_zone, "2_CutGridElements", TETRA_4, hexa_8_cells + 1, hexa_8_cells + tetra_4_cells, /* nbndry = */ 0, &tetra_4_elements[0], &index_section) != CG_OK) {
            std::cerr << cg_get_error() << std::endl;
            (void)cg_close(index_file);
            return 1;
        }
    }

    int index_bc = 0;
    // CGNS wants the element numbers of the faces which make up the boundary. Count the ones we are about to add, starting from the end of the 3D cells.
    int element_count = hexa_8_cells + tetra_4_cells;
    std::vector<int> quad_4_face_list;
    std::vector<int> tri_3_face_list;
    // Boundary cells
    // Left/right boundary
    for (V3NefIndex x = 0; x <= N_.shape()[0]; x += N_.shape()[0]) { // x = [0, NX]
        for (V3NefIndex y = 0; y < N_.shape()[1]; ++y) {
            for (V3NefIndex z = 0; z < N_.shape()[2]; ++z) {
                const int cx = x > 0 ? x - 1 : x, cy = y, cz = z; // It's helpful to refer to the cell_ and T arrays, which have an index one less than the upper boundary.
                if (cell_[cx][cy][cz].type() == Fluid) {
                    boundary_elements.push_back(QUAD_4);
                    for (unsigned dy = 0; dy < 2; ++dy) {
                        for (unsigned dz = 0; dz < 2; ++dz) {
                            Nef_polyhedron::Point_3 point(static_cast<double>(x), static_cast<double>(y + dy), static_cast<double>(z + (dy == 0 ? (1 - dz) : dz)));
                            #ifndef NDEBUG
                            std::cerr << "Fluid cell boundary point at " << point << std::endl;
                            #endif
                            boundary_elements.push_back(find_or_insert_index_of(points_3, point, xvec, yvec, zvec) + 1);
                        }
                    }
                    ++quad_4_faces;
                    quad_4_face_list.push_back(++element_count);
                }
                else if (cell_[cx][cy][cz].type() == Cut) {
                    // For a cut cell, we want to write TRI_3 elements of the points which are on
                    // the plane which represents the boundary. Pick the three points out of any of
                    // the tetrahedra which meet this condition.
                    Nef_polyhedron::Plane_3 plane(Nef_polyhedron::Point_3(static_cast<double>(x), static_cast<double>(y), static_cast<double>(z)),
                                                  Nef_polyhedron::Vector_3(x == 0 ? -1 : 1, 0, 0));
                    for (std::vector<Delaunay_triangulation>::iterator tri_it = T[cx][cy][cz].begin(); tri_it != T[cx][cy][cz].end(); ++tri_it) {
                        for (Delaunay_triangulation::Finite_cells_iterator it = tri_it->finite_cells_begin(); it != tri_it->finite_cells_end(); ++it) {
                            std::vector<Nef_polyhedron::Point_3> points;
                            for (int i = 0; i < 4; ++i) {
                                Nef_polyhedron::Point_3 point = it->vertex(i)->point();
                                if (plane.has_on(point)) {
                                    points.push_back(point);
                                }
                            }
                            if (points.size() == 3) {
                                boundary_elements.push_back(TRI_3);
                                for (std::vector<Nef_polyhedron::Point_3>::iterator pointit = points.begin(); pointit != points.end(); ++pointit) {
                                    boundary_elements.push_back(find_or_insert_index_of(points_3, *pointit, xvec, yvec, zvec) + 1); // indices are 1-based.
                                }
                                ++tri_3_faces;
                                tri_3_face_list.push_back(++element_count);
                            }
                        }
                    }
                }
            }
        }
        if (quad_4_faces > 0) {
            #ifndef NDEBUG
            std::cerr << "Writing " << (x == 0 ? "Ilo" : "Ihi") << " fluid boundary as list of size " << quad_4_face_list.size() << "." << std::endl;
            #endif
            if (cg_boco_write(index_file, index_base, index_zone, x == 0 ? "IloFluid" : "IhiFluid", BCWall, ElementList, quad_4_face_list.size(), &quad_4_face_list[0], &index_bc) != CG_OK) {
                std::cerr << cg_get_error() << std::endl;
                (void)cg_close(index_file);
                return 1;
            }
            quad_4_face_list.clear();
        }

        if (tri_3_faces > 0) {
            #ifndef NDEBUG
            std::cerr << "Writing " << (x == 0 ? "Ilo" : "Ihi") << " cut boundary as list of size " << tri_3_face_list.size() << "." << std::endl;
            #endif
            if (cg_boco_write(index_file, index_base, index_zone, x == 0 ? "IloCut" : "IhiCut", BCWall, ElementList, tri_3_face_list.size(), &tri_3_face_list[0], &index_bc) != CG_OK) {
                std::cerr << cg_get_error() << std::endl;
                (void)cg_close(index_file);
                return 1;
            }
            tri_3_face_list.clear();
        }

        tri_3_faces = 0;
        quad_4_faces = 0;
    }
    // Top/Bottom boundary
    for (V3NefIndex y = 0; y <= N_.shape()[1]; y += N_.shape()[1]) { // y = [0, NY]
        for (V3NefIndex x = 0; x < N_.shape()[0]; ++x) {
            for (V3NefIndex z = 0; z < N_.shape()[2]; ++z) {
                const int cx = x, cy = y > 0 ? y - 1 : y, cz = z; // It's helpful to refer to the cell_ and T arrays, which have an index one less than the upper boundary.
                if (cell_[cx][cy][cz].type() == Fluid) {
                    boundary_elements.push_back(QUAD_4);
                    for (unsigned dx = 0; dx < 2; ++dx) {
                        for (unsigned dz = 0; dz < 2; ++dz) {
                            Nef_polyhedron::Point_3 point(static_cast<double>(x + dx), static_cast<double>(y), static_cast<double>(z + (dx == 0 ? (1 - dz) : dz)));
                            #ifndef NDEBUG
                            std::cerr << "Fluid cell boundary point at " << point << std::endl;
                            #endif
                            boundary_elements.push_back(find_or_insert_index_of(points_3, point, xvec, yvec, zvec) + 1);
                        }
                    }
                    ++quad_4_faces;
                    quad_4_face_list.push_back(++element_count);
                }
                else if (cell_[cx][cy][cz].type() == Cut) {
                    // For a cut cell, we want to write TRI_3 elements of the points which are on
                    // the plane which represents the boundary. Pick the three points out of any of
                    // the tetrahedra which meet this condition.
                    Nef_polyhedron::Plane_3 plane(Nef_polyhedron::Point_3(static_cast<double>(x), static_cast<double>(y), static_cast<double>(z)),
                                                  Nef_polyhedron::Vector_3(0, y == 0 ? -1 : 1, 0));
                    for (std::vector<Delaunay_triangulation>::iterator tri_it = T[cx][cy][cz].begin(); tri_it != T[cx][cy][cz].end(); ++tri_it) {
                        for (Delaunay_triangulation::Finite_cells_iterator it = tri_it->finite_cells_begin(); it != tri_it->finite_cells_end(); ++it) {
                            std::vector<Nef_polyhedron::Point_3> points;
                            for (int i = 0; i < 4; ++i) {
                                Nef_polyhedron::Point_3 point = it->vertex(i)->point();
                                if (plane.has_on(point)) {
                                    points.push_back(point);
                                }
                            }
                            if (points.size() == 3) {
                                boundary_elements.push_back(TRI_3);
                                for (std::vector<Nef_polyhedron::Point_3>::iterator pointit = points.begin(); pointit != points.end(); ++pointit) {
                                    boundary_elements.push_back(find_or_insert_index_of(points_3, *pointit, xvec, yvec, zvec) + 1); // indices are 1-based.
                                }
                                ++tri_3_faces;
                                tri_3_face_list.push_back(++element_count);
                            }
                        }
                    }
                }
            }
        }
        if (quad_4_faces > 0) {
            #ifndef NDEBUG
            std::cerr << "Writing " << (y == 0 ? "Jlo" : "Jhi") << " fluid boundary as list of size " << quad_4_face_list.size() << "." << std::endl;
            #endif
            if (cg_boco_write(index_file, index_base, index_zone, y == 0 ? "JloFluid" : "JhiFluid", BCWall, ElementList, quad_4_face_list.size(), &quad_4_face_list[0], &index_bc) != CG_OK) {
                std::cerr << cg_get_error() << std::endl;
                (void)cg_close(index_file);
                return 1;
            }
            quad_4_face_list.clear();
        }

        if (tri_3_faces > 0) {
            #ifndef NDEBUG
            std::cerr << "Writing " << (y == 0 ? "Jlo" : "Jhi") << " cut boundary as list of size " << tri_3_face_list.size() << "." << std::endl;
            #endif
            if (cg_boco_write(index_file, index_base, index_zone, y == 0 ? "JloCut" : "JhiCut", BCWall, ElementList, tri_3_face_list.size(), &tri_3_face_list[0], &index_bc) != CG_OK) {
                std::cerr << cg_get_error() << std::endl;
                (void)cg_close(index_file);
                return 1;
            }
            tri_3_face_list.clear();
        }

        tri_3_faces = 0;
        quad_4_faces = 0;
    }
    // Near/Far boundary
    for (V3NefIndex z = 0; z <= N_.shape()[2]; z += N_.shape()[2]) { // z = [0, NZ]
        for (V3NefIndex x = 0; x < N_.shape()[0]; ++x) {
            for (V3NefIndex y = 0; y < N_.shape()[1]; ++y) {
                const int cx = x, cy = y, cz = z > 0 ? z - 1 : z; // It's helpful to refer to the cell_ and T arrays, which have an index one less than the upper boundary.
                if (cell_[cx][cy][cz].type() == Fluid) {
                    boundary_elements.push_back(QUAD_4);
                    for (unsigned dx = 0; dx < 2; ++dx) {
                        for (unsigned dy = 0; dy < 2; ++dy) {
                            Nef_polyhedron::Point_3 point(static_cast<double>(x + dx), static_cast<double>(y + (dx == 0 ? (1 - dy) : dy)), static_cast<double>(z));
                            #ifndef NDEBUG
                            std::cerr << "Fluid cell boundary point at " << point << std::endl;
                            #endif
                            boundary_elements.push_back(find_or_insert_index_of(points_3, point, xvec, yvec, zvec) + 1);
                        }
                    }
                    ++quad_4_faces;
                    quad_4_face_list.push_back(++element_count);
                }
                else if (cell_[cx][cy][cz].type() == Cut) {
                    // For a cut cell, we want to write TRI_3 elements of the points which are on
                    // the plane which represents the boundary. Pick the three points out of any of
                    // the tetrahedra which meet this condition.
                    Nef_polyhedron::Plane_3 plane(Nef_polyhedron::Point_3(static_cast<double>(x), static_cast<double>(y), static_cast<double>(z)),
                                                  Nef_polyhedron::Vector_3(0, 0, z == 0 ? -1 : 1));
                    for (std::vector<Delaunay_triangulation>::iterator tri_it = T[cx][cy][cz].begin(); tri_it != T[cx][cy][cz].end(); ++tri_it) {
                        for (Delaunay_triangulation::Finite_cells_iterator it = tri_it->finite_cells_begin(); it != tri_it->finite_cells_end(); ++it) {
                            std::vector<Nef_polyhedron::Point_3> points;
                            for (int i = 0; i < 4; ++i) {
                                Nef_polyhedron::Point_3 point = it->vertex(i)->point();
                                if (plane.has_on(point)) {
                                    points.push_back(point);
                                }
                            }
                            if (points.size() == 3) {
                                boundary_elements.push_back(TRI_3);
                                for (std::vector<Nef_polyhedron::Point_3>::iterator pointit = points.begin(); pointit != points.end(); ++pointit) {
                                    boundary_elements.push_back(find_or_insert_index_of(points_3, *pointit, xvec, yvec, zvec) + 1); // indices are 1-based.
                                }
                                ++tri_3_faces;
                                tri_3_face_list.push_back(++element_count);
                            }
                        }
                    }
                }
            }
        }
        if (quad_4_faces > 0) {
            #ifndef NDEBUG
            std::cerr << "Writing " << (z == 0 ? "Klo" : "Khi") << " fluid boundary as list of size " << quad_4_face_list.size() << "." << std::endl;
            #endif
            if (cg_boco_write(index_file, index_base, index_zone, z == 0 ? "KloFluid" : "KhiFluid", BCWall, ElementList, quad_4_face_list.size(), &quad_4_face_list[0], &index_bc) != CG_OK) {
                std::cerr << cg_get_error() << std::endl;
                (void)cg_close(index_file);
                return 1;
            }
            quad_4_face_list.clear();
        }

        if (tri_3_faces > 0) {
            #ifndef NDEBUG
            std::cerr << "Writing " << (z == 0 ? "Klo" : "Khi") << " cut boundary as list of size " << tri_3_face_list.size() << "." << std::endl;
            #endif
            if (cg_boco_write(index_file, index_base, index_zone, z == 0 ? "KloCut" : "KhiCut", BCWall, ElementList, tri_3_face_list.size(), &tri_3_face_list[0], &index_bc) != CG_OK) {
                std::cerr << cg_get_error() << std::endl;
                (void)cg_close(index_file);
                return 1;
            }
            tri_3_face_list.clear();
        }

        tri_3_faces = 0;
        quad_4_faces = 0;
    }
    if (boundary_elements.size() > 0) {
        #ifndef NDEBUG
        std::cerr << "Writing " << (element_count - hexa_8_cells - tetra_4_cells) << " Fluid Boundary face elements with " << boundary_elements.size() << " nodes..." << std::endl;
        #endif
        if (cg_section_write(index_file, index_base, index_zone, "3_BoundaryElements", MIXED, hexa_8_cells + tetra_4_cells + 1, element_count /*- hexa_8_cells - tetra_4_cells*/, /* nbndry = */ 0, &boundary_elements[0], &index_section) != CG_OK) {
            std::cerr << cg_get_error() << std::endl;
            (void)cg_close(index_file);
            return 1;
        }
    }
    #ifndef NDEBUG
    std::cerr << "Pre-allocated " << hexa_8_capacity << " sizeof(int)s for HEXA_8, actually used " << hexa_8_elements.size() << "." << std::endl;
    std::cerr << "Pre-allocated " << tetra_4_capacity << " sizeof(int)s for TETRA_4, actually used " << tetra_4_elements.size() << "." << std::endl;
    std::cerr << "Pre-allocated " << boundary_capacity << " sizeof(int)s for boundary, actually used " << boundary_elements.size() << "." << std::endl;
    std::cerr << "Pre-allocated " << points_capacity << " sizeof(int)s for each points vector, actually used " << xvec.size() << "." << std::endl;
    #endif

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
