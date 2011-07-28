/* Cartesian cut cell grid generator: Test driver
 *
 * Cut a test cube at (0.5, 0.5, 0.5) and output the CGNS representation
 * of the resulting cartesian cut cell grid.
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

int main() {
    const int NX = 5, NY = 5, NZ = 5;
    cutcell::Nef_polyhedron N1;

    std::istringstream in(cutcell::cube); // Load the definition of a unit cube from string
    in >> N1;                    // into the test cube
    assert(N1.number_of_vertices() == 8);
    assert(N1.number_of_facets() == 6);
    assert(N1.number_of_edges() == 12);
    assert(N1.number_of_volumes() == 2);

    // Put the test cube at an appropriate place
    cutcell::Aff_transformation Aff1(cutcell::TRANSLATION, cutcell::Vector(0.5, 0.5, 0.5));
    cutcell::Aff_transformation Aff2(cutcell::SCALING, 1.5);
    N1.transform(Aff1);
    N1.transform(Aff2);

    // Create the cutting object.
    cutcell::Grid grid(NX, NY, NZ);

    // Cut the grid to the test cube.
    grid.addSolid(N1);
    grid.cut();

    // Output the grid in NEF format.
    std::cout << cutcell::GridFormat(cutcell::GridFormat::OUTPUT_CGNS) << grid << std::endl;

    return 0;
}
