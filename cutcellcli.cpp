/* Cartesian cut cell grid generator: Command line interface
 *
 * Given an OFF file and grid parameters, cut a grid and output the CGNS
 * representation of the resulting cartesian cut cell grid.
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
#include <boost/program_options.hpp>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <iostream>
#include <fstream>
#include <string>
#include <string.h> // For strerror.
#include <errno.h>

int main(int argc, char **argv) {
    int NX, NY, NZ;
    double dx, dy, dz;
    double x0, y0, z0;
    double oX, oY, oZ, s;
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("size-x,X", boost::program_options::value<int>(&NX), "grid size in X dimension")
        ("size-y,Y", boost::program_options::value<int>(&NY), "grid size in Y dimension")
        ("size-z,Z", boost::program_options::value<int>(&NZ), "grid size in Z dimension")
        ("dx", boost::program_options::value<double>(&dx)->default_value(1.0), "cell size in X dimension")
        ("dy", boost::program_options::value<double>(&dy)->default_value(1.0), "cell size in Y dimension")
        ("dz", boost::program_options::value<double>(&dz)->default_value(1.0), "cell size in Z dimension")
        ("x0", boost::program_options::value<double>(&x0)->default_value(0.0), "X component of origin")
        ("y0", boost::program_options::value<double>(&y0)->default_value(0.0), "Y component of origin")
        ("z0", boost::program_options::value<double>(&z0)->default_value(0.0), "Z component of origin")
        ("offset-x,x", boost::program_options::value<double>(&oX)->default_value(0.0), "translation vector in X dimension")
        ("offset-y,y", boost::program_options::value<double>(&oY)->default_value(0.0), "translation vector in Y dimension")
        ("offset-z,z", boost::program_options::value<double>(&oZ)->default_value(0.0), "translation vector in Z dimension")
        ("scaling,s", boost::program_options::value<double>(&s)->default_value(1.0), "scale factor")
        ("input-file,I", boost::program_options::value<std::string>(), "Input file in OFF format (omit to read from stdin)")
        ("output-file,O", boost::program_options::value<std::string>(), "Output file in CGNS format (omit to write to stdout)")
    ;
    boost::program_options::positional_options_description p;
    p.add("input-file", -1);
    p.add("output-file", -1);
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
    boost::program_options::notify(vm);

    if (vm.count("help") || !vm.count("size-x") || !vm.count("size-y") || !vm.count("size-z")) {
        std::cout << desc << "\n";
        return 1;
    }

    cutcell::Polyhedron P;

    bool badbit = false;
    // Read the OFF format from the input file or stdin.
    if (vm.count("input-file")) {
        std::ifstream in(vm["input-file"].as<std::string>().c_str());
        if (!in) {
            std::cerr << "Failed to open input file. Unreliable error message: " << strerror(errno) << std::endl;
            badbit = true;
        }
        else {
            CGAL::scan_OFF(in, P, true);
            badbit = in.bad();
        }
    } else {
        CGAL::scan_OFF(std::cin, P, true);
        badbit = std::cin.bad();
    }
    if (badbit) {
        std::cerr << "Input file does not contain a permissible polyhedral surface." << std::endl;
        return 1;
    }

    // Convert Polyhedron to Nef_polyhedron
    cutcell::Nef_polyhedron N1(P);

    // Transform the object accordingly.
    cutcell::Aff_transformation Aff1(cutcell::TRANSLATION, cutcell::Vector(oX, oY, oZ));
    cutcell::Aff_transformation Aff2(cutcell::SCALING, s);
    N1.transform(Aff1);
    N1.transform(Aff2);

    // Create the cutting grid.
    cutcell::Grid grid(NX, NY, NZ);

    // Register the cell size and origin transform
    grid.addTransform(cutcell::Aff_transformation(dx, 0,  0,  x0,
                                                 0,  dy, 0,  y0,
                                                 0,  0,  dz, z0));

    // Cut the grid to the object.
    grid.addSolid(N1);
    grid.cut();

    // Output the grid in CGNS format.
    if (vm.count("output-file"))
        grid.output_cgns_file(vm["output-file"].as<std::string>());
    else
        std::cout << cutcell::GridFormat(cutcell::GridFormat::OUTPUT_CGNS) << grid << std::endl;

    return 0;
}
