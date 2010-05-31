/* Intersection Test
 * 
 * This program generates two cubes and computes the intersection between
 * them using Release 3.6 of the CGAL library.
 */
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Polyhedron_VRML_1_ostream.h>

#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Nef_polyhedron::Vector_3 Vector_3;
typedef Nef_polyhedron::Aff_transformation_3 Aff_transformation_3;

int main() {
  Nef_polyhedron UnitCube, N1;
  Polyhedron P;

  std::ifstream in("cube.nef3"); // Load the definition of a unit cube from file
  in >> UnitCube;

  // Now copy the Unit Cube and translate it slightly
  N1 = UnitCube;
  Aff_transformation_3 transl(CGAL::TRANSLATION, Vector_3(0.5, 0.5, 0.5));
  N1.transform(transl);

  // The result is the intersection of the UnitCube and the offset cube
  Nef_polyhedron result = UnitCube - N1;

  // Write out the NEF Polyhedron as a OFF file
  result.convert_to_polyhedron(P);
  CGAL::VRML_1_ostream out( std::cout);
  out << P;

}
