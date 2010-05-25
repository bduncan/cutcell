/* Simple Polyhedron
 * 
 * This program generates a simple polyhedron using the make_cube_3
 * function from the CGAL examples. This creates unit cubes with the
 * bottom left corner as the origin which are examples of the
 * Polyhedron_3 class.
 *
 * The main routine converts the cube from a simple Polyhedron to a
 * NEF_Polyhedron and then applies a simple translation before
 * outputing the resultin polyhedron as a .off file.
 *
 * This program uses the exact_predicates_exact_constructions kernal
 * from CGAL which is recomended for manipulations using <double>
 * coordinates and is efficient for geometric operations.
 */
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <iostream>

template <class Poly>
typename Poly::Halfedge_handle make_cube_3( Poly& P) {
  // appends a unit cube of size [0,1]^3 to the polyhedron P. The cube
  // is constructed using Euler operations and this routine is taken
  // from the example code in section 25.3.7 of the CGAL Manual.
  CGAL_precondition( P.is_valid());
  typedef typename Poly::Point_3 Point;
  typedef typename Poly::Halfedge_handle Halfedge_handle;
  Halfedge_handle h = P.make_tetrahedron( Point( 1, 0, 0),
					  Point( 0, 0, 1),
					  Point( 0, 0, 0),
					  Point( 0, 1, 0));

  Halfedge_handle g = h->next()->opposite()->next(); 
  P.split_edge( h->next());
  P.split_edge( g->next());
 
  P.split_edge( g);
  h->next()->vertex()->point() = Point( 1, 0, 1);
  g->next()->vertex()->point() = Point( 0, 1, 1);

  g->opposite()->vertex()->point() = Point( 1, 1, 0);

  Halfedge_handle f = P.split_facet( g->next(),
				     g->next()->next()->next());

  Halfedge_handle e = P.split_edge( f);
  e->vertex()->point() = Point( 1, 1, 1);

  P.split_facet( e, f->next()->next()); 
  CGAL_postcondition( P.is_valid());
  return h;
}

/* This is the start of the main program.  First we declaire the CGAL
   kernal we wish to use. In this case we will use the exact
   predicates kernel. The Exact predicates exact constructions kernel
   uses filtering. In non-degenerate scenarios it's faster than the
   Homogeneous kernel. The most important advantage of the filtered
   kernel is that it is a Cartesian kernel, which allows the proper
   handling of OFF files using floating-point coordinates. */

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

// Now create an example of the Polyhedron_3 class, together with the
// half edge handle. These are needed to create cubes.
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Polyhedron::Halfedge_handle Halfedge_handle;

// Now create an example of the Nef_Polygon_3 class
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;

// And finally the Vector_3 and Aff_transformation_3 classes are
// needed to be able to perform the Affine transformation to be
// applied to the unit cube.
typedef Nef_polyhedron::Vector_3 Vector_3;
typedef Nef_polyhedron::Aff_transformation_3 Aff_transformation_3;

int main() {
  Polyhedron P; // P is a polyhedron
  Halfedge_handle h = make_cube_3( P); // make a unit cube 

  Nef_polyhedron N(P); // construct N a Nef_Polyhedron from P

  // Apply a translation of (0.5,0.5,0.5) to the cube
  Aff_transformation_3 Aff(CGAL::TRANSLATION, Vector_3(0.5, 0.5, 0.5));
  N.transform(Aff);

  // Now check to see if it is a simple polyhedron and if it is
  // convert it back to the Polyhedron_3 class and output it, this
  // will write the result out in OFF format.
  if(N.is_simple()) {
    N.convert_to_polyhedron(P);
    std::cout << P;
  }
  else {
    std::cout << N;
  }

}
