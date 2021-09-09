#ifndef _INCLUDE_MY_TOOLS_H_
#define _INCLUDE_MY_TOOLS_H_

#include <deal.II/base/point.h>
#include <deal.II/distributed/tria.h>

namespace MyTools
{
  using namespace dealii;

  /****************************************************************************/
  /* Tools */

  // Function that rotates a point -90° or -(pi/2) around the y axis
  template <int dim>
  const Point<dim>
    roty(const Point<dim> &p);


  // Function that rotates a point +90° or +(pi/2) around the y axis
  template <int dim>
  const Point<dim>
    backroty(const Point<dim> &p);


  // Function that gives a specific subset of the surface a new face_id
  template <int dim>
  void set_dirichlet_id(const Point<dim> &p1,
                          const Point<dim> &p2,
                          const unsigned int face_id,
                          const unsigned int id,
                          parallel::distributed::Triangulation<dim> &triangulation);
}

#endif // _INCLUDE_MY_TOOLS_H_