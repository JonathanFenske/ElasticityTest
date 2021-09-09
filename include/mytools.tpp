#ifndef _INCLUDE_MY_TOOLS_TPP_
#define _INCLUDE_MY_TOOLS_TPP_

namespace MyTools
{
  using namespace dealii;

  /****************************************************************************/
  /* Tools */

  template <int dim>
  const Point<dim>
  roty(const Point<dim> &p)
  {
    Point<dim> q = p;
    if (dim == 3)
      {
        // Point<dim> axis = {0,1,0};
        // Tensor<2,dim> rot =
        // Physics::Transformations::Rotations::rotation_matrix_3d(axis,-M_PI);
        // return Point<dim>(rot*q);
        return Point<dim>(-q(2), q(1), q(0));
      }
    else
      {
        return q;
      }
  }


  template <int dim>
  const Point<dim>
  backroty(const Point<dim> &p)
  {
    Point<dim> q = p;
    if (dim == 3)
      {
        // Point<dim> axis = {0,1,0};
        // Tensor<2,dim> rot =
        // Physics::Transformations::Rotations::rotation_matrix_3d(axis,M_PI);
        // return Point<dim>(rot*q);
        return Point<dim>(q(2), q(1), -q(0));
      }
    else
      {
        return q;
      }
  }


  template <int dim>
  void
  set_dirichlet_id(const Point<dim> &p1,
                    const Point<dim> &p2,
                    const unsigned int face_id,
                    const unsigned int id,
                    parallel::distributed::Triangulation<dim> &triangulation)
  {
    for (auto &face : triangulation.active_face_iterators())
      {
        if (dim == 3)
          {
            if (face->at_boundary() && (face->boundary_id() == face_id) &&
                (face->center()[0] > p1[0]) && (face->center()[0] < p2[0]) &&
                (face->center()[1] > p1[1]) && (face->center()[1] < p2[1]) &&
                (face->center()[2] > p1[2]) && (face->center()[2] < p2[2]))
              face->set_boundary_id(id);
          }
        else if (face->at_boundary() && (face->boundary_id() == face_id) &&
                 (face->center()[0] > p1[0]) && (face->center()[0] < p2[0]) &&
                 (face->center()[1] > p1[1]) && (face->center()[1] < p2[1]))
          face->set_boundary_id(id);
      }
  }
}

#endif // _INCLUDE_MY_TOOLS_TPP_