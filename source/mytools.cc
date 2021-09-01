#include "mytools.h"
#include "mytools.tpp"

namespace MyTools
{
  template const Point<2> roty(const Point<2> &p);
  template const Point<3> roty(const Point<3> &p);

  template const Point<2> backroty(const Point<2> &p);
  template const Point<3> backroty(const Point<3> &p);

  template void
  set_dirichlet_id(const Point<2> &p1,
                    const Point<2> &p2,
                    const unsigned int face_id,
                    const unsigned int id,
                    parallel::distributed::Triangulation<2> &triangulation);

  template void
  set_dirichlet_id(const Point<3> &p1,
                    const Point<3> &p2,
                    const unsigned int face_id,
                    const unsigned int id,
                    parallel::distributed::Triangulation<3> &triangulation);
}