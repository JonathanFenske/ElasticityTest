#ifndef _INCLUDE_MY_TOOLS_TPP_
#define _INCLUDE_MY_TOOLS_TPP_

#include <deal.II/base/symmetric_tensor.h>

#include <math.h>

#include "mytools.h"

namespace MyTools
{
  using namespace dealii;

  /****************************************************************************/
  /* Tools */

  template <int dim>
  const std::vector<unsigned int>
  get_repetitions(const Point<dim> &p1, const Point<dim> &p2)
  {
    std::list<double> side_length_list;
    for (unsigned int i = 0; i < dim; ++i)
      side_length_list.push_back(p2[i] - p1[i]);
    double cell_length =
      *(std::min_element(side_length_list.begin(), side_length_list.end()));
    std::vector<unsigned int> repetitions;
    for (auto it : side_length_list)
      repetitions.push_back((unsigned int)(it / cell_length));

    return repetitions;
  }


  template <int dim>
  Rotation<dim>::Rotation(const Point<dim> init_p1,
                          const Point<dim> init_p2,
                          const double     angle)
    : Function<dim>(dim)
    , y((init_p2(1) - init_p1(1)) / 2)
    , z((init_p2(2) - init_p1(2)) / 2)
    , R(Physics::Transformations::Rotations::rotation_matrix_3d({1, 0, 0},
                                                                angle))
  {
    AssertDimension(dim, 3);
  }

  template <int dim>
  void
  Rotation<dim>::vector_value(const Point<dim> &p, Vector<double> &values) const
  {
    // first translate the point such that the x-axis of the body is at the
    // origin and then apply the inverted translation to the rotated point
    Tensor<1, dim> tmp(p);
    tmp[1] -= y;
    tmp[2] -= z;
    tmp = R * tmp;
    tmp[1] += y;
    tmp[2] += z;
    tmp -= p;
    for (unsigned int i = 0; i < dim; ++i)
      values[i] = tmp[i];
  }
} // namespace MyTools

#endif // _INCLUDE_MY_TOOLS_TPP_