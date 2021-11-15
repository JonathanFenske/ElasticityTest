#ifndef _INCLUDE_LAME_PRM_CONST_TPP_
#define _INCLUDE_LAME_PRM_CONST_TPP_

#include "lame_prm_const.h"

namespace Elasticity
{
  using namespace dealii;

  template <int dim>
  LamePrmConst<dim>::LamePrmConst(const double &c_value)
    : LamePrmBase<dim>()
    , c_value(c_value)
  {}

  template <int dim>
  double
  LamePrmConst<dim>::value(const Point<dim> &p, const unsigned int) const
  {
    return c_value;
  }
} // namespace Elasticity

#endif // _INCLUDE_LAME_PRM_CONST_TPP_