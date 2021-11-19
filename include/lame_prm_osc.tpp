#ifndef _INCLUDE_LAME_PRM_OSC_TPP_
#define _INCLUDE_LAME_PRM_OSC_TPP_

#include "lame_prm_osc.h"

namespace Elasticity
{
  using namespace dealii;

  template <int dim>
  LamePrmOsc<dim>::LamePrmOsc(const double     &fr,
                              const double     &mean,
                              const Point<dim> &init_p1,
                              const Point<dim> &init_p2)
    : LamePrmBase<dim>()
    , p11(init_p1(0))
    , mean(mean)
  {
    two_pi_fr_l = 2 * M_PI * fr / (init_p2[dim - 1] - init_p1[dim - 1]);
  }

  template <int dim>
  double
  LamePrmOsc<dim>::value(const Point<dim> &p, const unsigned int) const
  {
    return mean * (0.8 * std::sin(two_pi_fr_l * (p(0) - p11)) + 1);
  }
} // namespace Elasticity

#endif // _INCLUDE_LAME_PRM_OSC_TPP_