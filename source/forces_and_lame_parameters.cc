#include "forces_and_lame_parameters.h"

#include "forces_and_lame_parameters.tpp"

namespace Elasticity
{
  template class SurfaceForce<2>;
  template class SurfaceForce<3>;

  template class BodyForce<2>;
  template class BodyForce<3>;

  template class LamePrm<2>;
  template class LamePrm<3>;
} // namespace Elasticity