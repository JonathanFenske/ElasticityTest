#include "forces_and_parameters.h"
#include "forces_and_parameters.tpp"

namespace Elasticity
{
  template class SurfaceForce<2>;
  template class SurfaceForce<3>;

  template class BodyForce<2>;
  template class BodyForce<3>;

  template class mu<2>;
  template class mu<3>;

  template class lambda<2>;
  template class lambda<3>;
}