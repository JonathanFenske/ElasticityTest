#include "postprocessing.h"
#include "postprocessing.tpp"

namespace Elasticity
{
  /****************************************************************************/
  /* Postprocessing */

  template class StrainPostprocessor<2>;
  template class StrainPostprocessor<3>;

  template class StressPostprocessor<2>;
  template class StressPostprocessor<3>;
}