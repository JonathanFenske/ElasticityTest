#include "lame_prm_base.h"
#include "lame_prm_const.h"
#include "lame_prm_const.tpp"
#include "lame_prm_layers.h"
#include "lame_prm_layers.tpp"
#include "lame_prm_osc.h"
#include "lame_prm_osc.tpp"

namespace Elasticity
{
  template class LamePrmBase<2>;
  template class LamePrmBase<3>;

  template class LamePrmConst<2>;
  template class LamePrmConst<3>;

  template class LamePrmLayers<2>;
  template class LamePrmLayers<3>;

  template class LamePrmOsc<2>;
  template class LamePrmOsc<3>;
} // namespace Elasticity