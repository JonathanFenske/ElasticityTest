#ifndef _INCLUDE_LAME_PRM_CONST_H_
#define _INCLUDE_LAME_PRM_CONST_H_

#include "lame_prm_base.h"

namespace Elasticity
{
  using namespace dealii;

  /**
   * @brief Class for constant Lamé parameters
   *
   * @tparam dim Space dimension
   */
  template <int dim>
  class LamePrmConst : public LamePrmBase<dim>
  {
  public:
    /**
     * @brief Constructor
     *
     */
    LamePrmConst(const double &c_value);

    /**
     * @brief Returns the value of lambda at the point p.
     *
     * @param p Point
     * @param component Component from the output vector to be returned
     * @return double
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

  private:
    const double c_value;
  };

  // exernal template instantiations
  extern template class LamePrmConst<2>;
  extern template class LamePrmConst<3>;
} // namespace Elasticity

#endif // _INCLUDE_LAME_PRM_CONST_H_