#ifndef _INCLUDE_LAME_PRM_OSC_H_
#define _INCLUDE_LAME_PRM_OSC_H_

#include "lame_prm_base.h"

namespace Elasticity
{
  using namespace dealii;

  /**
   * @brief Class for oscillating Lamé parameters in x-direction
   *
   * @tparam dim Space dimension
   */
  template <int dim>
  class LamePrmOsc : public LamePrmBase<dim>
  {
  public:
    /**
     * @brief Constructor
     *
     */
    LamePrmOsc(const double     &fr_tmp,
               const double     &mean,
               const Point<dim> &init_p1,
               const Point<dim> &init_p2);

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
    /**
     * 2 * pi * fr / length
     */
    double two_pi_fr_l;

    /**
     * First component of init_p1
     */
    double p11;

    /**
     * Mean value
     *
     */
    double mean;
  };

  // exernal template instantiations
  extern template class LamePrmOsc<2>;
  extern template class LamePrmOsc<3>;
} // namespace Elasticity

#endif // _INCLUDE_LAME_PRM_OSC_H_