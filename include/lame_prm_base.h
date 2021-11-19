#ifndef _INCLUDE_LAME_PRM_H_
#define _INCLUDE_LAME_PRM_H_

#include <deal.II/base/function.h>

#include <deal.II/fe/mapping_q1_eulerian.h>

#include <map>

/**
 * @file LAME_PRM.h
 *
 * @brief Forces and parameters
 *
 */


namespace Elasticity
{
  using namespace dealii;

  /****************************************************************************/
  /* Forces and parameters */

  /**
   * @brief Base class for the Lamé parameters
   *
   * @tparam dim Space dimension
   */
  template <int dim>
  class LamePrmBase : public Function<dim>
  {
  public:
    LamePrmBase()
      : Function<dim>()
    {}

    LamePrmBase(const LamePrmBase<dim> &other) = default;

    /**
     * @brief Returns the value of lambda at the point p.
     *
     * @param p Point
     * @param component Component from the output vector to be returned
     * @return double
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const = 0;
  };

  // exernal template instantiations
  extern template class LamePrmBase<2>;
  extern template class LamePrmBase<3>;
} // namespace Elasticity

#endif // _INCLUDE_LAME_PRM_H_