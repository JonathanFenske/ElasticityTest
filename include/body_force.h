#ifndef _INCLUDE_BODY_FORCE_H_
#define _INCLUDE_BODY_FORCE_H_

#include <deal.II/base/function.h>

namespace Elasticity
{
  using namespace dealii;

  /**
   * @brief Class for the density of the body force in N/m^3
   *
   * @tparam dim Space dimension
   */
  template <int dim>
  class BodyForce : public Function<dim>
  {
  public:
    /**
     * @brief Construct a new Body Force object.
     *
     * @param rho Mass density
     */
    BodyForce(const double rho)
      : Function<dim>()
      , force_density(9.81 * rho)
    {}

    /**
     * @brief Returns one component of the body force at the point p.
     *
     * @param p Point
     * @param component Component from the output vector to be returned
     * @return double
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component) const override;

    /**
     * @brief Creates a list of body force values for one component.
     *
     * @param points Vector of points
     * @param values Vector of doubles to be overridden with body force values
     * @param component The body force component to be used.
     */
    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double>           &values,
               const unsigned int             component) const;

    /**
     * @brief Returns the body force at the point p.
     *
     * @param point Point
     * @param value Vector to be overridden with the body force.
     */
    virtual void
    vector_value(const Point<dim> &point, Vector<double> &value) const override;

    /**
     * @brief Creates a list of body force values.
     *
     * @param points Vector of points
     * @param value_list Vector of points to be overridden
     *                   with body force values
     */
    virtual void
    vector_value_list(const std::vector<Point<dim>> &points,
                      std::vector<Vector<double>>   &value_list) const override;

  private:
    /**
     * Body force density
     */
    const double force_density;
  };

  // exernal template instantiations
  extern template class BodyForce<2>;
  extern template class BodyForce<3>;
} // namespace Elasticity

#endif // _INCLUDE_BODY_FORCE_H_