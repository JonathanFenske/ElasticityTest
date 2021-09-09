#ifndef _INCLUDE_FORCES_AND_PARAMETERS_H_
#define _INCLUDE_FORCES_AND_PARAMETERS_H_

#include <deal.II/base/function.h>

#include <deal.II/fe/mapping_q1_eulerian.h>



namespace Elasticity
{
  using namespace dealii;

  /****************************************************************************/
  /* Forces and parameters */

  // Class for the density of the surface force in N/m^2
  template <int dim>
  class SurfaceForce : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;
  };


  // Class for the density of the surface force in N/m^2
  template <int dim>
  class BodyForce : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> &p, const unsigned int component) const override;

    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double> &          values,
               const unsigned int             component) const;

    virtual void
    vector_value(const Point<dim> &point, Vector<double> &value) const override;

    virtual void
    vector_value_list(const std::vector<Point<dim>> &points,
                      std::vector<Vector<double>> &  value_list) const override;

  private:
    /*!
     * Gravitational acceleration
     */
    const double grav = 9.81;

    /*!
     * Mass density (of steel)
     */
    const double rho = 7.85e3;
  };


  // Class for the first Lamé parameter of linear elasticity
  template <int dim>
  class lambda : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;
  };


  // Class for the second Lamé parameter/shear modulus in linear elasticity
  template <int dim>
  class mu : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;
  };

  // exernal template instantiations
  extern template class BodyForce<2>;
  extern template class SurfaceForce<2>;
  extern template class lambda<2>;
  extern template class mu<2>;

  extern template class BodyForce<3>;
  extern template class SurfaceForce<3>;
  extern template class lambda<3>;
  extern template class mu<3>;
} // namespace Elasticity

#endif // _INCLUDE_FORCES_AND_PARAMETERS_H_