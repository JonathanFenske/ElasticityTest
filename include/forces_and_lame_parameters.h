#ifndef _INCLUDE_FORCES_AND_LAME_PARAMETERS_H_
#define _INCLUDE_FORCES_AND_LAME_PARAMETERS_H_

#include <deal.II/base/function.h>

#include <deal.II/fe/mapping_q1_eulerian.h>

#include <map>

#include "process_parameter_file.h"



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
    SurfaceForce(const double force_value)
      : Function<dim>()
      , force_value(force_value)
    {}

    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

  private:
    const double force_value;
  };


  // Class for the density of the surface force in N/m^2
  template <int dim>
  class BodyForce : public Function<dim>
  {
  public:
    BodyForce(const double rho)
      : Function<dim>()
      , rho(rho)
      , force_value(grav * rho)
    {}

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
    const double rho;

    const double force_value;
  };


  // Class for the first Lamé parameter of linear elasticity
  template <int dim>
  class lambda : public Function<dim>
  {
  public:
    lambda(const GlobalParameters<dim> &global_parameters);

    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

  private:
    GlobalParameters<dim> global_parameters;

    double layer_size_inv;

    double height;
    double length;

    std::vector<double> values;
  };


  // Class for the second Lamé parameter/shear modulus in linear elasticity
  template <int dim>
  class mu : public Function<dim>
  {
  public:
    mu(const GlobalParameters<dim> &global_parameters);

    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

  private:
    GlobalParameters<dim> global_parameters;

    double layer_size_inv;

    double height;
    double length;

    std::vector<double> values;
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

#endif // _INCLUDE_FORCES_AND_LAME_PARAMETERS_H_