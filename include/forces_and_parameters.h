#ifndef _INCLUDE_FORCES_AND_PARAMETERS_H_
#define _INCLUDE_FORCES_AND_PARAMETERS_H_

#include <deal.II/base/function.h>

namespace Elasticity
{
  using namespace dealii;
  
  /****************************************************************************/
  /* Parameters and boundary conditions */

  // Currently parameter values for steel. Note that the mass density that is
  // needed for the force density of the body forces must also modified if
  // another material is simulated.
  const double E  = 210.e9; // Young modulus
  const double nu = 0.3;    // Poisson ratio

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
      vector_value_list(
        const std::vector<Point<dim>> &       points,
        std::vector<Vector<double>> &value_list) const override;
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
}

#endif // _INCLUDE_FORCES_AND_PARAMETERS_H_