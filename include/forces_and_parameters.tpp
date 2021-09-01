#ifndef _INCLUDE_FORCES_AND_PARAMETERS_TPP_
#define _INCLUDE_FORCES_AND_PARAMETERS_TPP_

#include "forces_and_parameters.h"

namespace Elasticity
{
  using namespace dealii;
  
  /****************************************************************************/
  /* Forces and parameters */

  // Currently parameter values for steel. Note that the mass density that is
  // needed for the force density of the body forces must also modified if
  // another material is simulated.
  const double E  = 210.e9; // Young modulus
  const double nu = 0.3;    // Poisson ratio


  template <int dim>
  double
  SurfaceForce<dim>::value(const Point<dim> &p, const unsigned int) const
  {
    // if(p(0) > 4)
    //   return -100;//-100 *fabs(std::sin(M_PI*p(0)/5));
    // else
    //   return 0.;
    return 0;
  }


  template <int dim>
  double
  BodyForce<dim>::value(const Point<dim> &p, const unsigned int component) const
  {   
    const double grav = 9.81; // Gravitational acceleration
    const double rho = 7.85e3; // Mass density (of steel)
    return ((component == dim - 1) ? (-grav * rho) : 0.);
  }


  template <int dim>
  void
  BodyForce<dim>::vector_value_list(
    const std::vector<Point<dim>> &       points,
    std::vector<Vector<double>> &value_list) const
  {
    const unsigned int n_points(points.size());
    for (unsigned int i = 0; i < n_points; ++i)
      for (unsigned int j = 0; j < dim; ++j)
        value_list[i][j] = value(points[i],j);
  }


  template <int dim>
  double
  lambda<dim>::value(const Point<dim> &p, const unsigned int) const
  {
    //int fr = 80;
    return E / (2 * (1 + nu)); // * (std::sin(2 * fr * M_PI * p(0) / 20) + 1);
    // return -(std::sin(M_PI*p(0)/15)+2);//(-0.1 * p(0) + 2.5) * 1e9;////*;
  }


  template <int dim>
  double
  mu<dim>::value(const Point<dim> &p, const unsigned int) const
  {
    int fr = 40;
    return E * nu / ((1 + nu) * (1 - 2 * nu)) *
           (std::sin(2 * fr * M_PI * p(0) / 20) + 1);
    // return E * nu / ((1 + nu) * (1 - 2 * nu)) * (0.5 * std::sin(2 * fr * M_PI
    // * p(0) / 20) +1);//* std::sin(2 * fr * M_PI * p(1) / 20) + 1);
  }
}

#endif // _INCLUDE_FORCES_AND_PARAMETERS_TPP_