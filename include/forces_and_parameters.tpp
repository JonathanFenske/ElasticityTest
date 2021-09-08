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
  SurfaceForce<dim>::value(const Point<dim> & /*p*/, const unsigned int) const
  {
    // if(p(0) > 4)
    //   return -100;//-100 *fabs(std::sin(M_PI*p(0)/5));
    // else
    //   return 0.;
    return 0;
  }



  template <int dim>
  double
  BodyForce<dim>::value(const Point<dim> & /*p*/,
                        const unsigned int component) const
  {
    if (component == dim - 1)
      {
        return -grav * rho;
      }
    else
      {
        return 0.0;
      }
  }


  template <int dim>
  void
  BodyForce<dim>::value_list(const std::vector<Point<dim>> &points,
                             std::vector<double> &          values,
                             const unsigned int             component) const
  {
    Assert(points.size() == values.size(),
           ExcDimensionMismatch(points.size(), values.size()));

    Assert(component <= dim, ExcDimensionMismatch(component, dim));

    if (component == dim - 1)
      {
        for (unsigned int i = 0; i < points.size(); ++i)
          {
            for (unsigned int j = 0; j < dim; ++j)
              {
                values[i] = -grav * rho;
              }
          }
      }
    else
      {
        std::fill(values.begin(), values.end(), 0.0);
      }
  }


  template <int dim>
  void
  BodyForce<dim>::vector_value(const Point<dim> & /*p*/,
                               Vector<double> &value) const
  {
    value          = 0.0;
    value[dim - 1] = -grav * rho;
  }


  template <int dim>
  void
  BodyForce<dim>::vector_value_list(const std::vector<Point<dim>> &points,
                                    std::vector<Vector<double>> &  values) const
  {
    Assert(points.size() == values.size(),
           ExcDimensionMismatch(points.size(), values.size()));

    for (unsigned int i = 0; i < points.size(); ++i)
      {
        vector_value(points[i], values[i]);
      }
  }



  template <int dim>
  double
  lambda<dim>::value(const Point<dim> & /*p*/, const unsigned int) const
  {
    // int fr = 80;
    return E / (2 * (1 + nu)); // * (std::sin(2 * fr * M_PI * p(0) / 20) + 1);
    // return -(std::sin(M_PI*p(0)/15)+2);//(-0.1 * p(0) + 2.5) * 1e9;////*;
  }


  template <int dim>
  double
  mu<dim>::value(const Point<dim> &p, const unsigned int) const
  {
    // int fr = 40;
    return E * nu / ((1 + nu) * (1 - 2 * nu)); // *
          //  (std::sin(2 * fr * M_PI * p(0) / 20) + 1);
    // return E * nu / ((1 + nu) * (1 - 2 * nu)) * (0.5 * std::sin(2 * fr * M_PI
    // * p(0) / 20) +1);//* std::sin(2 * fr * M_PI * p(1) / 20) + 1);
  }
} // namespace Elasticity

#endif // _INCLUDE_FORCES_AND_PARAMETERS_TPP_