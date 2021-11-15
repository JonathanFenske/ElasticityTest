#ifndef _INCLUDE_BODY_FORCE_TPP_
#define _INCLUDE_BODY_FORCE_TPP_

#include "body_force.h"

namespace Elasticity
{
  using namespace dealii;

  template <int dim>
  double
  BodyForce<dim>::value(const Point<dim> & /*p*/,
                        const unsigned int component) const
  {
    if (component == dim - 1)
      {
        return -force_density;
      }
    else
      {
        return 0.0;
      }
  }


  template <int dim>
  void
  BodyForce<dim>::value_list(const std::vector<Point<dim>> &points,
                             std::vector<double>           &values,
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
                values[i] = -force_density;
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
    value[dim - 1] = -force_density;
  }


  template <int dim>
  void
  BodyForce<dim>::vector_value_list(const std::vector<Point<dim>> &points,
                                    std::vector<Vector<double>>   &values) const
  {
    Assert(points.size() == values.size(),
           ExcDimensionMismatch(points.size(), values.size()));

    for (unsigned int i = 0; i < points.size(); ++i)
      {
        vector_value(points[i], values[i]);
      }
  }
} // namespace Elasticity

#endif // _INCLUDE_BODY_FORCE_TPP_