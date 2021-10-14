#ifndef _INCLUDE_FORCES_AND_LAME_PARAMETERS_TPP_
#define _INCLUDE_FORCES_AND_LAME_PARAMETERS_TPP_

#include "forces_and_lame_parameters.h"

namespace Elasticity
{
  using namespace dealii;

  /****************************************************************************/
  /* Forces and parameters */

  template <int dim>
  double
  SurfaceForce<dim>::value(const Point<dim> & /*p*/, const unsigned int) const
  {
    // if(p(0) > 4)
    //   return -100;//-100 *fabs(std::sin(M_PI*p(0)/5));
    // else
    //   return 0.;
    return force_value;
  }



  template <int dim>
  double
  BodyForce<dim>::value(const Point<dim> & /*p*/,
                        const unsigned int component) const
  {
    if (component == dim - 1)
      {
        return -force_value;
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
  LamePrm<dim>::LamePrm()
    : Function<dim>()
  {}

  template <int dim>
  LamePrm<dim>::LamePrm(const unsigned int &               n_x_layers,
                        const unsigned int &               n_y_layers,
                        const unsigned int &               n_z_layers,
                        const double &                     mean,
                        const std::vector<unsigned int> &  index_set,
                        const std::map<std::string, bool> &material_structure,
                        const Point<dim> &                 init_p1,
                        const Point<dim> &                 init_p2)
    : Function<dim>()
    , n_x_layers(n_x_layers)
    , n_y_layers(n_y_layers)
    , n_z_layers(n_z_layers)
    , material_structure(material_structure)
    , init_p1(init_p1)
    , init_p2(init_p2)
    , values(n_x_layers * n_y_layers * n_z_layers)
  {
    if (material_structure.at("horizontal layers") ||
        material_structure.at("vertical layers") ||
        material_structure.at("y-layers"))
      {
        unsigned int        n_values(n_x_layers * n_y_layers * n_z_layers);
        double              value_step_size = 2 * mean / (n_values + 1);
        std::vector<double> values_tmp;
        for (unsigned int i = 1; i < n_values + 1; ++i)
          {
            values_tmp.push_back(i * value_step_size);
          }

        for (unsigned int i = 0; i < n_values; ++i)
          {
            values[i] = values_tmp[index_set[i]];
          }

        if (material_structure.at("horizontal layers"))
          {
            double length    = init_p2(0) - init_p1(0);
            layer_size_inv_x = n_x_layers / length;
          }

        if (material_structure.at("vertical layers"))
          {
            double height    = init_p2(dim - 1) - init_p1(dim - 1);
            layer_size_inv_z = n_z_layers / height;
          }

        if (material_structure.at("y-layers"))
          {
            double depth     = init_p2(1) - init_p1(1);
            layer_size_inv_y = n_y_layers / depth;
          }
      }
  }

  template <int dim>
  LamePrm<dim>::LamePrm(const double &                     fr_tmp,
                        const double &                     mean,
                        const std::map<std::string, bool> &material_structure,
                        const Point<dim> &                 init_p1,
                        const Point<dim> &                 init_p2)
    : Function<dim>()
    , material_structure(material_structure)
    , init_p1(init_p1)
    , init_p2(init_p2)
    , values(1)
  {
    fr        = fr_tmp / (init_p2[dim - 1] - init_p1[dim - 1]);
    values[0] = mean;
  }


  template <int dim>
  double
  LamePrm<dim>::value(const Point<dim> &p, const unsigned int) const
  {
    // Case: oscillation in x-direction
    if (material_structure.at("oscillations"))
      return values[0] * (0.8 * std::sin(2 * fr * M_PI * p(0)) + 1);

    unsigned int x_layer = 0;
    unsigned int y_layer = 0;
    unsigned int z_layer = 0;

    // Case: horizontal layers
    if (material_structure.at("horizontal layers"))
      {
        x_layer = (p(0) - init_p1(0)) * layer_size_inv_x;
        x_layer = std::min(x_layer, n_x_layers - 1);
      }

    // Case: vertical layers
    if (material_structure.at("vertical layers"))
      {
        z_layer = (p(dim - 1) - init_p1(dim - 1)) * layer_size_inv_z;
        z_layer = std::min(z_layer, n_z_layers - 1);
      }

    // Case: layers in y-direction
    if (material_structure.at("y-layers"))
      {
        AssertDimension(dim, 3);
        y_layer = (p(1) - init_p1(1)) * layer_size_inv_y;
        y_layer = std::min(y_layer, n_y_layers - 1);
      }

    return (values[x_layer + n_x_layers * y_layer +
                   n_x_layers * n_y_layers * z_layer]);
  }
} // namespace Elasticity

#endif // _INCLUDE_FORCES_AND_LAME_PARAMETERS_TPP_