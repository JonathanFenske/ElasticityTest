#ifndef _INCLUDE_FORCES_AND_LAME_PARAMETERS_TPP_
#define _INCLUDE_FORCES_AND_LAME_PARAMETERS_TPP_

#include "forces_and_lame_parameters.h"

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
  lambda<dim>::lambda(const GlobalParameters<dim> &global_parameters)
    : Function<dim>()
    , global_parameters(global_parameters)
  {
    if (global_parameters.material_structure.at("horizontal layers") ||
        global_parameters.material_structure.at("vertical layers") ||
        global_parameters.material_structure.at("y-layers"))
      {
        double value_step_size =
          global_parameters.lambda / (0.5 * global_parameters.n_layers);
        for (unsigned int i = 1; i < global_parameters.n_layers + 1; ++i)
          {
            values.push_back(i * value_step_size);
          }

        if (global_parameters.material_structure.at("horizontal layers"))
          {
            double length =
              global_parameters.init_p2(0) - global_parameters.init_p1(0);
            layer_size_inv = global_parameters.n_layers / length;
          }

        if (global_parameters.material_structure.at("vertical layers"))
          {
            double height = global_parameters.init_p2(dim - 1) -
                            global_parameters.init_p1(dim - 1);
            layer_size_inv = global_parameters.n_layers / height;
          }

        if (global_parameters.material_structure.at("y-layers"))
          {
            double depth =
              global_parameters.init_p2(1) - global_parameters.init_p1(1);
            layer_size_inv = global_parameters.n_layers / depth;
          }
      }
  }


  template <int dim>
  double
  lambda<dim>::value(const Point<dim> &p, const unsigned int) const
  {
    // Case: oscillation in x-direction
    if (global_parameters.material_structure.at("oscillations"))
      return global_parameters.lambda *
             (0.8 * std::sin(2 * global_parameters.lambda_fr * M_PI * p(0)) +
              1);

    // Case: horizontal layers
    if (global_parameters.material_structure.at("horizontal layers"))
      {
        unsigned int layer =
          (p(0) - global_parameters.init_p1(0)) * layer_size_inv;
        layer = std::min(layer, global_parameters.n_layers - 1);
        if (layer == global_parameters.n_layers)
          std::cout << layer << std::endl;
        return values[layer];
      }

    // Case: vertical layers
    if (global_parameters.material_structure.at("vertical layers"))
      {
        unsigned int layer =
          (p(dim - 1) - global_parameters.init_p1(dim - 1)) * layer_size_inv;
        layer = std::min(layer, global_parameters.n_layers - 1);
        return values[layer];
      }

    // Case: layers in y-direction
    if (global_parameters.material_structure.at("y-layers"))
      {
        AssertDimension(dim, 3);
        unsigned int layer =
          (p(1) - global_parameters.init_p1(1)) * layer_size_inv;
        layer = std::min(layer, global_parameters.n_layers - 1);
        return values[layer];
      }

    std::cout << "The material structure was not declared." << std::endl;
    exit(1);
  }


  template <int dim>
  mu<dim>::mu(const GlobalParameters<dim> &global_parameters)
    : Function<dim>()
    , global_parameters(global_parameters)
  {
    if (global_parameters.material_structure.at("horizontal layers") ||
        global_parameters.material_structure.at("vertical layers") ||
        global_parameters.material_structure.at("y-layers"))
      {
        double value_step_size =
          global_parameters.mu / (0.5 * global_parameters.n_layers);
        for (unsigned int i = 1; i < global_parameters.n_layers + 1; ++i)
          {
            values.push_back(i * value_step_size);
          }

        if (global_parameters.material_structure.at("horizontal layers"))
          {
            double length =
              global_parameters.init_p2(0) - global_parameters.init_p1(0);
            layer_size_inv = global_parameters.n_layers / length;
          }

        if (global_parameters.material_structure.at("vertical layers"))
          {
            double height = global_parameters.init_p2(dim - 1) -
                            global_parameters.init_p1(dim - 1);
            layer_size_inv = global_parameters.n_layers / height;
          }

        if (global_parameters.material_structure.at("y-layers"))
          {
            double depth =
              global_parameters.init_p2(1) - global_parameters.init_p1(1);
            layer_size_inv = global_parameters.n_layers / depth;
          }
      }
  }


  template <int dim>
  double
  mu<dim>::value(const Point<dim> &p, const unsigned int) const
  {
    // Case: oscillation in x-direction
    if (global_parameters.material_structure.at("oscillations"))
      return global_parameters.mu *
             (0.8 * std::sin(2 * global_parameters.mu_fr * M_PI * p(0)) + 1);

    // Case: horizontal layers
    if (global_parameters.material_structure.at("horizontal layers"))
      {
        unsigned int layer =
          (p(0) - global_parameters.init_p1(0)) * layer_size_inv;
        layer = std::min(layer, global_parameters.n_layers - 1);
        return values[layer];
      }

    // Case: vertical layers
    if (global_parameters.material_structure.at("vertical layers"))
      {
        unsigned int layer =
          (p(dim - 1) - global_parameters.init_p1(dim - 1)) * layer_size_inv;
        layer = std::min(layer, global_parameters.n_layers - 1);
        return values[layer];
      }

    // Case: layers in y-direction
    if (global_parameters.material_structure.at("y-layers"))
      {
        AssertDimension(dim, 3);
        unsigned int layer =
          (p(1) - global_parameters.init_p1(1)) * layer_size_inv;
        layer = std::min(layer, global_parameters.n_layers - 1);
        return values[layer];
      }

    std::cout << "The material structure was not declared." << std::endl;
    exit(1);
  }
} // namespace Elasticity

#endif // _INCLUDE_FORCES_AND_LAME_PARAMETERS_TPP_