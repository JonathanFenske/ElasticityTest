#ifndef _INCLUDE_LAME_PRM_LAYERS_TPP_
#define _INCLUDE_LAME_PRM_LAYERS_TPP_

#include "lame_prm_layers.h"

namespace Elasticity
{
  using namespace dealii;

  template <>
  LamePrmLayers<3>::LamePrmLayers(const double &                   mean,
                                  const std::vector<unsigned int> &index_set,
                                  const std::vector<bool> &        layers,
                                  const std::vector<unsigned int> &n_layers,
                                  const Point<3> &                 init_p1,
                                  const Point<3> &                 init_p2)
    : LamePrmBase<3>()
    , layers(layers)
    , n_layers(n_layers)
    , layer_size_inv(3)
    , init_p1(init_p1)
    , init_p2(init_p2)
    , values(n_layers[0] * n_layers[1] * n_layers[2])
  {
    unsigned int        n_values(n_layers[0] * n_layers[1] * n_layers[2]);
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

    if (layers[0])
      {
        double length     = init_p2(0) - init_p1(0);
        layer_size_inv[0] = n_layers[0] / length;
      }
    else
      {
        layer_size_inv[0] = 1;
      }


    if (layers[1])
      {
        double depth      = init_p2(1) - init_p1(1);
        layer_size_inv[1] = n_layers[1] / depth;
      }
    else
      {
        layer_size_inv[1] = 1;
      }

    if (layers[2])
      {
        double height     = init_p2(2) - init_p1(2);
        layer_size_inv[2] = n_layers[2] / height;
      }
    else
      {
        layer_size_inv[2] = 1;
      }
  }

  template <>
  double
  LamePrmLayers<3>::value(const Point<3> &p, const unsigned int) const
  {
    unsigned int x_layer = 0;
    unsigned int y_layer = 0;
    unsigned int z_layer = 0;

    // Case: x-layers
    if (layers[0])
      {
        x_layer = (p(0) - init_p1(0)) * layer_size_inv[0];
        x_layer = std::min(x_layer, n_layers[0] - 1);
      }

    // Case: y-layers
    if (layers[1])
      {
        y_layer = (p(1) - init_p1(1)) * layer_size_inv[1];
        y_layer = std::min(y_layer, n_layers[1] - 1);
      }

    // Case: z-layers
    if (layers[2])
      {
        z_layer = (p(2) - init_p1(2)) * layer_size_inv[2];
        z_layer = std::min(z_layer, n_layers[2] - 1);
      }

    return (values[x_layer + n_layers[0] * y_layer +
                   n_layers[0] * n_layers[1] * z_layer]);
  }

  template <>
  LamePrmLayers<2>::LamePrmLayers(const double &                   mean,
                                  const std::vector<unsigned int> &index_set,
                                  const std::vector<bool> &        layers,
                                  const std::vector<unsigned int> &n_layers,
                                  const Point<2> &                 init_p1,
                                  const Point<2> &                 init_p2)
    : LamePrmBase<2>()
    , layers(layers)
    , n_layers(n_layers)
    , layer_size_inv(2)
    , init_p1(init_p1)
    , init_p2(init_p2)
    , values(n_layers[0] * n_layers[1])
  {
    unsigned int        n_values(n_layers[0] * n_layers[1]);
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

    if (layers[0])
      {
        double length     = init_p2(0) - init_p1(0);
        layer_size_inv[0] = n_layers[0] / length;
      }
    else
      {
        layer_size_inv[0] = 1;
      }


    if (layers[1])
      {
        double depth      = init_p2(1) - init_p1(1);
        layer_size_inv[1] = n_layers[1] / depth;
      }
    else
      {
        layer_size_inv[1] = 1;
      }
  }

  template <>
  double
  LamePrmLayers<2>::value(const Point<2> &p, const unsigned int) const
  {
    unsigned int x_layer = 0;
    unsigned int y_layer = 0;

    // Case: x-layers
    if (layers[0])
      {
        x_layer = (p(0) - init_p1(0)) * layer_size_inv[0];
        x_layer = std::min(x_layer, n_layers[0] - 1);
      }

    // Case: y-layers
    if (layers[1])
      {
        y_layer = (p(1) - init_p1(1)) * layer_size_inv[1];
        y_layer = std::min(y_layer, n_layers[1] - 1);
      }

    return (values[x_layer + n_layers[0] * y_layer]);
  }
} // namespace Elasticity

#endif // _INCLUDE_LAME_PRM_LAYERS_TPP_