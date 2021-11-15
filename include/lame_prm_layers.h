#ifndef _INCLUDE_LAME_PRM_LAYERS_H_
#define _INCLUDE_LAME_PRM_LAYERS_H_

#include "lame_prm_base.h"

namespace Elasticity
{
  using namespace dealii;

  /**
   * @brief Class for the Lamé parameters in layered materials
   *
   * @tparam dim Space dimension
   */
  template <int dim>
  class LamePrmLayers : public LamePrmBase<dim>
  {
  public:
    LamePrmLayers();
    /**
     * @brief Constructor
     */
    LamePrmLayers(const double                      &mean,
                  const std::vector<unsigned int>   &index_set,
                  const std::map<std::string, bool> &material_structure,
                  const Point<dim>                  &init_p1,
                  const Point<dim>                  &init_p2,
                  const unsigned int                &n_x_layers,
                  const unsigned int                &n_y_layers,
                  const unsigned int                &n_z_layers = 1);

    /**
     * @brief Returns the value of lambda at the point p.
     *
     * @param p Point
     * @param component Component from the output vector to be returned
     * @return double
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

  private:
    /**
     * Number of layers in x-,y- and in the
     * 3-dimensional case z-direction.
     */
    std::vector<unsigned int> n_layers;

    /**
     * The first entry is true if the material is layered in x-direction,
     * the second is true if it is layered in y-direction (and in the
     * 3-dimensional case the third entry is true if it is layered in
     * z-direction)
     */
    std::vector<bool> layers;

    /**
     * The inverse of the size of the layers in x-,y- and in the
     *  3-dimensional case z-direction.
     */
    std::vector<double> layer_size_inv;

    /**
     * @see ElaParameters::init_p1
     */
    Point<dim> init_p1;
    /**
     * @see ElaParameters::init_p2
     */
    Point<dim> init_p2;

    /**
     * Contains the values for each layer.
     */
    std::vector<double> values;
  };



  // declare specializations
  template <>
  LamePrmLayers<3>::LamePrmLayers(
    const double                      &mean,
    const std::vector<unsigned int>   &index_set,
    const std::map<std::string, bool> &material_structure,
    const Point<3>                    &init_p1,
    const Point<3>                    &init_p2,
    const unsigned int                &n_x_layers,
    const unsigned int                &n_y_layers,
    const unsigned int                &n_z_layers);

  template <>
  double
  LamePrmLayers<3>::value(const Point<3> &p, const unsigned int) const;

  template <>
  LamePrmLayers<2>::LamePrmLayers(
    const double                      &mean,
    const std::vector<unsigned int>   &index_set,
    const std::map<std::string, bool> &material_structure,
    const Point<2>                    &init_p1,
    const Point<2>                    &init_p2,
    const unsigned int                &n_x_layers,
    const unsigned int                &n_y_layers,
    const unsigned int                &n_z_layers);

  template <>
  double
  LamePrmLayers<2>::value(const Point<2> &p, const unsigned int) const;

  // exernal template instantiations
  extern template class LamePrmLayers<2>;
  extern template class LamePrmLayers<3>;
} // namespace Elasticity

#endif // _INCLUDE_LAME_PRM_LAYERS_H_