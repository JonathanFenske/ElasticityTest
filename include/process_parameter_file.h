#ifndef _INCLUDE_PARAMETER_FILE_H_
#define _INCLUDE_PARAMETER_FILE_H_

#include <deal.II/base/exceptions.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>

#include <deal.II/grid/cell_id.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "include_lame_prm.h"
#include "mytools.h"


/**
 * @file process_parameter_file.h
 *
 * @brief Classes that can process parameter files.
 */


namespace Elasticity
{
  using namespace dealii;

  DeclExceptionMsg(
    ExcLayers,
    "The number of layers in each used direction must be at least 1.");

  /**
   * @brief Get the space dimesion from parameter files.
   */
  struct Dimension
  {
    /**
     * @brief Construct a new Dimension object.
     *
     * @param parameter_filename Filename of the parameter file
     *
     * This constructor creates the Dimension object and uses
     * declare_parameters() and parse_parameters() to
     * get the dimension #dim from a parameter file.
     */
    Dimension(const std::string &parameter_filename);

    /**
     * @brief Declare dimension
     *
     * @param prm ParameterHandler
     *
     * Declare dimension the parameter #dim for the <a href=
     * "https://www.dealii.org/9.2.0/doxygen/deal.II/
     * classParameterHandler.html">ParameterHandler</a>.
     */
    static void
    declare_parameters(ParameterHandler &prm);

    /**
     * @brief Parse the parameter
     *
     * @param prm ParameterHandler
     *
     * Parse the parameter #dim with the <a href=
     * "https://www.dealii.org/9.2.0/doxygen/deal.II/
     * classParameterHandler.html">ParameterHandler</a>.
     */
    void
    parse_parameters(ParameterHandler &prm);

    /**
     * Space dimension
     */
    int dim;
  };


  /**
   * @brief Collection of globally needed parameters
   *
   * @tparam dim Space dimension
   */
  template <int dim>
  struct ElaParameters
  {
    /**
     * @brief Construct a new (empty) ElaParameters object
     */
    ElaParameters()
    {}

    /**
     * @brief Construct a new ElaParameters object
     *
     * @param parameter_filename Path to parameter file
     *
     * This constructor creates the ElaParameters object and uses
     * declare_parameters() and parse_parameters() to
     * get the needed parameters from a parameter file.
     */
    ElaParameters(const std::string &parameter_filename);

    /**
     * @brief Copy constructor for ElaParameters
     *
     * @param other Other ElaParameters object
     */
    ElaParameters(const ElaParameters<dim> &other) = default;

    /**
     * @brief Declare parameters
     *
     * @param prm ParameterHandler
     *
     * Declare the needed parameters for the <a href=
     * "https://www.dealii.org/9.2.0/doxygen/deal.II/
     * classParameterHandler.html">ParameterHandler</a>.
     */
    static void
    declare_parameters(ParameterHandler &prm);

    /**
     * @brief Parse the parameters
     *
     * @param prm ParameterHandler
     *
     * Parse the needed parameters with the <a href=
     * "https://www.dealii.org/9.2.0/doxygen/deal.II/
     * classParameterHandler.html">ParameterHandler</a>.
     */
    void
    parse_parameters(ParameterHandler &prm);

    void
    parse_material_parameters(ParameterHandler &prm);

    /**
     * @brief First vertex of the body
     *
     * The first vertex that defines the body that is in this case
     * a <a href=
     * "https://www.dealii.org/9.2.0/doxygen/deal.II/
     * namespaceGridGenerator.html#ac76417d7404b75cf53c732f456e6e971">
     * subdivided_hyper_rectangle</a>.
     *
     * #init_p1 must be componentwise smaller than #init_p2.
     */
    Point<dim> init_p1;

    /**
     * @brief Second vertex of the body
     *
     * The second vertex that defines the body that is in this case
     * a <a href=
     * "https://www.dealii.org/9.2.0/doxygen/deal.II/
     * namespaceGridGenerator.html#ac76417d7404b75cf53c732f456e6e971">
     * subdivided_hyper_rectangle</a>.
     *
     * #init_p1 must be componentwise smaller than #init_p2.
     */
    Point<dim> init_p2;

    /**
     * First Lamé parameter
     */
    std::shared_ptr<LamePrmBase<dim>> lambda;

    /**
     * Second Lamé parameter/shear modulus
     */
    std::shared_ptr<LamePrmBase<dim>> mu;

    /**
     * Mass density of the body
     */
    double rho;

    /**
     * If true a direct solver will be used. This is often better
     * for 2D problems. If false an iterative solver will be used.
     */
    bool direct_solver_std;

    /**
     * If true a direct solver will be used. This is often better
     * for 2D problems. If false an iterative solver will be used.
     */
    bool direct_solver_ms;

    /**
     * If true a direct solver will be used. This is often better
     * for 2D problems. If false an iterative solver will be used.
     */
    bool direct_solver_basis;

    /**
     * Prevent output of basis functions.
     */
    bool prevent_basis_output;
    /**
     * Number of refinements on the coarse scale.
     */
    unsigned int coarse_refinements;

    /**
     * Number of refinements on the fine scale.
     */
    unsigned int fine_refinements;

    /**
     * True if one end shall be rotated.
     */
    bool rotate;

    /**
     * Rotation angle
     */
    double angle;
  };

  template <>
  void
  ElaParameters<2>::parse_material_parameters(ParameterHandler &prm);

  template <>
  void
  ElaParameters<3>::parse_material_parameters(ParameterHandler &prm);

  extern template struct ElaParameters<2>;
  extern template struct ElaParameters<3>;
} // namespace Elasticity

#endif /* _INCLUDE_PARAMETER_FILE_H_ */