#ifndef _INCLUDE_PARAMETER_FILE_H_
#define _INCLUDE_PARAMETER_FILE_H_

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>

#include <deal.II/grid/cell_id.h>

#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "forces_and_lame_parameters.h"
#include "mytools.h"


/**
 * @file process_parameter_file.h
 *
 * @brief Classes that can process parameter files.
 */


namespace Elasticity
{
  using namespace dealii;

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
  struct GlobalParameters
  {
    /**
     * @brief Construct a new (empty) GlobalParameters object
     */
    GlobalParameters()
    {}

    /**
     * @brief Construct a new GlobalParameters object
     *
     * @param parameter_filename Path to parameter file
     *
     * This constructor creates the GlobalParameters object and uses
     * declare_parameters() and parse_parameters() to
     * get the needed parameters from a parameter file.
     */
    GlobalParameters(const std::string &parameter_filename);

    /**
     * @brief Copy constructor for GlobalParameters
     *
     * @param other Other GlobalParameters object
     */
    GlobalParameters(const GlobalParameters<dim> &other) = default;

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

    /**
     * @brief Reset face id for dirichlet boundary condition?
     *
     * True if the dirichlet boundary condition not only
     * applies to the face with the smallest x-values
     *
     * @see MyTools::set_dirichlet_id()
     */
    bool other_dirichlet_id;

    /**
     * True if a Neumann boundary condition are applied.
     */
    bool neumann_bc;

    bool rotate;

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
     * @brief First vertex of area where the Dirichlet boundary
     *        condition is applied.
     *
     *  The first vertex of area where the Dirichlet boundary
     *  condition is applied.
     *
     * #dirichlet_p1 must be componentwise smaller than #dirichlet_p2.
     */
    Point<dim> dirichlet_p1;

    /**
     * @brief Second vertex of area where the Dirichlet boundary
     *        condition is applied.
     *
     *  The second vertex of area where the Dirichlet boundary
     *  condition is applied.
     *
     * #dirichlet_p1 must be componentwise smaller than #dirichlet_p2.
     */
    Point<dim> dirichlet_p2;

    /**
     * @brief First vertex of area where the Neumann boundary
     *        condition is applied.
     *
     *  The first vertex of area where the Neumann boundary
     *  condition is applied.
     *
     * #neumann_p1 must be componentwise smaller than #neumann_p2.
     */
    Point<dim> neumann_p1;

    /**
     * @brief Second vertex of area where the Neumann boundary
     *        condition is applied.
     *
     *  The second vertex of area where the Neumann boundary
     *  condition is applied.
     *
     * #neumann_p1 must be componentwise smaller than #neumann_p2.
     */
    Point<dim> neumann_p2;

    /**
     * First Lamé parameter
     */
    LamePrm<dim> lambda;

    /**
     * Second Lamé parameter/shear modulus
     */
    LamePrm<dim> mu;

    /**
     * Mass density of the body
     */
    double rho;

    /**
     * Value of the surface force. Only has to be declared if
     * there is a Neumann boundary condition.
     */
    double surface_force;
  };


  /**
   * @brief Parameters needed for the class ElaStd
   */
  struct ParametersStd
  {
    /**
     * @brief Construct a new ParametersStd object
     *
     * @param parameter_filename Path to parameter file
     *
     * This constructor creates the ParametersStd object and uses
     * declare_parameters() and parse_parameters() to
     * get the needed parameters from a parameter file.
     */
    ParametersStd(const std::string &parameter_filename);

    /**
     * @brief Copy constructor for ParametersStd
     *
     * @param other ParametersStd
     */
    ParametersStd(const ParametersStd &other) = default;

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

    /**
     * Verbose
     */
    bool verbose;

    /**
     * If true a direct solver will be used. This is often better
     * for 2D problems. If false an iterative solver will be used.
     */
    bool direct_solver;

    /**
     * Number of refinements in the first cycle
     */
    unsigned int n_refine;

    /**
     * Number of cycles
     */
    unsigned int n_cycles;
  };


  /**
   * @brief Parameters needed for the class ElaMs
   *
   */
  struct ParametersMs
  {
    /**
     * @brief Construct a new ParametersStd object
     *
     * @param parameter_filename Path to parameter file
     *
     * This constructor creates the ParametersMs object and uses
     * declare_parameters() and parse_parameters() to
     * get the needed parameters from a parameter file.
     */
    ParametersMs(const std::string &parameter_filename);

    /**
     * @brief Copy constructor for ParametersMs
     *
     * @param other ParametersMs
     */
    ParametersMs(const ParametersMs &other) = default;

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

    /**
     * Verbose
     */
    bool verbose;

    /**
     * If true a direct solver will be used. This is often better
     * for 2D problems. If false an iterative solver will be used.
     */
    bool direct_solver;

    /**
     * Number of refinements on the coarse level
     */
    unsigned int n_refine;

    /**
     * Number of cycles
     */
    unsigned int n_cycles;
  };


  struct ParametersBasis
  {
    /**
     * @brief Construct a new (empty) ParametersBasis object
     */
    ParametersBasis()
    {}

    /**
     * @brief Construct a new ParametersBasis object
     *
     * @param parameter_filename Path to parameter file
     *
     * This constructor creates the ParametersBasis object and uses
     * declare_parameters() and parse_parameters() to
     * get the needed parameters from a parameter file.
     */
    ParametersBasis(const std::string &parameter_filename);

    /**
     * @brief Copy contructor for ParametersBasis
     *
     * @param other Other ParametersBasis object
     */
    ParametersBasis(const ParametersBasis &other) = default;

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

    /**
     * Verbose
     */
    bool verbose;

    /**
     * If true ElaBasis will use a direct solver. This is often better
     * for 2D problems. If false an iterative solver will be used.
     */
    bool direct_solver;

    /**
     * If true, the output of the locally constructed basis function
     * to vtu files is prevented.
     */
    bool prevent_output;

    /**
     * Number of refinements on the fine level
     */
    unsigned int n_refine;
  };


  extern template struct GlobalParameters<2>;
  extern template struct GlobalParameters<3>;
} // namespace Elasticity

#endif /* _INCLUDE_PARAMETER_FILE_H_ */