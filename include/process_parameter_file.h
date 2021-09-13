#ifndef _INCLUDE_PARAMETER_FILE_H_
#define _INCLUDE_PARAMETER_FILE_H_

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>

#include <deal.II/grid/cell_id.h>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>


namespace Elasticity
{
  using namespace dealii;

  struct Dimension
  {
    Dimension(const std::string &parameter_filename);

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);

    int dim;
  };


  template <int dim>
  struct GlobalParameters
  {
    GlobalParameters(const std::string &parameter_filename);
    GlobalParameters(
      const GlobalParameters<dim> &other); // This the the copy constructor

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);

    bool other_dirichlet_id;
    bool neumann_bc;
    bool use_E_and_nu;

    Point<dim> init_p1;
    Point<dim> init_p2;

    Point<dim> dirichlet_p1;
    Point<dim> dirichlet_p2;

    Point<dim> neumann_p1;
    Point<dim> neumann_p2;

    double E;
    double nu;

    double       mu;
    unsigned int mu_fr;

    double       lambda;
    unsigned int lambda_fr;

    double rho;

    double surface_force;
  };


  struct ParametersStd
  {
    ParametersStd(const std::string &parameter_filename);
    ParametersStd(const ParametersStd &other); // This the the copy constructor

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);

    bool verbose;
    bool direct_solver; /* This is often better for 2D problems. */

    unsigned int n_refine;
    unsigned int n_cycles;
  };


  struct ParametersMs
  {
    ParametersMs(const std::string &parameter_filename);
    ParametersMs(const ParametersMs &other); // This the the copy constructor

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);

    bool verbose;
    bool direct_solver; /* This is often better for 2D problems. */

    unsigned int n_refine;
  };


  struct ParametersBasis
  {
    ParametersBasis(const std::string &        parameter_filename,
                    const GlobalParameters<2> &global_parameters);
    ParametersBasis(const std::string &        parameter_filename,
                    const GlobalParameters<3> &global_parameters);
    ParametersBasis(
      const ParametersBasis &other); // This the the copy constructor

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);

    bool verbose;
    bool direct_solver; /* This is often better for 2D problems. */

    bool prevent_output;

    unsigned int n_refine;

    double       mu;
    unsigned int mu_fr;

    double       lambda;
    unsigned int lambda_fr;

    double rho;
  };


  extern template struct GlobalParameters<2>;
  extern template struct GlobalParameters<3>;
} // namespace Elasticity

#endif /* _INCLUDE_PARAMETER_FILE_H_ */