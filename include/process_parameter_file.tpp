#ifndef _INCLUDE_PARAMETER_FILE_TPP_
#define _INCLUDE_PARAMETER_FILE_TPP_

#include <process_parameter_file.h>

/*!
 * Namespace for standard Q Multiscale problems.
 */
namespace Elasticity
{
  using namespace dealii;

  template <int dim>
  GlobalParameters<dim>::GlobalParameters(const std::string &parameter_filename)
  {
    ParameterHandler prm;

    declare_parameters(prm);

    std::ifstream parameter_file(parameter_filename);
    if (!parameter_file)
      {
        parameter_file.close();
        std::ofstream parameter_out(parameter_filename);
        prm.print_parameters(parameter_out, ParameterHandler::Text);
        AssertThrow(
          false,
          ExcMessage(
            "Input parameter file <" + parameter_filename +
            "> not found. Creating a template file of the same name."));
      }

    prm.parse_input(parameter_file,
                    /* filename = */ "generated_parameter.in",
                    /* last_line = */ "",
                    /* skip_undefined = */ true);
    parse_parameters(prm);
  }


  template <int dim>
  GlobalParameters<dim>::GlobalParameters(const GlobalParameters<dim> &other)
    : other_dirichlet_id(other.other_dirichlet_id)
    , neumann_bc(other.neumann_bc)
    , use_E_and_nu(other.use_E_and_nu)
    , init_p1(other.init_p1)
    , init_p2(other.init_p2)
    , dirichlet_p1(other.dirichlet_p1)
    , dirichlet_p2(other.dirichlet_p2)
    , neumann_p1(other.neumann_p1)
    , neumann_p2(other.neumann_p2)
    , E(other.E)
    , nu(other.nu)
    , mu(other.mu)
    , mu_fr(other.mu_fr)
    , lambda(other.lambda)
    , lambda_fr(other.lambda_fr)
    , rho(other.rho)
    , surface_force(other.surface_force)
  {}


  template <int dim>
  void
  GlobalParameters<dim>::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Global parameters");
    {
      prm.enter_subsection("Bools");
      {
        prm.declare_entry(
          "other dirichlet id",
          "false",
          Patterns::Bool(),
          "Choose whether to set the Dirichlet boundary condition to a place other than the default value.");
        prm.declare_entry(
          "neumann boundary condition",
          "false",
          Patterns::Bool(),
          "Choose whether a Neumann boundary condition is to be declared");
        prm.declare_entry("use E and nu",
                          "true",
                          Patterns::Bool(),
                          "Use Young's modulus and the Poisson ratio"
                          " to declare the Lamé parameters.");
      }
      prm.leave_subsection();

      std::string dim_str = std::to_string(dim) + "D";

      prm.enter_subsection(dim_str);
      {
        for (unsigned int i = 0; i < dim; ++i)
          {
            std::string tmp(std::string("init p1[") + std::to_string(i) + "]");
            prm.declare_entry(tmp,
                              "10",
                              Patterns::Double(),
                              "Declare the first vertex of the body.");

            prm.declare_entry(std::string("init p2[") + std::to_string(i) + "]",
                              "10",
                              Patterns::Double(),
                              "Declare the second vertex of the body.");

            prm.declare_entry(
              std::string("dirichlet p1[") + std::to_string(i) + "]",
              "10",
              Patterns::Double(),
              "Declare the first vertex of the rectangle,"
              " in which we apply the dirichlet boundary condition.");

            prm.declare_entry(
              std::string("dirichlet p2[") + std::to_string(i) + "]",
              "10",
              Patterns::Double(),
              "Declare the second vertex of the rectangle,"
              " in which we apply the dirichlet boundary condition.");

            prm.declare_entry(
              std::string("neumann p1[") + std::to_string(i) + "]",
              "10",
              Patterns::Double(),
              "Declare the first vertex of the rectangle,"
              " in which we apply the neumann boundary condition.");

            prm.declare_entry(
              std::string("neumann p2[") + std::to_string(i) + "]",
              "10",
              Patterns::Double(),
              "Declare the second vertex of the rectangle,"
              " in which we apply the dirichlet boundary condition.");
          }
      }
      prm.leave_subsection();

      prm.enter_subsection("Material parameters");
      {
        prm.declare_entry("E", "1", Patterns::Double(), "Set Young's modulus.");
        prm.declare_entry("nu",
                          "0.3",
                          Patterns::Double(),
                          "Set the Poisson ratio.");

        prm.declare_entry("mu",
                          "1.",
                          Patterns::Double(),
                          "Set the first Lamé parameter.");
        prm.declare_entry("mu frequency",
                          "0",
                          Patterns::Integer(),
                          "Set the frequency per body length of "
                          "the second Lamé parameter.");

        prm.declare_entry("lambda",
                          "1.",
                          Patterns::Double(),
                          "Set the first Lamé parameter.");
        prm.declare_entry("lambda frequency",
                          "0",
                          Patterns::Integer(),
                          "Set the frequency per body length of "
                          "the second Lamé parameter.");


        prm.declare_entry("rho",
                          "9.81",
                          Patterns::Double(),
                          "Set the mass density.");
      }
      prm.leave_subsection();

      prm.enter_subsection("Forces");
      {
        prm.declare_entry("surface force",
                          "1",
                          Patterns::Double(),
                          "Set the surface force.");
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }

  template <int dim>
  void
  GlobalParameters<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Global parameters");
    {
      prm.enter_subsection("Bools");
      {
        other_dirichlet_id = prm.get_bool("other dirichlet id");
        neumann_bc         = prm.get_bool("neumann boundary condition");
        use_E_and_nu       = prm.get_bool("use E and nu");
      }
      prm.leave_subsection();

      std::string dim_str = std::to_string(dim) + "D";

      prm.enter_subsection(dim_str);
      {
        for (unsigned int i = 0; i < dim; ++i)
          {
            init_p1[i] =
              prm.get_double(std::string("init p1[") + std::to_string(i) + "]");
            init_p2[i] =
              prm.get_double(std::string("init p2[") + std::to_string(i) + "]");

            dirichlet_p1[i] = prm.get_double(std::string("dirichlet p1[") +
                                             std::to_string(i) + "]");
            dirichlet_p2[i] = prm.get_double(std::string("dirichlet p2[") +
                                             std::to_string(i) + "]");

            neumann_p1[i] = prm.get_double(std::string("neumann p1[") +
                                           std::to_string(i) + "]");
            neumann_p2[i] = prm.get_double(std::string("neumann p2[") +
                                           std::to_string(i) + "]");
          }
      }
      prm.leave_subsection();

      prm.enter_subsection("Material parameters");
      {
        E  = prm.get_double("E");
        nu = prm.get_double("nu");

        if (use_E_and_nu)
          {
            mu     = E * nu / ((1 + nu) * (1 - 2 * nu));
            lambda = E / (2 * (1 + nu));
          }
        else
          {
            mu     = prm.get_double("mu");
            lambda = prm.get_double("lambda");
          }

        mu_fr = prm.get_integer("mu frequency");
        mu_fr = mu_fr / (init_p2[dim - 1] - init_p1[dim - 1]);

        lambda_fr = prm.get_integer("lambda frequency");
        lambda_fr = lambda_fr / (init_p2[dim - 1] - init_p1[dim - 1]);

        rho = prm.get_double("rho");
      }
      prm.leave_subsection();

      prm.enter_subsection("Forces");
      {
        surface_force = prm.get_double("surface force");
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }
} // namespace Elasticity

#endif /* _INCLUDE_PARAMETER_FILE_TPP_ */