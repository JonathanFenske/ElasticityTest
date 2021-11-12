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
        prm.declare_entry(
          "use E and nu",
          "true",
          Patterns::Bool(),
          "Choose whether to use Young's modulus and the Poisson ratio"
          " to declare Lamé's parameters");
        prm.declare_entry("oscillations",
                          "false",
                          Patterns::Bool(),
                          "Declare if the Lamé parameters oscillate.");
        prm.declare_entry("horizontal layers",
                          "false",
                          Patterns::Bool(),
                          "Declare if the material "
                          "consists of horizontal layers.");
        prm.declare_entry("vertical layers",
                          "false",
                          Patterns::Bool(),
                          "Declare if the material "
                          "consists of vertical layers.");
        prm.declare_entry("y-layers",
                          "false",
                          Patterns::Bool(),
                          "Declare if the material "
                          "is structured with layers in y-direction.");
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
        prm.declare_entry("n x layers",
                          "1",
                          Patterns::Integer(),
                          "Number of horizontal layers.");

        prm.declare_entry("n y layers",
                          "1",
                          Patterns::Integer(),
                          "Number of y layers.");

        prm.declare_entry("n z layers",
                          "1",
                          Patterns::Integer(),
                          "Number of vertical layers.");

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
                          Patterns::Double(),
                          "Set the frequency per body length of "
                          "the second Lamé parameter.");

        prm.declare_entry("lambda",
                          "1.",
                          Patterns::Double(),
                          "Set the first Lamé parameter.");
        prm.declare_entry("lambda frequency",
                          "0",
                          Patterns::Double(),
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

      prm.enter_subsection("Refinements");
      {
        prm.declare_entry("coarse refinements",
                          "1",
                          Patterns::Integer(),
                          "Set the number of refinements on the coarse scale.");
        prm.declare_entry("fine refinements",
                          "1",
                          Patterns::Integer(),
                          "Set the number of refinements on the fine scale.");
      }
      prm.leave_subsection();

      prm.enter_subsection("Rotations");
      {
        prm.declare_entry("rotate",
                          "false",
                          Patterns::Bool(),
                          "True if an end of the body shall be rotated.");
        prm.declare_entry("rotation angle",
                          "0",
                          Patterns::Double(),
                          "Set the rotation angle.");
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }

  template <int dim>
  void
  GlobalParameters<dim>::parse_parameters(ParameterHandler &prm)
  {
    bool                        use_E_and_nu;
    std::map<std::string, bool> material_structure;

    prm.enter_subsection("Global parameters");
    {
      prm.enter_subsection("Bools");
      {
        other_dirichlet_id = prm.get_bool("other dirichlet id");
        neumann_bc         = prm.get_bool("neumann boundary condition");

        // True if E and nu shall be used to declare mu and lambda.
        use_E_and_nu = prm.get_bool("use E and nu");

        int m = 0;
        material_structure.insert(
          std::make_pair("oscillations", prm.get_bool("oscillations")));

        material_structure.insert(
          std::make_pair("horizontal layers",
                         prm.get_bool("horizontal layers")));
        if (material_structure["horizontal layers"] == true)
          ++m;

        material_structure.insert(
          std::make_pair("vertical layers", prm.get_bool("vertical layers")));
        if (material_structure["vertical layers"] == true)
          ++m;

        material_structure.insert(
          std::make_pair("y-layers", prm.get_bool("y-layers")));
        if (material_structure["y-layers"] == true)
          ++m;

        if (material_structure["oscillations"] && m != 0)
          {
            std::cout << "The material can only be structured as either"
                      << "oscillations or layers but not both." << std::endl;
            exit(1);
          }
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
        // Mean value of the first Lamé parameter
        double lambda_mean;
        // Mean value of the second Lamé parameter
        double mu_mean;


        if (use_E_and_nu)
          {
            // Young's modulus/elastic modulus
            double E = prm.get_double("E");
            // Poisson ratio
            double nu = prm.get_double("nu");

            mu_mean     = E * nu / ((1 + nu) * (1 - 2 * nu));
            lambda_mean = E / (2 * (1 + nu));
          }
        else
          {
            mu_mean     = prm.get_double("mu");
            lambda_mean = prm.get_double("lambda");
          }

        double mu_fr = prm.get_double("mu frequency");

        double lambda_fr = prm.get_double("lambda frequency");

        if (material_structure.at("oscillations"))
          {
            lambda = LamePrm<dim>(
              lambda_fr, lambda_mean, material_structure, init_p1, init_p2);

            mu = LamePrm<dim>(
              mu_fr, mu_mean, material_structure, init_p1, init_p2);
          }
        else
          {
            unsigned int n_x_layers = prm.get_integer("n x layers");
            unsigned int n_y_layers = prm.get_integer("n y layers");
            unsigned int n_z_layers = prm.get_integer("n z layers");
            unsigned int n_layers   = n_x_layers * n_y_layers * n_z_layers;

            std::vector<unsigned int> index_set(n_layers);
            std::iota(index_set.begin(), index_set.end(), 0);
            std::seed_seq seq{1, 2, 3, 4, 5};
            std::mt19937  rd(seq);
            std::shuffle(index_set.begin(), index_set.end(), rd);

            lambda = LamePrm<dim>(n_x_layers,
                                  n_y_layers,
                                  n_z_layers,
                                  lambda_mean,
                                  index_set,
                                  material_structure,
                                  init_p1,
                                  init_p2);

            mu = LamePrm<dim>(n_x_layers,
                              n_y_layers,
                              n_z_layers,
                              mu_mean,
                              index_set,
                              material_structure,
                              init_p1,
                              init_p2);
          }

        rho = prm.get_double("rho");
      }
      prm.leave_subsection();

      prm.enter_subsection("Forces");
      {
        surface_force = prm.get_double("surface force");
      }
      prm.leave_subsection();

      prm.enter_subsection("Refinements");
      {
        coarse_refinements = prm.get_integer("coarse refinements");
        fine_refinements   = prm.get_integer("fine refinements");
      }
      prm.leave_subsection();

      prm.enter_subsection("Rotations");
      {
        rotate = prm.get_bool("rotate");
        angle  = prm.get_double("rotation angle") * M_PI;
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }
} // namespace Elasticity

#endif /* _INCLUDE_PARAMETER_FILE_TPP_ */