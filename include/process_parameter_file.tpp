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
  ElaParameters<dim>::ElaParameters(const std::string &parameter_filename)

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
  ElaParameters<dim>::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Elasticity Parameters");
    {
      std::string dim_str = std::to_string(dim) + "D";

      prm.enter_subsection(dim_str);
      {
        std::string tmp;
        for (unsigned int i = 0; i < dim; ++i)
          {
            tmp = std::string("init p1[") + std::to_string(i) + ']';
            prm.declare_entry(tmp,
                              "10",
                              Patterns::Double(),
                              "Declare the first vertex of the body.",
                              /* has_to_be_set = */ true);

            tmp = std::string("init p2[") + std::to_string(i) + "]";
            prm.declare_entry(tmp,
                              "10",
                              Patterns::Double(),
                              "Declare the second vertex of the body.",
                              /* has_to_be_set = */ true);
          }
      }
      prm.leave_subsection();

      prm.enter_subsection("Material Parameters");
      {
        prm.declare_entry(
          "use E and nu",
          "true",
          Patterns::Bool(),
          "Choose whether to use Young's modulus and the Poisson ratio"
          " to declare Lamé's parameters",
          /* has_to_be_set = */ true);

        prm.declare_entry("oscillations",
                          "false",
                          Patterns::Bool(),
                          "Declare if the Lamé parameters oscillate.",
                          /* has_to_be_set = */ true);
        prm.declare_entry("x-layers",
                          "false",
                          Patterns::Bool(),
                          "Declare if the material "
                          "is structured with layers in x-direction.",
                          /* has_to_be_set = */ true);
        prm.declare_entry("y-layers",
                          "false",
                          Patterns::Bool(),
                          "Declare if the material "
                          "is structured with layers in y-direction.",
                          /* has_to_be_set = */ true);
        prm.declare_entry("z-layers",
                          "false",
                          Patterns::Bool(),
                          "Declare if the material "
                          "is structured with layers in z-direction.");

        prm.declare_entry("no. of x-layers",
                          "1",
                          Patterns::Integer(),
                          "Number of x-layers.");
        prm.declare_entry("no. of y-layers",
                          "1",
                          Patterns::Integer(),
                          "Number of y-layers.");
        prm.declare_entry("no. of z-layers",
                          "1",
                          Patterns::Integer(),
                          "Number of z-layers.");

        prm.declare_entry("E", "1", Patterns::Double(), "Set Young's modulus.");
        prm.declare_entry("nu",
                          "0.3",
                          Patterns::Double(),
                          "Set the Poisson ratio.");

        prm.declare_entry("mu",
                          "1.",
                          Patterns::Double(),
                          "Set the first Lamé parameter.");
        prm.declare_entry("lambda",
                          "1.",
                          Patterns::Double(),
                          "Set the first Lamé parameter.");

        prm.declare_entry("mu frequency",
                          "0",
                          Patterns::Double(),
                          "Set the frequency per body length of "
                          "the second Lamé parameter.");
        prm.declare_entry("lambda frequency",
                          "0",
                          Patterns::Double(),
                          "Set the frequency per body length of "
                          "the second Lamé parameter.");


        prm.declare_entry("rho",
                          "1.",
                          Patterns::Double(),
                          "Set the mass density.",
                          /* has_to_be_set = */ true);
      }
      prm.leave_subsection();

      prm.enter_subsection("Other Parameters");
      {
        prm.declare_entry("verbose",
                          "true",
                          Patterns::Bool(),
                          "Verbose",
                          /* has_to_be_set = */ true);

        prm.declare_entry(
          "direct solver std",
          "false",
          Patterns::Bool(),
          "Declare if a direct solver shall be used for the standard FEM",
          /* has_to_be_set = */ true);
        prm.declare_entry(
          "direct solver ms",
          "false",
          Patterns::Bool(),
          "Declare if a direct solver shall be used for the coarse scale "
          "of the MsFEM",
          /* has_to_be_set = */ true);
        prm.declare_entry(
          "direct solver basis",
          "false",
          Patterns::Bool(),
          "Declare if a direct solver shall be used for the MsFEM "
          "on the fine scale in each cell.",
          /* has_to_be_set = */ true);

        prm.declare_entry(
          "prevent basis output",
          "true",
          Patterns::Bool(),
          "Declare the output of the basis functions of the MsFEM"
          " should be prevented.",
          /* has_to_be_set = */ true);

        prm.declare_entry("coarse refinements",
                          "1",
                          Patterns::Integer(),
                          "Set the number of refinements on the coarse scale.",
                          /* has_to_be_set = */ true);
        prm.declare_entry("fine refinements",
                          "1",
                          Patterns::Integer(),
                          "Set the number of refinements on the fine scale.",
                          /* has_to_be_set = */ true);

        prm.declare_entry("rotate",
                          "false",
                          Patterns::Bool(),
                          "True if an end of the body shall be rotated.",
                          /* has_to_be_set = */ true);
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
  ElaParameters<dim>::parse_parameters(ParameterHandler &prm)
  {
    bool                        use_E_and_nu;
    std::map<std::string, bool> material_structure;

    prm.enter_subsection("Elasticity Parameters");
    {
      std::string dim_str = std::to_string(dim) + "D";

      prm.enter_subsection(dim_str);
      {
        for (unsigned int i = 0; i < dim; ++i)
          {
            init_p1[i] =
              prm.get_double(std::string("init p1[") + std::to_string(i) + "]");
            init_p2[i] =
              prm.get_double(std::string("init p2[") + std::to_string(i) + "]");
          }
      }
      prm.leave_subsection();

      prm.enter_subsection("Material Parameters");
      {
        // True if E and nu shall be used to declare mu and lambda.
        use_E_and_nu = prm.get_bool("use E and nu");

        int m = 0;
        material_structure.insert(
          std::make_pair("oscillations", prm.get_bool("oscillations")));

        material_structure.insert(
          std::make_pair("x-layers", prm.get_bool("x-layers")));
        if (material_structure["x-layers"] == true)
          ++m;

        material_structure.insert(
          std::make_pair("y-layers", prm.get_bool("y-layers")));
        if (material_structure["y-layers"] == true)
          ++m;

        material_structure.insert(
          std::make_pair("z-layers", prm.get_bool("z-layers")));
        if ((material_structure["z-layers"] == true) && (dim == 3))
          ++m;

        if (material_structure["oscillations"] && m != 0)
          {
            std::cout << "The material can only depend on either "
                      << "oscillations or layers but not both." << std::endl;
            exit(1);
          }

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
            lambda =
              new LamePrmOsc<dim>(lambda_fr, lambda_mean, init_p1, init_p2);

            mu = new LamePrmOsc<dim>(mu_fr, mu_mean, init_p1, init_p2);
          }
        else
          {
            if (dim == 3)
              {
                if (material_structure.at("x-layers") ||
                    material_structure.at("y-layers") ||
                    material_structure.at("z-layers"))
                  {
                    unsigned int n_x_layers = 1;
                    unsigned int n_y_layers = 1;
                    unsigned int n_z_layers = 1;
                    if (material_structure.at("x-layers"))
                      {
                        unsigned int n_x_layers =
                          prm.get_integer("no. of x-layers");
                        if (n_x_layers < 1)
                          {
                            std::cout
                              << "The number of x-layers must be greater than 0 "
                              << "if x-layers are used" << std::endl;
                            exit(-1);
                          }
                      }
                    if (material_structure.at("y-layers"))
                      {
                        unsigned int n_y_layers =
                          prm.get_integer("no. of y-layers");
                        if (n_y_layers < 1)
                          {
                            std::cout
                              << "The number of y-layers must be greater than 0 "
                              << "if y-layers are used" << std::endl;
                            exit(-1);
                          }
                      }
                    if (material_structure.at("z-layers"))
                      {
                        unsigned int n_z_layers =
                          prm.get_integer("no. of z-layers");
                        if (n_z_layers < 1)
                          {
                            std::cout
                              << "The number of z-layers must be greater than 0 "
                              << "if z-layers are used" << std::endl;
                            exit(-1);
                          }
                      }

                    unsigned int n_layers =
                      n_x_layers * n_y_layers * n_z_layers;

                    std::vector<unsigned int> index_set(n_layers);
                    std::iota(index_set.begin(), index_set.end(), 0);
                    std::seed_seq seq{1, 2, 3, 4, 5};
                    std::mt19937  rd(seq);
                    std::shuffle(index_set.begin(), index_set.end(), rd);

                    lambda = new LamePrmLayers<dim>(lambda_mean,
                                                    index_set,
                                                    material_structure,
                                                    init_p1,
                                                    init_p2,
                                                    n_x_layers,
                                                    n_y_layers,
                                                    n_z_layers);

                    mu = new LamePrmLayers<dim>(mu_mean,
                                                index_set,
                                                material_structure,
                                                init_p1,
                                                init_p2,
                                                n_x_layers,
                                                n_y_layers,
                                                n_z_layers);
                  }
                else
                  {
                    lambda = new LamePrmConst<dim>(lambda_mean);
                    mu     = new LamePrmConst<dim>(mu_mean);
                  }
              }
            else
              {
                AssertDimension(dim, 2);
                if (material_structure.at("x-layers") ||
                    material_structure.at("y-layers"))
                  {
                    unsigned int n_x_layers = 1;
                    unsigned int n_y_layers = 1;
                    if (material_structure.at("x-layers"))
                      {
                        unsigned int n_x_layers =
                          prm.get_integer("no. of x-layers");
                        if (n_x_layers < 1)
                          std::cout
                            << "The number of x-layers must be greater than 0 "
                            << "if x-layers are used" << std::endl;
                        exit(-1);
                      }
                    if (material_structure.at("y-layers"))
                      {
                        unsigned int n_y_layers =
                          prm.get_integer("no. of y-layers");
                        if (n_y_layers < 1)
                          std::cout
                            << "The number of y-layers must be greater than 0 "
                            << "if y-layers are used" << std::endl;
                        exit(-1);
                      }

                    unsigned int n_layers = n_x_layers * n_y_layers;

                    std::vector<unsigned int> index_set(n_layers);
                    std::iota(index_set.begin(), index_set.end(), 0);
                    std::seed_seq seq{1, 2, 3, 4, 5};
                    std::mt19937  rd(seq);
                    std::shuffle(index_set.begin(), index_set.end(), rd);

                    lambda = new LamePrmLayers<dim>(lambda_mean,
                                                    index_set,
                                                    material_structure,
                                                    init_p1,
                                                    init_p2,
                                                    n_x_layers,
                                                    n_y_layers);

                    mu = new LamePrmLayers<dim>(mu_mean,
                                                index_set,
                                                material_structure,
                                                init_p1,
                                                init_p2,
                                                n_x_layers,
                                                n_y_layers);
                  }
                else
                  {
                    lambda = new LamePrmConst<dim>(lambda_mean);
                    mu     = new LamePrmConst<dim>(mu_mean);
                  }
              }
          }


        rho = prm.get_double("rho");
      }
      prm.leave_subsection();

      prm.enter_subsection("Other Parameters");
      {
        verbose = prm.get_bool("verbose");

        direct_solver_std   = prm.get_bool("direct solver std");
        direct_solver_ms    = prm.get_bool("direct solver ms");
        direct_solver_basis = prm.get_bool("direct solver basis");

        prevent_basis_output = prm.get_bool("prevent basis output");

        coarse_refinements = prm.get_integer("coarse refinements");
        fine_refinements   = prm.get_integer("fine refinements");

        rotate = prm.get_bool("rotate");
        angle  = prm.get_double("rotation angle") * M_PI;
      }
    }
    prm.leave_subsection();
  }
} // namespace Elasticity

#endif /* _INCLUDE_PARAMETER_FILE_TPP_ */