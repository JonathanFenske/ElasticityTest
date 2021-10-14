#include "process_parameter_file.h"

#include "process_parameter_file.tpp"

namespace Elasticity
{
  template struct GlobalParameters<2>;
  template struct GlobalParameters<3>;

  using namespace dealii;

  Dimension::Dimension(const std::string &parameter_filename)
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


  void
  Dimension::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Dimension");
    {
      prm.declare_entry("dim", "3", Patterns::Integer(), "Set the Dimension.");
    }
    prm.leave_subsection();
  }


  void
  Dimension::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Dimension");
    {
      dim = prm.get_integer("dim");
    }
    prm.leave_subsection();
  }


  ParametersStd::ParametersStd(const std::string &parameter_filename)
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


  void
  ParametersStd::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Standard method parameters");
    {
      prm.enter_subsection("Bools");
      {
        prm.declare_entry("verbose",
                          "true",
                          Patterns::Bool(),
                          "Choose whether to verbose.");
        prm.declare_entry("use direct solver",
                          "true",
                          Patterns::Bool(),
                          "Choose whether to use a direct solver.");
      }
      prm.leave_subsection();

      prm.enter_subsection("Mesh");
      {
        prm.declare_entry("refinements",
                          "3",
                          Patterns::Integer(0, 10),
                          "Number of initial mesh refinements.");
        prm.declare_entry("cycles",
                          "3",
                          Patterns::Integer(1, 10),
                          "Number of cycles that the problems runs through.");
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }


  void
  ParametersStd::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Standard method parameters");
    {
      prm.enter_subsection("Bools");
      {
        verbose       = prm.get_bool("verbose");
        direct_solver = prm.get_bool("use direct solver");
      }
      prm.leave_subsection();

      prm.enter_subsection("Mesh");
      {
        n_refine = prm.get_integer("refinements");
        n_cycles = prm.get_integer("cycles");
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }


  ParametersMs::ParametersMs(const std::string &parameter_filename)
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


  void
  ParametersMs::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Multiscale method parameters");
    {
      prm.enter_subsection("Coarse scale");
      {
        prm.enter_subsection("Bools");
        {
          prm.declare_entry("verbose",
                            "true",
                            Patterns::Bool(),
                            "Choose whether to verbose.");
          prm.declare_entry("use direct solver",
                            "true",
                            Patterns::Bool(),
                            "Choose whether to use a direct solver.");
        }
        prm.leave_subsection();

        prm.enter_subsection("Mesh");
        {
          prm.declare_entry("refinements",
                            "3",
                            Patterns::Integer(0, 10),
                            "Number of initial mesh refinements.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }

  void
  ParametersMs::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Multiscale method parameters");
    {
      prm.enter_subsection("Coarse scale");
      {
        prm.enter_subsection("Bools");
        {
          verbose       = prm.get_bool("verbose");
          direct_solver = prm.get_bool("use direct solver");
        }
        prm.leave_subsection();

        prm.enter_subsection("Mesh");
        {
          n_refine = prm.get_integer("refinements");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }


  ParametersBasis::ParametersBasis(const std::string &parameter_filename)
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


  void
  ParametersBasis::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Multiscale method parameters");
    {
      prm.enter_subsection("Fine scale");
      {
        prm.enter_subsection("Bools");
        {
          prm.declare_entry("verbose",
                            "true",
                            Patterns::Bool(),
                            "Choose whether to verbose.");
          prm.declare_entry("use direct solver",
                            "true",
                            Patterns::Bool(),
                            "Choose whether to use a direct solver.");
          prm.declare_entry(
            "prevent output",
            "true",
            Patterns::Bool(),
            "Choose whether to prevent the output on the fine scale.");
        }
        prm.leave_subsection();

        prm.enter_subsection("Mesh");
        {
          prm.declare_entry("refinements",
                            "3",
                            Patterns::Integer(1, 10),
                            "Number of initial mesh refinements.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }


  void
  ParametersBasis::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Multiscale method parameters");
    {
      prm.enter_subsection("Fine scale");
      {
        prm.enter_subsection("Bools");
        {
          verbose        = prm.get_bool("verbose");
          direct_solver  = prm.get_bool("use direct solver");
          prevent_output = prm.get_bool("prevent output");
        }
        prm.leave_subsection();

        prm.enter_subsection("Mesh");
        {
          n_refine = prm.get_integer("refinements");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }
} // namespace Elasticity