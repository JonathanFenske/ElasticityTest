#include "process_parameter_file.h"

#include "myexceptions.h"
#include "process_parameter_file.tpp"

namespace Elasticity
{
  template struct ElaParameters<2>;
  template struct ElaParameters<3>;

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
      AssertThrow((dim == 2) || (dim == 3), MyExcDims(dim, 2, 3));
    }
    prm.leave_subsection();
  }
} // namespace Elasticity