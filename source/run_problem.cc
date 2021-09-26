#include "run_problem.h"


namespace Elasticity
{
  void
  run_2d_problem(const std::string &input_file)
  {
    GlobalParameters<2> global_parameters(input_file);

    {
      ParametersStd parameters_std(input_file);
      ElaStd<2>     ela_std(global_parameters, parameters_std);
      ela_std.run();
    }

    {
      ParametersMs    parameters_ms(input_file);
      ParametersBasis parameters_basis(input_file);
      ElaMs<2> ela_ms(global_parameters, parameters_ms, parameters_basis);
      ela_ms.run();
    }
  }

  void
  run_3d_problem(const std::string &input_file)
  {
    GlobalParameters<3> global_parameters(input_file);

    {
      ParametersStd parameters_std(input_file);
      ElaStd<3>     ela_std(global_parameters, parameters_std);
      ela_std.run();
    }

    {
      ParametersMs    parameters_ms(input_file);
      ParametersBasis parameters_basis(input_file);
      ElaMs<3> ela_ms(global_parameters, parameters_ms, parameters_basis);
      ela_ms.run();
    }
  }
} // namespace Elasticity