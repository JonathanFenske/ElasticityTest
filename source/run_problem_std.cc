#include "run_problem_std.h"


namespace Elasticity
{
  void
  run_2d_problem_std(const std::string &input_file)
  {
    ElaParameters<2> ela_parameters(input_file);

    ElaStd<2> ela_std(ela_parameters);
    ela_std.run();
  }

  void
  run_3d_problem_std(const std::string &input_file)
  {
    ElaParameters<3> ela_parameters(input_file);

    ElaStd<3> ela_std(ela_parameters);
    ela_std.run();
  }
} // namespace Elasticity