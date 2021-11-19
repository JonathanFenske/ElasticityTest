#include "run_problem_ms.h"


namespace Elasticity
{
  void
  run_2d_problem_ms(const std::string &input_file)
  {
    ElaParameters<2> ela_parameters(input_file);

    ElaMs<2> ela_ms(ela_parameters);
    ela_ms.run();
  }

  void
  run_3d_problem_ms(const std::string &input_file)
  {
    ElaParameters<3> ela_parameters(input_file);

    ElaMs<3> ela_ms(ela_parameters);
    ela_ms.run();
  }
} // namespace Elasticity