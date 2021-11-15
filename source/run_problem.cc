#include "run_problem.h"


namespace Elasticity
{
  void
  run_2d_problem(const std::string &input_file)
  {
    ElaParameters<2> ela_parameters(input_file);

    {
      ElaStd<2> ela_std(ela_parameters);
      ela_std.run();
    }

    {
      ElaMs<2> ela_ms(ela_parameters);
      ela_ms.run();
    }
  }

  void
  run_3d_problem(const std::string &input_file)
  {
    ElaParameters<3> ela_parameters(input_file);

    {
      ElaStd<3> ela_std(ela_parameters);
      ela_std.run();
    }

    {
      ElaMs<3> ela_ms(ela_parameters);
      ela_ms.run();
    }
  }
} // namespace Elasticity