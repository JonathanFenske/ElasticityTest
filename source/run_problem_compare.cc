#include "run_problem_compare.h"


namespace Elasticity
{
  void
  run_2d_problem_compare(const std::string &input_file)
  {
    ElaParameters<2> ela_parameters(input_file, true);

    Vector<double> coarse_solution;
    Vector<double> fine_solution;

    {
      ElaStd<2> ela_std(ela_parameters);
      ela_std.run();

      ela_std.get_solutions(coarse_solution, fine_solution);
    }

    {
      ElaMs<2> ela_ms(ela_parameters);
      ela_ms.run();

      ela_ms.compute_errors(coarse_solution, fine_solution);
    }
  }

  void
  run_3d_problem_compare(const std::string &input_file)
  {
    ElaParameters<3> ela_parameters(input_file, true);

    Vector<double> coarse_solution;
    Vector<double> fine_solution;

    {
      ElaStd<3> ela_std(ela_parameters);
      ela_std.run();

      ela_std.get_solutions(coarse_solution, fine_solution);
    }

    {
      ElaMs<3> ela_ms(ela_parameters);
      ela_ms.run();

      ela_ms.compute_errors(coarse_solution, fine_solution);
    }
  }
} // namespace Elasticity