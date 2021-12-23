#include "run_problem_compare.h"


namespace Elasticity
{
  void
  run_2d_problem_compare(const std::string &input_file)
  {
    ElaParameters<2> ela_parameters(input_file, true);

    dealii::Vector<double> ms_solution;

    {
      ElaMs<2> ela_ms(ela_parameters);
      ela_ms.run();

      ms_solution = ela_ms.get_fine_solution();
    }

    {
      ElaStd<2> ela_std(ela_parameters);
      ela_std.run();

      ela_std.compare_solutions(ms_solution);
    }
  }

  void
  run_3d_problem_compare(const std::string &input_file)
  {
    ElaParameters<3> ela_parameters(input_file, true);

    dealii::Vector<double> ms_solution;

    {
      ElaMs<3> ela_ms(ela_parameters);
      ela_ms.run();

      ms_solution = ela_ms.get_fine_solution();
    }

    {
      ElaStd<3> ela_std(ela_parameters);
      ela_std.run();

      ela_std.compare_solutions(ms_solution);
    }
  }
} // namespace Elasticity