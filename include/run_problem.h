#ifndef _RUN_PROBLEM_H_
#define _RUN_PROBLEM_H_

#include "ela_ms.h"
#include "ela_std.h"


namespace Elasticity
{
  void
  run_2d_problem(const std::string &input_file);

  void
  run_3d_problem(const std::string &input_file);
} // namespace Elasticity

#endif // _RUN_PROBLEM_H_