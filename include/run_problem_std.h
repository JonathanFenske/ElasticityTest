#ifndef _RUN_PROBLEM_STD_H_
#define _RUN_PROBLEM_STD_H_

#include "ela_std.h"

/**
 * @namespace Elasticity
 *
 * @brief Contains all functions and classes
 *        related to linear elasticity problems.
 */
namespace Elasticity
{
  /**
   * @brief Run the 2D problem with the Standard FEM.
   *
   * @param input_file Path to parameter file
   */
  void
  run_2d_problem_std(const std::string &input_file);

  /**
   * @brief Run the 3D problem with the Standard FEM.
   *
   * @param input_file Path to parameter file
   */
  void
  run_3d_problem_std(const std::string &input_file);
} // namespace Elasticity

#endif // _RUN_PROBLEM_STD_H_