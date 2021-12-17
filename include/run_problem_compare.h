#ifndef _RUN_PROBLEM_COMPARE_H_
#define _RUN_PROBLEM_COMPARE_H_

#include "ela_ms.h"
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
   * @brief Run the 2D problem with the MsFEM.
   *
   * @param input_file Path to parameter file
   */
  void
  run_2d_problem_compare(const std::string &input_file);

  /**
   * @brief Run the 3D problem with the MsFEM.
   *
   * @param input_file Path to parameter file
   */
  void
  run_3d_problem_compare(const std::string &input_file);
} // namespace Elasticity

#endif // _RUN_PROBLEM_COMPARE_H_