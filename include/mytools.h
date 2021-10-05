#ifndef _INCLUDE_MY_TOOLS_H_
#define _INCLUDE_MY_TOOLS_H_

#include <deal.II/base/point.h>

#include <deal.II/distributed/tria.h>

#include <sys/stat.h>

#include <stdexcept>
#include <string>

#include "process_parameter_file.h"

/**
 * @namespace MyTools
 *
 * @brief Different tools that could also be used
 *        by functions unrelated to elasticity.
 */
namespace MyTools
{
  using namespace dealii;

  /****************************************************************************/
  /* Tools */

  //
  /**
   * @brief Sets a new boundary id for subset of the surface.
   *
   * @tparam dim Space dimension
   * @param p1 First point
   * @param p2 Second point
   * @param face_id Old face id
   * @param id New face id
   * @param triangulation Triangulation
   *
   * This is a function that gives a specific subset of the surface
   * located on the face with the id face_id a new face id.
   */
  template <int dim>
  void
  set_dirichlet_id(const Point<dim> &                         p1,
                   const Point<dim> &                         p2,
                   const unsigned int                         face_id,
                   const unsigned int                         id,
                   parallel::distributed::Triangulation<dim> &triangulation);


  /**
   * @brief Get the repetitions for a subdivided_hyper_rectangle
   *
   * @tparam dim Space dimension
   * @param p1 First vertex point
   * @param p2 Second vertex point
   * @return const std::vector<unsigned int>
   *
   * This function computes the number of repeated refinements in x- and y-
   * (and if dim=3 in z-)direction such that the cells resemble squares/cubes
   * as close as possible.
   */
  template <int dim>
  const std::vector<unsigned int>
  get_repetitions(const Point<dim> &p1, const Point<dim> &p2);


  /*!
   * @brief Creates a directory with a name from a given string
   *
   * @param dir_name Path to a directory
   */
  void
  create_data_directory(const char *dir_name);
} // namespace MyTools

#endif // _INCLUDE_MY_TOOLS_H_