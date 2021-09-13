#ifndef _INCLUDE_MY_TOOLS_H_
#define _INCLUDE_MY_TOOLS_H_

#include <deal.II/base/point.h>

#include <deal.II/distributed/tria.h>

#include <sys/stat.h>

#include <stdexcept>
#include <string>

#include "process_parameter_file.h"

namespace MyTools
{
  using namespace dealii;

  /****************************************************************************/
  /* Tools */

  // Function that gives a specific subset of the surface a new face_id
  template <int dim>
  void
  set_dirichlet_id(const Point<dim> &                         p1,
                   const Point<dim> &                         p2,
                   const unsigned int                         face_id,
                   const unsigned int                         id,
                   parallel::distributed::Triangulation<dim> &triangulation);


  template <int dim>
  const std::vector<unsigned int>
  get_repetitions(const Point<dim> &p1, const Point<dim> &p2);


  /*!
   * @brief Creates a directory with a name from a given string
   *
   * @param dir_name
   */
  void
  create_data_directory(const char *dir_name);


  auto
  get_global_parameters(int dim, const std::string &parameter_filename);
} // namespace MyTools

#endif // _INCLUDE_MY_TOOLS_H_