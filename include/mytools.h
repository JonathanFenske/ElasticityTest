#ifndef _INCLUDE_MY_TOOLS_H_
#define _INCLUDE_MY_TOOLS_H_

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/physics/transformations.h>

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
  set_dirichlet_id(const Point<dim>                          &p1,
                   const Point<dim>                          &p2,
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


  template <int dim>
  class Rotation : public Function<dim>
  {
  public:
    Rotation(Point<dim> init_p1, Point<dim> init_p2);
    virtual void
    vector_value(const Point<dim> &p, Vector<double> &values) const override;

  private:
    double       y;
    double       z;
    Tensor<2, 3> R;
  };


  /*
   * Generate a random unsigned integer number between [a,b)
   */
  // class RandomNumberUInt
  // {
  // public:
  //   RandomNumberUInt(const unsigned int b, const bool same_on_all_ranks =
  //   true);

  //   RandomNumberUInt(const unsigned int a,
  //                    const unsigned int b,
  //                    const bool         same_on_all_ranks = true);

  //   void
  //   reinit();

  //   double
  //   generate();

  // private:
  //   unsigned int a, b;

  //   bool same_on_all_ranks;

  //   uint64_t timeSeed;

  //   std::seed_seq seed_sequence;

  //   std::mt19937_64 rng;
  // };
} // namespace MyTools

#endif // _INCLUDE_MY_TOOLS_H_