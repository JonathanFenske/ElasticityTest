#ifndef _INCLUDE_FINE_MS_SOLUTION_FUNCTION_H_
#define _INCLUDE_FINE_MS_SOLUTION_FUNCTION_H_

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/numerics/vector_tools.h>

#include <map>

#include "ela_basis.h"


namespace Elasticity
{
  using namespace dealii;

  /**
   * @brief Function based on the fine scale MsFEM solution
   *
   * @tparam dim Space dimension
   */
  template <int dim>
  class FineMsSolution : public Function<dim>
  {
  public:
    // FineMsSolution()
    //   : Function<dim>(dim)
    //   , triangulation(typename Triangulation<dim>::MeshSmoothing(
    //       Triangulation<dim>::smoothing_on_refinement |
    //       Triangulation<dim>::smoothing_on_coarsening))
    // {}

    FineMsSolution(const Triangulation<dim>                &other_triangulation,
                   const std::map<CellId, Function<dim> *> &function_map,
                   MPI_Comm                                 mpi_communicator)
      : Function<dim>(dim)
      , mpi_communicator(mpi_communicator)
      , triangulation(mpi_communicator,
                      typename Triangulation<dim>::MeshSmoothing(
                        Triangulation<dim>::smoothing_on_refinement |
                        Triangulation<dim>::smoothing_on_coarsening))
      , function_map(function_map)
    {
      triangulation.copy_triangulation(other_triangulation);
    }

    FineMsSolution(const FineMsSolution<dim> &other)
      : Function<dim>(dim)
      , mpi_communicator(other.mpi_communicator)
      , triangulation(mpi_communicator,
                      typename Triangulation<dim>::MeshSmoothing(
                        Triangulation<dim>::smoothing_on_refinement |
                        Triangulation<dim>::smoothing_on_coarsening))
      , function_map(other.function_map)
    {
      triangulation.copy_triangulation(other.triangulation);
    }

    // FineMsSolution<dim>
    // operator=(const FineMsSolution<dim> &other);

    virtual void
    vector_value(const Point<dim> &p,
                 Vector<double>   &vector_value) const override;

    virtual void
    vector_gradient(const Point<dim>            &p,
                    std::vector<Tensor<1, dim>> &gradients) const override;

  private:
    MPI_Comm                             mpi_communicator;
    parallel::shared::Triangulation<dim> triangulation;
    std::map<CellId, Function<dim> *>    function_map;
  };

  // exernal template instantiations
  extern template class FineMsSolution<2>;
  extern template class FineMsSolution<3>;
} // namespace Elasticity

#endif // _INCLUDE_FINE_MS_SOLUTION_FUNCTION_H_