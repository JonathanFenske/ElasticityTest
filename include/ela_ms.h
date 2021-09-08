// Deal.II
#ifndef _INCLUDE_ELA_MS_H_
#define _INCLUDE_ELA_MS_H_

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/physics/transformations.h>

#include "forces_and_parameters.h"
#include "postprocessing.h"
#include "mytools.h"
#include "ela_basis.h"

// STL
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <filesystem>


namespace Elasticity
{
  using namespace dealii;

  /****************************************************************************/
 /* Class for the coarse scale part of the multiscale implementation for 
    linear elasticity problems */

  // class that enables solving linear elasticity problems in parallel
  // which is based on the step-8 and step-40 tutorials of deal.ii
  template <int dim>
  class ElaMs
  {
  public:
    ElaMs(const bool direct_solver, const bool neumann_bc);
    void
    run();

  private:
    const std::tuple<Point<dim>, Point<dim>>
    get_init_vert(std::vector<double> p) const;
    void
    setup_system();
    void
    initialize_and_compute_basis();
    void
    assemble_system();
    void
    solve();
    void
    send_global_weights_to_cell();
    void
    refine_grid();
    void
    output_results();

    MPI_Comm                                  mpi_communicator;
    parallel::distributed::Triangulation<dim> triangulation;
    FESystem<dim>                             fe;
    DoFHandler<dim>                           dof_handler;
    IndexSet                                  locally_owned_dofs;
    IndexSet                                  locally_relevant_dofs;
    AffineConstraints<double>                 constraints;
    TrilinosWrappers::SparseMatrix            system_matrix;
    TrilinosWrappers::SparseMatrix            preconditioner_matrix;
    TrilinosWrappers::MPI::Vector             locally_relevant_solution;
    TrilinosWrappers::MPI::Vector             system_rhs;
    std::map<CellId, ElaBasis<dim>>           cell_basis_map;
    ConditionalOStream                        pcout;
    TimerOutput                               computing_timer;
    const bool                                direct_solver;
    const bool                                neumann_bc;
  };


  // template <int dim>
  // class ElaMs : public ElaMs<dim>
  // {
  // };
} // namespace Elasticity

#endif // _INCLUDE_ELA_MS_H_