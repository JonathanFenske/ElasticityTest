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

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
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

#include "ela_basis.h"
// #include "fine_ms_solution_function.h"
#include "mytools.h"
#include "postprocessing.h"
#include "process_parameter_file.h"

// include headers that implement a archive in simple text format
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/array.hpp>

// STL
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <utility>


namespace Elasticity
{
  using namespace dealii;

  /**
   * @brief Class for the coarse scale of the MsFEM for linear elasticity.
   *
   * @tparam dim Space dimension
   *
   * This class computes and outputs the global solution of
   * the given linear elasticity problem.
   *
   * Based on the implementations in the repository
   * https://github.com/konsim83/MPI-MSFEC/ and the step-40
   * tutorial of deal.ii.
   */
  template <int dim>
  class ElaMs
  {
  public:
    /**
     * @brief Construct a new ElaMs object.
     */
    ElaMs(const ElaParameters<dim> &ela_parameters);

    /**
     * @brief Function that runs the problem.
     *
     * This function runs all the member functions
     * that are necessary to solve the problem.
     */
    void
    run();

    /**
     * @brief This function computes the errors of the MsFEM and the
     * standard FEM solution. The fine scale standard FEM solution is
     * used as reference.
     *
     * @param coarse_solution coarse scale standard FEM solution
     * @param fine_solution fine scale standard FEM solution
     */
    void
    compute_errors(Vector<double> &coarse_solution,
                   Vector<double> &fine_solution);

  private:
    /**
     * @brief Sets up the system.
     *
     * This function sets up the system, i.e. it creates the needed sparse
     * matrices with the correct sparsity pattern and the vectors.
     *
     * Moreover, it creates the constraint vector that contains
     * the global Dirichlet boundary conditions
     * and the hanging node constraints.
     */
    void
    setup_system();

    /**
     * @brief Initializes and computes the basis functions.
     *
     * This function initializes and computes the local basis functions
     * on each cell with the MsFEM by creating an ElaBasis object for
     * each cell and computing the basis functions with these objects.
     */
    void
    initialize_and_compute_basis();

    /**
     * @brief Assembles the system.
     *
     * This function assembles the #system_matrix and the #system_rhs by
     * assembling the contributions of all the ElaBasis objects and adding
     * the Neumann boundary condition.
     */
    void
    assemble_system();

    /**
     * @brief Solves the global problem.
     *
     * This function solves the global problem.
     *
     * In #parameters_ms, it can be specified if a direct
     * or iterative (CG-method with AMG preconditioner) shall
     * be used.
     */
    void
    solve();

    /**
     * @brief Sends global weights to cell.
     *
     * For each cell, this function sends the respective local part of the
     * solution vector to the corresponding ElaBasis object such that the
     * computed solution vector can be applied to the local basis functions
     * of the MsFEM.
     */
    void
    send_global_weights_to_cell();

    /**
     * @brief Adaptively refines grid.
     */
    void
    refine_grid();

    /**
     * @brief Outputs the solutions in pvtu files.
     *
     * This function outputs the global solution with the coarse basis functions
     * as well as with the constructed multiscale basis functions as pvtu files.
     *
     * In the first case, it creates vtu files for every subdomain
     * (corresponding to the respective processor) and combines them into
     * a single pvtu file.
     *
     * In the latter case, it lets each ElaBasis object output the
     * local solution as vtu files and combines all of them into a
     * single pvtu file.
     */
    void
    output_results();

    /**
     * @brief Outputs vtu files and a pvtu file for fine_solution which must
     * live on dof_handler_fine.
     */
    void
    output_fine_solution(
      parallel::shared::Triangulation<dim> &triangulation_fine,
      DoFHandler<dim>                      &dof_handler_fine,
      const Vector<double>                 &fine_solution,
      std::string                           name);

    /**
     * @brief Assembles the solution vector of the MsFEM on the fine scale.
     *
     */
    const Vector<double>
    get_fine_solution(DoFHandler<dim> &dof_handler_fine);

    MPI_Comm                             mpi_communicator;
    parallel::shared::Triangulation<dim> triangulation;
    FESystem<dim>                        fe;
    DoFHandler<dim>                      dof_handler;
    IndexSet                             locally_owned_dofs;
    IndexSet                             locally_relevant_dofs;
    AffineConstraints<double>            constraints;
    TrilinosWrappers::SparseMatrix       system_matrix;
    TrilinosWrappers::SparseMatrix       preconditioner_matrix;
    TrilinosWrappers::MPI::Vector        locally_relevant_solution;
    TrilinosWrappers::MPI::Vector        system_rhs;
    CellId                               first_cell_id;
    std::map<CellId, ElaBasis<dim>>      cell_basis_map;
    const ElaParameters<dim>             ela_parameters;
    bool                                 processor_is_used;
    /**< True if this processor is assigned at least one coarse cell. */

    ConditionalOStream pcout;
    TimerOutput        computing_timer;
  };
} // namespace Elasticity

#endif // _INCLUDE_ELA_MS_H_