#ifndef _INCLUDE_ELA_BASIS_H_
#define _INCLUDE_ELA_BASIS_H_

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

#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
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

#include "basis_funs.h"
#include "body_force.h"
#include "mytools.h"
#include "postprocessing.h"
#include "process_parameter_file.h"

// STL
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>


namespace Elasticity
{
  using namespace dealii;

  /**
   * @brief Fine-scale basis class for MsFEMs in linear elasticity.
   *
   * @tparam dim Space dimension
   *
   * This class solves the fine-scale problem of the two-scale
   * MsFEM for linear elasticity problems.
   *
   * It gets a cell from the coarse mesh of the class ElaMs,
   * solves the problem on a finer scale and returns the solution.
   *
   * Furthermore, if given the solution of ElaMs, it can output
   * the solution with the basis functions that are generated here.
   *
   * Based on the implementations in the repository
   * https://github.com/konsim83/MPI-MSFEC/
   */
  template <int dim>
  class ElaBasis
  {
  public:
    /**
     * @brief Construct a new ElaBasis object.
     *
     * @param global_cell The cell id of the cell for which
     *                    this class constructs basis functions.
     * @param first_cell The first cell on this processor.
     * @param local_subdomain The id of the local subdomain, i.e. the id of the
     *                        processor that solves this problem.
     */
    ElaBasis(typename Triangulation<dim>::active_cell_iterator &global_cell,
             const CellId                                      &first_cell_id,
             unsigned int                                       local_subdomain,
             MPI_Comm                  mpi_communicator,
             const ElaParameters<dim> &ela_parameters);

    /**
     * @brief Copy Constructor for Ela Basis.
     *
     * @param other The ElaBasis object that we want to copy.
     */
    ElaBasis(const ElaBasis<dim> &other);

    /**
     * @brief Function that runs the problem.
     *
     * This member function runs all the member functions
     * that are necessary to construct the basis functions
     * for the MsFEM.
     *
     * If this is the first cell on this processor and verbose of
     * #parameters_basis is true, an output for the constructed
     * basis function will be created.
     */
    void
    run();

    /**
     * @brief Returns the #global_element_matrix.
     *
     * @return const FullMatrix<double>&
     */
    const FullMatrix<double> &
    get_global_element_matrix() const;

    /**
     * @brief Returns the #global_element_rhs.
     *
     * @return const Vector<double>&
     */
    const Vector<double> &
    get_global_element_rhs() const;

    /**
     * @brief Set the weights for the local contribution
     *        to the solution of the global MsFEM problem.
     *
     * @param global_weights
     */
    void
    set_global_weights(const std::vector<double> &global_weights);
    void

    /**
     * @brief Creates a .vtu output file with the local contribution
     *        to the global solution with the local basis functions.
     */
    output_global_solution_in_cell();

    /**
     * @brief Return the filename for the local contribution to
     *        global solution.
     *
     * @return const std::string
     */
    const std::string
    get_filename() const;

    /**
     * @brief Get the local contribution to the solution vector
     * on the fine scale.
     */
    const std::vector<Vector<double>>
    get_global_solution();

  private:
    /**
     * @brief Sets up the system.
     *
     * This function sets up the system, i.e. it creates the needed sparse
     * matrices with the correct sparsity pattern and the vectors.
     *
     * Moreover, it creates the constraint vector that contains the Dirichlet
     * boundary conditions (,which are, in this case, the point values of the
     * basis functions on the coarse mesh at the quadrature points of this cell)
     * and the hanging node constraints.
     */
    void
    setup_system();

    /**
     * @brief Assembles the system.
     *
     * This function assembles the #assembled_cell_matrix and the
     * #assembled_cell_rhs for this linear elasticity problem.
     *
     * Based on the
     * <a href="https://www.dealii.org/9.2.0/doxygen/deal.II/step_8.html">
     * step-8</a> tutorial program of deal.ii.
     */
    void
    assemble_system();

    /**
     * @brief Solves the problem at a quadrature point.
     *
     * @param q_point The index of the quadrature point.
     *
     * This function solves the problem at a quadrature point
     * with the corresponding constraints and thus constructs
     * the basis functions for this cell.
     *
     * In #parameters_basis, it can be specified if a direct
     * or iterative (CG-method with SSOR preconditioner) shall
     * be used.
     */
    void
    solve(unsigned int q_point);

    /**
     * @brief Assembles the local contribution to the global system matrix
     *        in ElaMs.
     */
    void
    assemble_global_element_matrix();

    /**
     * @brief Outputs the constructed basis functions of the local cell.
     */
    void
    output_basis();

    MPI_Comm                               mpi_communicator;
    Triangulation<dim>                     triangulation;
    FESystem<dim>                          fe;
    DoFHandler<dim>                        dof_handler;
    std::vector<AffineConstraints<double>> constraints_vector;
    std::vector<Point<dim>>                corner_points;
    std::vector<Vector<double>>            solution_vector;
    SparsityPattern                        sparsity_pattern;
    Vector<double>                         assembled_cell_rhs;
    SparseMatrix<double>                   assembled_cell_matrix;
    Vector<double>                         global_element_rhs;
    FullMatrix<double>                     global_element_matrix;
    std::vector<double>                    global_weights;
    Vector<double>                         system_rhs;
    SparseMatrix<double>                   system_matrix;
    Vector<double>                         global_solution;
    const CellId                           global_cell_id;
    const CellId                           first_cell_id;
    const unsigned int                     local_subdomain;
    const ElaParameters<dim>               ela_parameters;
    std::string                            filename;
    BasisFun::BasisQ1<dim>                 basis_q1;
  };
} // namespace Elasticity

#endif // _INCLUDE_ELA_BASIS_H_