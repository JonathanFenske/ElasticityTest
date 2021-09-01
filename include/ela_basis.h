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

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

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
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
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
#include "basis_funs.h"

// STL
#include <cmath>
#include <fstream>
#include <iostream>


namespace Elasticity
{
  using namespace dealii;

  /****************************************************************************/
  /* The main class*/

  // class that enables solving linear elasticity problems in parallel
  // which is based on the step-8 and step-40 tutorials of deal.ii
  template <int dim>
  class ElaBasis
  {
  public:
    ElaBasis(typename Triangulation<dim>::active_cell_iterator &global_cell,
              unsigned int local_subdomain,
              MPI_Comm     mpi_communicator,
              const bool   direct_solver);
    ElaBasis(const ElaBasis<dim> &other);

    void
    run();
    const FullMatrix<double> &
    get_global_element_matrix() const;
    const Vector<double> &
    get_global_element_rhs() const;

  private:
    void
    setup_system();
    void
    assemble_system();
    void
    solve(unsigned int q_point);
    void
    assemble_global_element_matrix();

    MPI_Comm                                  mpi_communicator;
    Triangulation<dim>                        triangulation;
    FESystem<dim>                             fe;
    DoFHandler<dim>                           dof_handler;
    std::vector<AffineConstraints<double>>    constraints_vector;
    std::vector<Point<dim>>                   corner_points;
    std::vector<Vector<double>>               solution_vector;
    SparsityPattern                           sparsity_pattern;
    Vector<double>                            assembled_cell_rhs;
    SparseMatrix<double>                      assembled_cell_matrix;
    Vector<double>                            global_element_rhs;
    FullMatrix<double>                        global_element_matrix;
    Vector<double>                            system_rhs;
    SparseMatrix<double>                      system_matrix;
    Vector<double>                            global_solution;
    const CellId                              global_cell_id;
    const unsigned int                        local_subdomain;
    BasisFun::BasisQ1Grad<dim>                basis_q1_grad;
    ConditionalOStream                        pcout;
    const bool                                direct_solver;        
  };


  // template <int dim>
  // class ElaBasis : public ElaBasis<dim>
  // {
  // };
} // namespace Elasticity

#endif // _INCLUDE_ELA_BASIS_H_