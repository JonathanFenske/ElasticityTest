#ifndef _INCLUDE_ELA_BASIS_TPP_
#define _INCLUDE_ELA_BASIS_TPP_

#include "ela_basis.h"


namespace Elasticity
{
  using namespace dealii;

  /****************************************************************************/
  /* Class for the fine scale part of the multiscale implementation for 
     linear elasticity problems */

  // The constructor
  template <int dim>
  ElaBasis<dim>::ElaBasis(typename Triangulation<dim>::active_cell_iterator 
                          &global_cell,
                          unsigned int local_subdomain,
                          MPI_Comm     mpi_communicator,
                          const bool   direct_solver)
    : mpi_communicator(mpi_communicator)
    , triangulation()
    , fe(FE_Q<dim>(1), dim)
    , dof_handler(triangulation)
    , constraints_vector(GeometryInfo<dim>::vertices_per_cell)
    , corner_points(GeometryInfo<dim>::vertices_per_cell)
    , solution_vector(GeometryInfo<dim>::vertices_per_cell)
    , global_element_rhs(fe.dofs_per_cell)
    , global_element_matrix(fe.dofs_per_cell, fe.dofs_per_cell)
    , global_cell_id(global_cell->id())
    , local_subdomain(local_subdomain)
    , basis_q1_grad(global_cell)
    , pcout(std::cout,
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
    // If true, a direct solver will be used to solve the problem.
    , direct_solver(direct_solver)
  {
    // set corner points
    for (unsigned int vertex_n = 0;
         vertex_n < GeometryInfo<dim>::vertices_per_cell;
         ++vertex_n)
      {
        corner_points[vertex_n] = global_cell->vertex(vertex_n);
      }
  }


  template <int dim>
  ElaBasis<dim>::ElaBasis(const ElaBasis<dim> &other)
    : mpi_communicator(other.mpi_communicator)
    , triangulation()
    , fe(FE_Q<dim>(1), dim)
    , dof_handler(triangulation)
    , constraints_vector(other.constraints_vector)
    , corner_points(other.corner_points)
    , solution_vector(other.solution_vector)
    , global_element_rhs(other.global_element_rhs)
    , global_element_matrix(other.global_element_matrix)
    , global_cell_id(other.global_cell_id)
    , local_subdomain(other.local_subdomain)
    , basis_q1_grad(other.basis_q1_grad)
    , pcout(std::cout,
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
    // If true, a direct solver will be used to solve the problem.
    , direct_solver(other.direct_solver)
  {}


  template <int dim>
  void
  ElaBasis<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);
    DynamicSparsityPattern     dsp(dof_handler.n_dofs());

    for (unsigned int q_point = 0;
         q_point < GeometryInfo<dim>::vertices_per_cell;
         ++q_point)
      {
        basis_q1_grad.set_q_point(q_point);

        constraints_vector[q_point].clear();
        DoFTools::make_hanging_node_constraints(
          dof_handler, constraints_vector[q_point]);

        VectorTools::interpolate_boundary_values(
          dof_handler,
          /*boundary id*/ 0,
          basis_q1_grad,
          constraints_vector[q_point]);
        constraints_vector[q_point].close();
      }

    DoFTools::make_sparsity_pattern(
      dof_handler,
      dsp,
      constraints_vector[0], // sparsity pattern is the same for each basis
      /*keep_constrained_dofs =*/true); // for time stepping this is essential
                                        // to be true
    sparsity_pattern.copy_from(dsp);

    assembled_cell_matrix.reinit(sparsity_pattern);
    system_matrix.reinit(sparsity_pattern);

    for (unsigned int q_point = 0;
         q_point < GeometryInfo<dim>::vertices_per_cell;
         ++q_point)
      {
        solution_vector[q_point].reinit(dof_handler.n_dofs());
      }
    assembled_cell_rhs.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
  }


  template <int dim>
  void
  ElaBasis<dim>::assemble_system()
  {
    const QGauss<dim>     quadrature_formula(fe.degree + 1);
    FEValues<dim>         fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);
    const unsigned int    dofs_per_cell   = fe.n_dofs_per_cell();
    const unsigned int    n_q_points      = quadrature_formula.size();
    
    lambda<dim>           lambda;
    mu<dim>               mu;
    BodyForce<dim>        body_force;
    std::vector<Vector<double>>   body_force_values(n_q_points);
    for (unsigned int i = 0; i < n_q_points; ++i)
      body_force_values[i].reinit(dim);
    std::vector<double>   lambda_values(n_q_points), mu_values(n_q_points);

    FullMatrix<double>    local_cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>        local_cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        local_cell_matrix = 0.;
        local_cell_rhs    = 0.;
        fe_values.reinit(cell);
        lambda.value_list(fe_values.get_quadrature_points(), lambda_values);
        mu.value_list(fe_values.get_quadrature_points(), mu_values);
        body_force.vector_value_list(fe_values.get_quadrature_points()
                                      , body_force_values);
        for (unsigned int q_index = 0; q_index < n_q_points; ++q_index)
          {
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                const unsigned int component_i =
                  fe.system_to_component_index(i).first;
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  {
                    const unsigned int component_j =
                      fe.system_to_component_index(j).first;
                    local_cell_matrix(i, j) +=
                      ((fe_values.shape_grad(i, q_index)[component_i] *
                        fe_values.shape_grad(j, q_index)[component_j] *
                        lambda_values[q_index]) +
                        (fe_values.shape_grad(i, q_index)[component_j] *
                        fe_values.shape_grad(j, q_index)[component_i] *
                        mu_values[q_index]) +
                        ((component_i == component_j) ?
                          (fe_values.shape_grad(i, q_index) *
                            fe_values.shape_grad(j, q_index) *
                            mu_values[q_index]) :
                          0)) *
                      fe_values.JxW(q_index);
                  }
                local_cell_rhs(i) +=
                  fe_values.shape_value_component(i, q_index, component_i) *
                  body_force_values[q_index][component_i] *
                  fe_values.JxW(q_index);
              }
          }
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
          {
            assembled_cell_matrix.add(local_dof_indices[i],
                              local_dof_indices[j],
                              local_cell_matrix(i,j));
          }
          assembled_cell_rhs(local_dof_indices[i]) += local_cell_rhs(i);
        }
              
      }
  }


  template <int dim>
  void
  ElaBasis<dim>::solve(unsigned int q_point)
  {
    if (direct_solver)
      {
        SparseDirectUMFPACK A_inv;
        A_inv.initialize(system_matrix);
        A_inv.vmult(solution_vector[q_point], system_rhs);

        constraints_vector[q_point].distribute(solution_vector[q_point]);
      }
    else
      {
        unsigned int  n_iterations     = dof_handler.n_dofs();
        const double  solver_tolerance = 1e-8 * assembled_cell_rhs.l2_norm();
        SolverControl solver_control(
          /* n_max_iter */ n_iterations,
          solver_tolerance,
          /* log_history */ true,
          /* log_result */ true);

        SolverCG<> solver(solver_control);

        PreconditionSSOR<> preconditioner;
        preconditioner.initialize(system_matrix, 1.6);

        try
          {
            solver.solve(system_matrix,
                 solution_vector[q_point],
                 system_rhs,
                 preconditioner);
          }
        catch (std::exception &e)
          {
            Assert(false, ExcMessage(e.what()));
          }

        std::cout << "   "
                << "(cell   " << global_cell_id.to_string() << ") "
                << "(basis   " << q_point << ")   "
                << solver_control.last_step()
                << " fine CG iterations needed to obtain convergence."
                << std::endl;
      }
  }


  template<int dim>
  void
  ElaBasis<dim>::assemble_global_element_matrix()
  {
    // First, reset.
    global_element_matrix = 0;

    // Get lengths of tmp vectors for assembly
    const unsigned int dofs_per_cell = fe.dofs_per_cell;

    Vector<double> tmp(dof_handler.n_dofs());

    // This assembles the local contribution to the global global matrix
    // with an algebraic trick. It uses the local system matrix stored in
    // the respective basis object.
    for (unsigned int i_test = 0; i_test < dofs_per_cell; ++i_test)
      {
        // set an alias name
        const Vector<double> &test_vec = solution_vector[i_test];

        for (unsigned int i_trial = 0; i_trial < dofs_per_cell; ++i_trial)
          {
            // set an alias name
            const Vector<double> &trial_vec = solution_vector[i_trial];

            // tmp = system_matrix*trial_vec

            assembled_cell_matrix.vmult(tmp, trial_vec);

            // global_element_matrix = test_vec*tmp
            global_element_matrix(i_test, i_trial) += (test_vec * tmp);

            // reset
            tmp = 0;
          } // end for i_trial

        global_element_rhs(i_test) += test_vec * assembled_cell_rhs;

      } // end for i_test
  }


  template <int dim>
  void
  ElaBasis<dim>::run()
  {
    Timer timer;

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int  name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    std::string proc_name(processor_name, name_len);

    std::cout << "	Solving for basis in cell   "
              << global_cell_id.to_string() << "   [machine: " << proc_name
              << " | rank: "
              << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
              << "]   ....." << std::endl;
    timer.restart();

    GridGenerator::general_cell(triangulation,
                                  corner_points,
                                  /* colorize faces */ false);
    triangulation.refine_global(1);
    
    setup_system();
    
    assemble_system();
    
    for (unsigned int q_point = 0;
          q_point < GeometryInfo<dim>::vertices_per_cell;
          ++q_point)
      {
        
        system_rhs.reinit(solution_vector[q_point].size());
        system_matrix.reinit(sparsity_pattern);

        system_matrix.copy_from(assembled_cell_matrix);

        constraints_vector[q_point].condense(system_matrix, system_rhs);
                      
        solve(q_point);
        
      }
    
    assemble_global_element_matrix();

    // pcout << "in cell " << global_cell_id << 
    //     "number of entries: " << global_element_rhs.size() << std::endl;       
    //     for (unsigned int j = 0; j < global_element_rhs.size(); ++j)
    //       if (global_element_rhs(j) != 0)
    //         pcout << j << " : "
    //         << global_element_rhs(j) <<std::endl;

    // pcout << "in cell " << global_cell_id << 
    //     "number of rows: " << global_element_matrix.m() <<
    //     "number of colomns: " << global_element_matrix.n() <<  std::endl;
    //     for (unsigned int i = 0; i < global_element_matrix.m(); ++i)
    //       {
    //         for (unsigned int j = 0; j < global_element_matrix.n(); ++j)
    //           if (global_element_matrix(i,j) != 0)
    //             pcout << "(" << i << "," << j << "): "
    //             << global_element_matrix(i,j) <<std::endl;
          // }

    {
      // Free memory as much as possible
      system_matrix.clear();
      for (unsigned int i = 0; i < GeometryInfo<3>::vertices_per_cell; ++i)
      {
        constraints_vector[i].clear();
      }

      if (true)
      {
        timer.stop();

        std::cout << "done in   " << timer.cpu_time() << "   seconds."
                  << std::endl;
      }
    }
  }

  


  template <int dim>
  const FullMatrix<double> &
  ElaBasis<dim>::get_global_element_matrix() const
  {
    return global_element_matrix;
  }


  template <int dim>
  const Vector<double> &
  ElaBasis<dim>::get_global_element_rhs() const
  {
    return global_element_rhs;
  }
} // namespace Elasticity

#endif // _INCLUDE_ELA_BASIS_TPP_