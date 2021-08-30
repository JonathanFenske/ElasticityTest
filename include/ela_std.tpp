#ifndef _INCLUDE_ELA_STD_TPP_
#define _INCLUDE_ELA_STD_TPP_

#include "ela_std.h"

namespace Elasticity
{
  using namespace dealii;

  /****************************************************************************/
  /* The main class*/

  // The constructor
  template <int dim>
  ElaStd<dim>::ElaStd()
    : mpi_communicator(MPI_COMM_WORLD)
    , triangulation(mpi_communicator,
                    typename Triangulation<dim>::MeshSmoothing(
                      Triangulation<dim>::smoothing_on_refinement |
                      Triangulation<dim>::smoothing_on_coarsening))
    , fe(FE_Q<dim>(1), dim)
    , dof_handler(triangulation)
    , pcout(std::cout,
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
    , computing_timer(mpi_communicator,
                      pcout,
                      TimerOutput::summary,
                      TimerOutput::wall_times)
    , direct_solver(
        true) // If true, a direct solver will be used to solve the problem.
  {}


  template <int dim>
  std::tuple<Point<dim>, Point<dim>>
  ElaStd<dim>::get_init_vert(std::vector<double> p) const
  {
    Point<dim> p1_tmp, p2_tmp;
    for (int i = 0; i < dim; ++i)
      p1_tmp(i) = p[i], p2_tmp(i) = p[dim + i];
    return {p1_tmp, p2_tmp};
  }


  template <int dim>
  void
  ElaStd<dim>::setup_system()
  {
    TimerOutput::Scope t(computing_timer, "setup");
    dof_handler.distribute_dofs(fe);
    locally_owned_dofs = dof_handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
    locally_relevant_solution.reinit(locally_owned_dofs,
                                     locally_relevant_dofs,
                                     mpi_communicator);
    system_rhs.reinit(locally_owned_dofs, mpi_communicator);
    constraints.clear();
    constraints.reinit(locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);

    // The next part is used if a Dirichlet boundary condition is only applied
    // on a part of a face.

    // {
    //   Point<dim> p1 = {-1, -1, -0.1}, p2 = {1, 1, 0.1};
    //   MyTools::set_dirichlet_id(p1, p2, 4, 100);
    // }

    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<dim>(dim),
                                             constraints);

    constraints.close();
    DynamicSparsityPattern dsp(locally_relevant_dofs);
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
    SparsityTools::distribute_sparsity_pattern(
      dsp,
      Utilities::MPI::all_gather(mpi_communicator,
                                 dof_handler.n_locally_owned_dofs()),
      mpi_communicator,
      locally_relevant_dofs);
    system_matrix.reinit(locally_owned_dofs,
                         locally_owned_dofs,
                         dsp,
                         mpi_communicator);

    preconditioner_matrix.clear();
    preconditioner_matrix.reinit(locally_owned_dofs, dsp, mpi_communicator);
  }


  template <int dim>
  void
  ElaStd<dim>::assemble_system()
  {
    TimerOutput::Scope    t(computing_timer, "assembly");
    const QGauss<dim>     quadrature_formula(fe.degree + 1);
    const QGauss<dim - 1> face_quadrature_formula(fe.degree + 1);
    FEValues<dim>         fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);
    FEFaceValues<dim>     fe_face_values(fe,
                                     face_quadrature_formula,
                                     update_values | update_quadrature_points |
                                       update_normal_vectors |
                                       update_JxW_values);
    const unsigned int    dofs_per_cell   = fe.n_dofs_per_cell();
    const unsigned int    n_q_points      = quadrature_formula.size();
    const unsigned int    n_face_q_points = face_quadrature_formula.size();
    std::vector<double>   lambda_values(n_q_points), mu_values(n_q_points);
    std::vector<Vector<double>>   body_force_values(n_q_points);
    for (unsigned int i = 0; i < n_q_points; ++i)
      body_force_values[i].reinit(dim);
    lambda<dim>           lambda;
    mu<dim>               mu;
    BodyForce<dim>        body_force;
    SurfaceForce<dim>     surface_force;
    FullMatrix<double>    cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>        cell_rhs(dofs_per_cell);
    Vector<double>        cell_rhs_tmp(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            cell_matrix = 0.;
            cell_rhs    = 0.;
            fe_values.reinit(cell);
            lambda.value_list(fe_values.get_quadrature_points(), lambda_values);
            mu.value_list(fe_values.get_quadrature_points(), mu_values);
            body_force.vector_value_list(fe_values.get_quadrature_points()
                                          , body_force_values);
            for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
              {
                for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
                  {
                    const unsigned int component_i =
                      fe.system_to_component_index(i).first;
                    for (unsigned int j = 0; j < fe.dofs_per_cell; ++j)
                      {
                        const unsigned int component_j =
                          fe.system_to_component_index(j).first;
                        cell_matrix(i, j) +=
                          ((fe_values.shape_grad(i, q_point)[component_i] *
                            fe_values.shape_grad(j, q_point)[component_j] *
                            lambda_values[q_point]) +
                           (fe_values.shape_grad(i, q_point)[component_j] *
                            fe_values.shape_grad(j, q_point)[component_i] *
                            mu_values[q_point]) +
                           ((component_i == component_j) ?
                              (fe_values.shape_grad(i, q_point) *
                               fe_values.shape_grad(j, q_point) *
                               mu_values[q_point]) :
                              0)) *
                          fe_values.JxW(q_point);
                      }
                    cell_rhs(i) +=
                      fe_values.shape_value_component(i, q_point, component_i) *
                      body_force_values[q_point][component_i] *
                      fe_values.JxW(q_point);
                  }
                cell_rhs_tmp = cell_rhs;
              }
            for (const auto &face : cell->face_iterators())
              if (face->at_boundary() && (face->boundary_id() == 1))
                {
                  fe_face_values.reinit(cell, face);
                  for (unsigned int q_point = 0; q_point < n_face_q_points;
                       ++q_point)
                    {
                      const double surface_force_value = surface_force.value(
                        fe_face_values.quadrature_point(q_point)); // *
                      // fe_face_values.normal_vector(q_point));
                      for (unsigned int i = 0; i < dofs_per_cell; ++i)
                        cell_rhs(i) +=
                          (fe_face_values.shape_value(i,
                                                      q_point) * // phi_i(x_q)
                           surface_force_value *                 // g(x_q)
                           fe_face_values.JxW(q_point));         // dx
                    }
                }
            cell->get_dof_indices(local_dof_indices);
            constraints.distribute_local_to_global(cell_matrix,
                                                   cell_rhs,
                                                   local_dof_indices,
                                                   system_matrix,
                                                   system_rhs);
          }
      }
    system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);
  }


  template <int dim>
  void
  ElaStd<dim>::solve()
  {
    if (direct_solver)
      {
        TimerOutput::Scope t(computing_timer,
                             "parallel sparse direct solver (MUMPS)");

        pcout << "   Using direct solver..." << std::endl;

        TrilinosWrappers::MPI::Vector completely_distributed_solution(
          locally_owned_dofs, mpi_communicator);
        SolverControl                  solver_control;
        TrilinosWrappers::SolverDirect solver(solver_control);
        solver.initialize(system_matrix);
        solver.solve(system_matrix,
                     completely_distributed_solution,
                     system_rhs);

        pcout << "   Solved in with parallel sparse direct solver (MUMPS)."
              << std::endl;

        constraints.distribute(completely_distributed_solution);

        locally_relevant_solution = completely_distributed_solution;
      }
    else
      {
        TimerOutput::Scope t(computing_timer, "solve");

        pcout << "   Using iterative solver..." << std::endl;

        TrilinosWrappers::MPI::Vector completely_distributed_solution(
          locally_owned_dofs, mpi_communicator);

        unsigned int  n_iterations     = dof_handler.n_dofs();
        const double  solver_tolerance = 1e-8 * system_rhs.l2_norm();
        SolverControl solver_control(
          /* n_max_iter */ n_iterations,
          solver_tolerance,
          /* log_history */ true,
          /* log_result */ true);

        TrilinosWrappers::SolverCG solver(solver_control);

        TrilinosWrappers::PreconditionSSOR                 preconditioner;
        TrilinosWrappers::PreconditionSSOR::AdditionalData data(
          /* over relaxation */ 1.2);

        // TrilinosWrappers::PreconditionBlockJacobi preconditioner;
        // TrilinosWrappers::PreconditionBlockJacobi::AdditionalData data;

        // TrilinosWrappers::PreconditionAMG                 preconditioner;
        // TrilinosWrappers::PreconditionAMG::AdditionalData data;

        preconditioner.initialize(system_matrix, data);

        try
          {
            solver.solve(system_matrix,
                         completely_distributed_solution,
                         system_rhs,
                         preconditioner);
          }
        catch (std::exception &e)
          {
            Assert(false, ExcMessage(e.what()));
          }

        pcout << "   Solved (iteratively) in " << solver_control.last_step()
              << " iterations." << std::endl;
        constraints.distribute(completely_distributed_solution);
        locally_relevant_solution = completely_distributed_solution;
      }
  }


  // Adaptive refinement of the grid
  template <int dim>
  void
  ElaStd<dim>::refine_grid()
  {
    TimerOutput::Scope t(computing_timer, "refine");
    Vector<float>      estimated_error_per_cell(triangulation.n_active_cells());
    KellyErrorEstimator<dim>::estimate(
      dof_handler,
      QGauss<dim - 1>(fe.degree + 1),
      std::map<types::boundary_id, const Function<dim> *>(),
      locally_relevant_solution,
      estimated_error_per_cell);
    parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(
      triangulation, estimated_error_per_cell, 0.3, 0.03);
    triangulation.execute_coarsening_and_refinement();
  }


  template <int dim>
  void
  ElaStd<dim>::output_results(const unsigned int cycle) const
  {
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);

    // add the displacement to the output
    std::vector<std::string> solution_name(dim, "displacement");
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      interpretation(dim,
                     DataComponentInterpretation::component_is_part_of_vector);

    data_out.add_data_vector(locally_relevant_solution,
                             solution_name,
                             DataOut<dim>::type_dof_data,
                             interpretation);

    // add the partition of the domain to the output
    Vector<float> subdomain(triangulation.n_active_cells());
    for (unsigned int i = 0; i < subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();
    data_out.add_data_vector(subdomain, "subdomain");

    // add the linearized strain tensor to the output
    StrainPostprocessor<dim> strain_postproc;
    data_out.add_data_vector(locally_relevant_solution, strain_postproc);

    // add the linearized stress tensor to the output
    StressPostprocessor<dim> stress_postproc;
    data_out.add_data_vector(locally_relevant_solution, stress_postproc);

    data_out.build_patches();

    // write the output files
    const std::string filename =
      ("solution-" + Utilities::int_to_string(cycle, 2) + "." +
       Utilities::int_to_string(triangulation.locally_owned_subdomain(), 4) +
       ".vtu");
    std::ofstream output(filename);
    data_out.write_vtu(output);

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      {
        std::vector<std::string> filenames;
        for (unsigned int i = 0;
             i < Utilities::MPI::n_mpi_processes(mpi_communicator);
             ++i)
          filenames.push_back("solution-" + Utilities::int_to_string(cycle, 2) +
                              "." + Utilities::int_to_string(i, 4) + ".vtu");

        std::ofstream master_output(
          "solution-" + Utilities::int_to_string(cycle, 2) + ".pvtu");
        data_out.write_pvtu_record(master_output, filenames);
      }
  }


  template <int dim>
  void
  ElaStd<dim>::run()
  {
    pcout << "Running with "
          << "Trilinos"
          << " on " << Utilities::MPI::n_mpi_processes(mpi_communicator)
          << " MPI rank(s)..." << std::endl;

    const auto [p1, p2]         = get_init_vert({-10, 0, 0, 10, 1, 1});
    const unsigned int n_cycles = 3;
    for (unsigned int cycle = 0; cycle < n_cycles; ++cycle)
      {
        pcout << "Cycle " << cycle << ':' << std::endl;
        if (cycle == 0)
          {
            std::list<double> side_length_list;
            for (unsigned int i = 0; i < dim; ++i)
              side_length_list.push_back(p2[i] - p1[i]);
            double cell_length = *(std::min_element(side_length_list.begin(),
                                                    side_length_list.end()));
            std::vector<unsigned int> repetitions;
            for (auto it : side_length_list)
              repetitions.push_back((unsigned int)(it / cell_length));

            GridGenerator::subdivided_hyper_rectangle(
              triangulation, repetitions, p1, p2, true);

            // GridGenerator::cylinder(triangulation, 10., 0.1);

            triangulation.refine_global(1);
          }
        else
          {
            // GridTools::transform(MyTools::backroty<dim>, triangulation);
            // refine_grid();
            triangulation.refine_global(1);
          }

        // GridTools::transform(MyTools::roty<dim>, triangulation);
        setup_system();
        pcout << "   Number of active cells:       "
              << triangulation.n_global_active_cells() << std::endl
              << "   Number of degrees of freedom: " << dof_handler.n_dofs()
              << std::endl;
        assemble_system();
        solve();
        if (Utilities::MPI::n_mpi_processes(mpi_communicator) <= 32)
          {
            TimerOutput::Scope t(computing_timer, "output");
            output_results(cycle);
          }
        computing_timer.print_summary();
        computing_timer.reset();
        pcout << std::endl;
      }
  }
} // namespace Elasticity

#endif // _INCLUDE_ELA_STD_TPP_