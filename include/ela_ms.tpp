#ifndef _INCLUDE_ELA_MS_TPP_
#define _INCLUDE_ELA_MS_TPP_

#include "ela_ms.h"

namespace Elasticity
{
  using namespace dealii;

  /****************************************************************************/
  /* Class for the coarse scale part of the multiscale implementation for
     linear elasticity problems */

  // The constructor
  template <int dim>
  ElaMs<dim>::ElaMs(const GlobalParameters<dim> &global_parameters,
                    const ParametersMs &         parameters_ms,
                    const ParametersBasis &      parameters_basis)
    : mpi_communicator(MPI_COMM_WORLD)
    , triangulation(mpi_communicator,
                    typename Triangulation<dim>::MeshSmoothing(
                      Triangulation<dim>::smoothing_on_refinement |
                      Triangulation<dim>::smoothing_on_coarsening))
    , fe(FE_Q<dim>(1), dim)
    , dof_handler(triangulation)
    , cell_basis_map()
    , global_parameters(global_parameters)
    , parameters_ms(parameters_ms)
    , parameters_basis(parameters_basis)
    , pcout(std::cout,
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
    , computing_timer(mpi_communicator,
                      pcout,
                      TimerOutput::summary,
                      TimerOutput::wall_times)
  {}


  template <int dim>
  void
  ElaMs<dim>::setup_system()
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
    unsigned int dirichlet_id = 0;

    // The next part is used if a Dirichlet boundary condition is only applied
    // on a part of a face.

    if (global_parameters.other_dirichlet_id)
      {
        const Point<dim> p1(global_parameters.dirichlet_p1),
          p2(global_parameters.dirichlet_p2);
        MyTools::set_dirichlet_id<dim>(p1, p2, 4, 100, triangulation);
        dirichlet_id = 100;
      }

    VectorTools::interpolate_boundary_values(dof_handler,
                                             dirichlet_id,
                                             Functions::ZeroFunction<dim>(dim),
                                             constraints);

    constraints.close();
    DynamicSparsityPattern dsp(locally_relevant_dofs);
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, true);
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

    // std::filesystem::create_directories("output/basis_output/");
    // std::filesystem::create_directory("output/global_basis_output/");
    // std::filesystem::create_directory("output/coarse/");

    try
      {
        MyTools::create_data_directory("output/");
        MyTools::create_data_directory("output/basis_output/");
        MyTools::create_data_directory("output/global_basis_output/");
        MyTools::create_data_directory("output/coarse/");
      }
    catch (std::runtime_error &e)
      {
        // No exception handling here.
      }


    // preconditioner_matrix.clear();
    // preconditioner_matrix.reinit(locally_owned_dofs, dsp, mpi_communicator);
  }


  template <int dim>
  void
  ElaMs<dim>::initialize_and_compute_basis()
  {
    TimerOutput::Scope t(computing_timer,
                         "basis initialization and computation");

    typename Triangulation<dim>::active_cell_iterator first_cell,
      cell = dof_handler.begin_active(), endc = dof_handler.end();

    unsigned int i = 0;

    for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
          {
            if (i++ == 0)
              first_cell = cell;
            ElaBasis<dim> current_cell_problem(
              cell,
              first_cell,
              triangulation.locally_owned_subdomain(),
              mpi_communicator,
              parameters_basis);

            std::pair<typename std::map<CellId, ElaBasis<dim>>::iterator, bool>
              result;

            result = cell_basis_map.insert(
              std::make_pair(cell->id(), current_cell_problem));

            Assert(result.second,
                   ExcMessage(
                     "Insertion of local basis problem into std::map failed. "
                     "Problem with copy constructor?"));
          }
      } // end ++cell


    /*
     * Now each node possesses a set of basis objects.
     * We need to compute them on each node and do so in
     * a locally threaded way.
     */
    typename std::map<CellId, ElaBasis<dim>>::iterator it_basis =
                                                         cell_basis_map.begin(),
                                                       it_endbasis =
                                                         cell_basis_map.end();
    for (; it_basis != it_endbasis; ++it_basis)
      {
        (it_basis->second).run();
      }
  }


  template <int dim>
  void
  ElaMs<dim>::assemble_system()
  {
    TimerOutput::Scope    t(computing_timer, "assembly");
    const QGauss<dim - 1> face_quadrature_formula(fe.degree + 1);
    FEFaceValues<dim>     fe_face_values(fe,
                                     face_quadrature_formula,
                                     update_values | update_quadrature_points |
                                       update_normal_vectors |
                                       update_JxW_values);

    const unsigned int  n_face_q_points = face_quadrature_formula.size();
    SurfaceForce<dim>   surface_force(global_parameters.surface_force);
    std::vector<double> surface_force_values(n_face_q_points);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);
    Vector<double>     cell_rhs_tmp(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            typename std::map<CellId, ElaBasis<dim>>::iterator it_basis =
              cell_basis_map.find(cell->id());

            cell_matrix = 0.;
            cell_rhs    = 0.;

            cell_matrix = (it_basis->second).get_global_element_matrix();
            cell_rhs    = (it_basis->second).get_global_element_rhs();

            if (global_parameters.neumann_bc)
              {
                for (const auto &face : cell->face_iterators())
                  if (face->at_boundary() && (face->boundary_id() == 1))
                    {
                      std::vector<double> surface_force_values(n_face_q_points);
                      fe_face_values.reinit(cell, face);
                      surface_force.value_list(
                        fe_face_values.get_quadrature_points(),
                        surface_force_values);
                      for (unsigned int q_point = 0; q_point < n_face_q_points;
                           ++q_point)
                        {
                          for (unsigned int i = 0; i < dofs_per_cell; ++i)
                            cell_rhs(i) +=
                              (fe_face_values.shape_value(
                                 i,
                                 q_point) *                    // phi_i(x_q)
                               surface_force_values[q_point] * // g(x_q)
                               fe_face_values.JxW(q_point));   // dx
                        }
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
  ElaMs<dim>::solve()
  {
    if (parameters_ms.direct_solver)
      {
        TimerOutput::Scope t(computing_timer,
                             "parallel sparse direct solver (MUMPS)");

        if (parameters_ms.verbose)
          {
            pcout << "   Using direct solver..." << std::endl;
          }

        TrilinosWrappers::MPI::Vector completely_distributed_solution(
          locally_owned_dofs, mpi_communicator);
        SolverControl                  solver_control;
        TrilinosWrappers::SolverDirect solver(solver_control);

        solver.initialize(system_matrix);
        solver.solve(system_matrix,
                     completely_distributed_solution,
                     system_rhs);

        if (parameters_ms.verbose)
          {
            pcout << "   Solved in with parallel sparse direct solver (MUMPS)."
                  << std::endl;
          }

        constraints.distribute(completely_distributed_solution);

        locally_relevant_solution = completely_distributed_solution;
      }
    else
      {
        TimerOutput::Scope t(computing_timer, "solve");

        if (parameters_ms.verbose)
          {
            pcout << "   Using iterative solver..." << std::endl;
          }

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

        if (parameters_ms.verbose)
          {
            pcout << "   Solved (iteratively) in " << solver_control.last_step()
                  << " iterations." << std::endl;
          }

        constraints.distribute(completely_distributed_solution);
        locally_relevant_solution = completely_distributed_solution;
      }
  }


  template <int dim>
  void
  ElaMs<dim>::send_global_weights_to_cell()
  {
    // For each cell we get dofs_per_cell values
    const unsigned int                   dofs_per_cell = fe.dofs_per_cell;
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell =
                                                     dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
          {
            cell->get_dof_indices(local_dof_indices);
            std::vector<double> extracted_weights(dofs_per_cell, 0);
            locally_relevant_solution.extract_subvector_to(local_dof_indices,
                                                           extracted_weights);

            typename std::map<CellId, ElaBasis<dim>>::iterator it_basis =
              cell_basis_map.find(cell->id());
            (it_basis->second).set_global_weights(extracted_weights);
          }
      } // end ++cell
  }


  // Adaptive refinement of the grid
  template <int dim>
  void
  ElaMs<dim>::refine_grid()
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
  ElaMs<dim>::output_results()
  {
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);

    if (locally_relevant_solution.size() != 0)
      {
        // add the displacement to the output
        std::vector<std::string> solution_name(dim, "displacement");
        std::vector<DataComponentInterpretation::DataComponentInterpretation>
          interpretation(
            dim, DataComponentInterpretation::component_is_part_of_vector);

        data_out.add_data_vector(locally_relevant_solution,
                                 solution_name,
                                 DataOut<dim>::type_dof_data,
                                 interpretation);

        // add the linearized strain tensor to the output
        StrainPostprocessor<dim> strain_postproc;
        data_out.add_data_vector(locally_relevant_solution, strain_postproc);

        // add the linearized stress tensor to the output
        StressPostprocessor<dim> stress_postproc(global_parameters);
        data_out.add_data_vector(locally_relevant_solution, stress_postproc);

        data_out.build_patches();

        // write the output files
        const std::string coarse_filename =
          ("ms_solution-" +
           Utilities::int_to_string(triangulation.locally_owned_subdomain(),
                                    4) +
           ".vtu");
        std::ofstream output("output/coarse/" + coarse_filename);
        data_out.write_vtu(output);
      }

    std::vector<std::string> basis_filenames;

    typename std::map<CellId, ElaBasis<dim>>::iterator it_basis =
                                                         cell_basis_map.begin(),
                                                       it_endbasis =
                                                         cell_basis_map.end();

    for (; it_basis != it_endbasis; ++it_basis)
      {
        (it_basis->second).output_global_solution_in_cell();
        basis_filenames.push_back((it_basis->second).get_filename());
      }

    std::vector<std::vector<std::string>> gathered_basis_filenames =
      Utilities::MPI::gather(mpi_communicator, basis_filenames);

    std::vector<std::string> ordered_basis_filenames;
    for (unsigned int i = 0; i < gathered_basis_filenames.size(); ++i)
      for (unsigned int j = 0; j < gathered_basis_filenames[i].size(); ++j)
        {
          ordered_basis_filenames.push_back(gathered_basis_filenames[i][j]);
        }

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      {
        std::vector<std::string> coarse_filenames;
        std::string              tmp_filename, tmp_filename2;
        struct stat              info;
        for (unsigned int i = 0;
             i < Utilities::MPI::n_mpi_processes(mpi_communicator);
             ++i)
          {
            tmp_filename =
              "coarse/ms_solution-" + Utilities::int_to_string(i, 4) + ".vtu";
            tmp_filename2            = "output/" + tmp_filename;
            const char *tmp_filechar = tmp_filename2.c_str();
            if (stat(tmp_filechar, &info) == 0)
              coarse_filenames.push_back(tmp_filename);
          }
        std::ofstream master_output("output/ms_solution.pvtu");
        data_out.write_pvtu_record(master_output, coarse_filenames);

        std::ofstream fine_master_output("output/fine_ms_solution.pvtu");
        data_out.write_pvtu_record(fine_master_output, ordered_basis_filenames);
      }
  }


  template <int dim>
  void
  ElaMs<dim>::run()
  {
    if (parameters_ms.verbose)
      {
        pcout << "Running with "
              << "Trilinos"
              << " on " << Utilities::MPI::n_mpi_processes(mpi_communicator)
              << " MPI rank(s)..." << std::endl;
      }

    const Point<dim> p1 = global_parameters.init_p1,
                     p2 = global_parameters.init_p2;

    const std::vector<unsigned int> repetitions =
      MyTools::get_repetitions(p1, p2);

    GridGenerator::subdivided_hyper_rectangle(
      triangulation, repetitions, p1, p2, true);

    // GridGenerator::cylinder(triangulation, 10., 0.1);

    triangulation.refine_global(parameters_ms.n_refine);

    // GridTools::transform(MyTools::roty<dim>, triangulation);
    setup_system();

    if (parameters_ms.verbose)
      {
        pcout << "   Number of active cells:       "
              << triangulation.n_global_active_cells() << std::endl
              << "   Number of degrees of freedom: " << dof_handler.n_dofs()
              << std::endl;
      }

    initialize_and_compute_basis();

    assemble_system();

    solve();

    send_global_weights_to_cell();

    {
      TimerOutput::Scope t(computing_timer, "output");
      output_results();
    }

    computing_timer.print_summary();
    computing_timer.reset();
    pcout << std::endl;
  }
} // namespace Elasticity

#endif // _INCLUDE_ELA_MS_TPP_