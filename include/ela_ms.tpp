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
  ElaMs<dim>::ElaMs(const ElaParameters<dim> &ela_parameters)
    : mpi_communicator(MPI_COMM_WORLD)
    , triangulation(mpi_communicator,
                    typename Triangulation<dim>::MeshSmoothing(
                      Triangulation<dim>::smoothing_on_refinement |
                      Triangulation<dim>::smoothing_on_coarsening))
    , fe(FE_Q<dim>(1), dim)
    , dof_handler(triangulation)
    , first_cell_id()
    , cell_basis_map()
    , ela_parameters(ela_parameters)
    , processor_is_used(false)
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

    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<dim>(dim),
                                             constraints);

    if (dim == 3)
      {
        if (ela_parameters.rotate)
          {
            VectorTools::interpolate_boundary_values(
              dof_handler,
              1,
              MyTools::Rotation(ela_parameters.init_p1,
                                ela_parameters.init_p2,
                                ela_parameters.angle),
              constraints);
          }
      }


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
        MyTools::create_data_directory("output/fine_ms_partitioned");
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

    cell_basis_map.clear();

    typename Triangulation<dim>::active_cell_iterator cell = dof_handler
                                                               .begin_active(),
                                                      endc = dof_handler.end();

    for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
          {
            if (!processor_is_used)
              {
                processor_is_used = true;
                std::vector<CellId> first_cells =
                  Utilities::MPI::all_gather(mpi_communicator, cell->id());
                first_cell_id = first_cells[1];
              }
            else
              {
                Assert(first_cell_id.get_coarse_cell_id() !=
                         numbers::invalid_coarse_cell_id,
                       ExcMessage("first cell id not initialized"));
              }

            ElaBasis<dim> current_cell_problem(
              cell,
              first_cell_id,
              triangulation.locally_owned_subdomain(),
              mpi_communicator,
              ela_parameters);

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
    TimerOutput::Scope t(computing_timer, "assembly");

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
    if (ela_parameters.direct_solver_ms)
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

        unsigned int n_iterations = dof_handler.n_dofs();
        const double solver_tolerance =
          std::max(1.e-10, 1e-8 * system_rhs.l2_norm());
        SolverControl solver_control(
          /* n_max_iter */ n_iterations,
          solver_tolerance,
          /* log_history */ true,
          /* log_result */ true);

        TrilinosWrappers::SolverCG solver(solver_control);

        ////////////////////////////////////////////////////////////////
        ///////////////////////////////////
        // TrilinosWrappers::PreconditionSSOR                 preconditioner;
        // TrilinosWrappers::PreconditionSSOR::AdditionalData data(
        //   /* over relaxation */ 1.2);
        ///////////////////////////////////

        ///////////////////////////////////
        std::vector<std::vector<bool>> constant_modes;
        FEValuesExtractors::Vector     displacement_components(0);
        DoFTools::extract_constant_modes(dof_handler,
                                         fe.component_mask(
                                           displacement_components),
                                         constant_modes);
        TrilinosWrappers::PreconditionAMG                 preconditioner;
        TrilinosWrappers::PreconditionAMG::AdditionalData data;
        data.constant_modes        = constant_modes;
        data.elliptic              = true;
        data.higher_order_elements = false;
        data.smoother_sweeps       = 2;
        data.smoother_type         = "ML symmetric Gauss-Seidel";
        data.aggregation_threshold = 0.002;
        ///////////////////////////////////
        ////////////////////////////////////////////////////////////////

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
    GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                    estimated_error_per_cell,
                                                    0.3,
                                                    0.03);
    triangulation.execute_coarsening_and_refinement();
  }


  template <int dim>
  void
  ElaMs<dim>::output_results()
  {
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);

    if (processor_is_used)
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
        StressPostprocessor<dim> stress_postproc(ela_parameters);
        data_out.add_data_vector(locally_relevant_solution, stress_postproc);

        data_out.build_patches();

        // write the output files
        const std::string coarse_filename =
          (std::string("ms_solution.") +
           Utilities::int_to_string(triangulation.locally_owned_subdomain(),
                                    4) +
           std::string(".vtu"));
        std::ofstream output("output/coarse/" + coarse_filename);
        data_out.write_vtu(output);
      }

    std::vector<bool> used_processors =
      Utilities::MPI::gather(mpi_communicator, processor_is_used);

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

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      {
        std::vector<std::string> ordered_basis_filenames;
        for (unsigned int i = 0; i < gathered_basis_filenames.size(); ++i)
          for (unsigned int j = 0; j < gathered_basis_filenames[i].size(); ++j)
            {
              ordered_basis_filenames.push_back(gathered_basis_filenames[i][j]);
            }

        std::vector<std::string> coarse_filenames;
        for (unsigned int i = 0;
             i < Utilities::MPI::n_mpi_processes(mpi_communicator);
             ++i)
          if (used_processors[i])
            coarse_filenames.push_back(std::string("coarse/ms_solution-") +
                                       "." + Utilities::int_to_string(i, 4) +
                                       std::string(".vtu"));

        std::ofstream master_output(std::string("output/ms_solution") +
                                    std::string(".pvtu"));
        data_out.write_pvtu_record(master_output, coarse_filenames);

        std::ofstream fine_master_output(
          std::string("output/fine_ms_solution") + (".pvtu"));
        data_out.write_pvtu_record(fine_master_output, ordered_basis_filenames);
      }
  }


  template <int dim>
  void
  ElaMs<dim>::compute_errors(Vector<double> &coarse_solution,
                             Vector<double> &fine_solution)
  {
    {
      TimerOutput::Scope t(computing_timer, "computing errors");

      parallel::shared::Triangulation<dim> triangulation_fine(
        mpi_communicator,
        typename Triangulation<dim>::MeshSmoothing(
          Triangulation<dim>::smoothing_on_refinement |
          Triangulation<dim>::smoothing_on_coarsening));
      triangulation_fine.copy_triangulation(triangulation);
      triangulation_fine.refine_global(ela_parameters.fine_refinements);
      DoFHandler<dim> dof_handler_fine;
      dof_handler_fine.initialize(triangulation_fine, fe);

      Vector<double> difference_per_cell(triangulation_fine.n_active_cells());

      double L2error_ms, H1error_ms, L2error_coarse, H1error_coarse,
        fine_solution_l2_inv;

      VectorTools::integrate_difference(dof_handler_fine,
                                        fine_solution,
                                        ms_solution_function,
                                        difference_per_cell,
                                        QGauss<dim>(fe.degree + 1),
                                        VectorTools::L2_norm);

      {
        pcout << "Assembling the fine scale MsFEM solution..." << std::endl;
        Vector<double> ms_solution = get_fine_solution(dof_handler_fine);
        output_fine_solution(triangulation_fine, dof_handler_fine, ms_solution);

        pcout << "Computing the L2-error of the MsFEM solution..." << std::endl;
        Functions::FEFieldFunction<dim> ms_solution_function(dof_handler_fine,
                                                             ms_solution);

        VectorTools::integrate_difference(dof_handler_fine,
                                          fine_solution,
                                          ms_solution_function,
                                          difference_per_cell,
                                          QGauss<dim>(fe.degree + 1),
                                          VectorTools::L2_norm);

        L2error_ms = VectorTools::compute_global_error(triangulation_fine,
                                                       difference_per_cell,
                                                       VectorTools::L2_norm);

        pcout
          << "Computing the error of the MsFEM solution in the H1-seminorm..."
          << std::endl;
        VectorTools::integrate_difference(dof_handler_fine,
                                          fine_solution,
                                          ms_solution_function,
                                          difference_per_cell,
                                          QGauss<dim>(fe.degree + 1),
                                          VectorTools::H1_seminorm);

        H1error_ms =
          VectorTools::compute_global_error(triangulation_fine,
                                            difference_per_cell,
                                            VectorTools::H1_seminorm);
      }

      {
        pcout << "Computing the L2-error of the coarse scale standard FEM"
                 " solution..."
              << std::endl;
        Functions::FEFieldFunction<dim> coarse_solution_function(
          dof_handler, coarse_solution);

        VectorTools::integrate_difference(dof_handler_fine,
                                          fine_solution,
                                          coarse_solution_function,
                                          difference_per_cell,
                                          QGauss<dim>(fe.degree + 1),
                                          VectorTools::L2_norm);

        L2error_coarse =
          VectorTools::compute_global_error(triangulation_fine,
                                            difference_per_cell,
                                            VectorTools::L2_norm);

        pcout << "Computing the error of the coarse scale standard FEM"
                 " solution in the H1-seminorm..."
              << std::endl;
        VectorTools::integrate_difference(dof_handler_fine,
                                          fine_solution,
                                          coarse_solution_function,
                                          difference_per_cell,
                                          QGauss<dim>(fe.degree + 1),
                                          VectorTools::H1_seminorm);

        H1error_coarse =
          VectorTools::compute_global_error(triangulation_fine,
                                            difference_per_cell,
                                            VectorTools::H1_seminorm);
      }

      if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
        {
          std::string header("+");
          header += std::string(17, '-');
          header += '+';
          header += std::string(17, '-');
          header += '+';
          header += std::string(17, '-');
          header += '+';

          std::ios oldState(nullptr);
          oldState.copyfmt(std::cout);

          std::cout << "\n\n\n" << header << std::endl;
          std::cout << "| " << std::left << std::setw(15) << "Error"
                    << " | " << std::right << std::setw(15) << "MsFEM"
                    << " | " << std::right << std::setw(15) << "Standard FEM"
                    << " |" << std::endl;
          std::cout << header << std::endl;
          std::cout << "| " << std::left << std::setw(15) << "L2-norm"
                    << " | " << std::right << std::setw(15) << std::scientific
                    << L2error_ms << " | " << std::right << std::setw(15)
                    << L2error_coarse << " |" << std::endl;
          std::cout << "| " << std::left << std::setw(15) << std::scientific
                    << "H1-seminorm"
                    << " | " << std::right << std::setw(15) << H1error_ms
                    << " | " << std::right << std::setw(15) << H1error_coarse
                    << " |" << std::endl;
          std::cout << header << std::endl;

          std::cout.copyfmt(oldState);
        }
    }
  }


  template <int dim>
  const Vector<double>
  ElaMs<dim>::get_fine_solution(DoFHandler<dim> &dof_handler_fine)
  {
    IndexSet locally_relevant_dofs_fine;
    DoFTools::extract_locally_relevant_dofs(dof_handler_fine,
                                            locally_relevant_dofs_fine);

    IndexSet locally_owned_dofs_fine;
    locally_owned_dofs_fine = dof_handler_fine.locally_owned_dofs();

    TrilinosWrappers::MPI::Vector locally_owned_solution_fine;
    locally_owned_solution_fine.reinit(locally_owned_dofs_fine,
                                       mpi_communicator);

    std::vector<types::global_dof_index> local_dof_indices(
      fe.n_dofs_per_cell());

    std::vector<types::global_dof_index> cell_dof_indices;

    std::vector<bool> dof_index_indicator(dof_handler_fine.n_dofs(), false);

    std::map<unsigned int, IndexSet> index_set_to_send;
    for (unsigned int i =
           Utilities::MPI::this_mpi_process(mpi_communicator) + 1;
         i < Utilities::MPI::n_mpi_processes(mpi_communicator);
         ++i)
      {
        index_set_to_send.insert(std::make_pair(i, locally_owned_dofs_fine));
      }

    std::map<unsigned int, IndexSet> locally_owned_dofs_map_other(
      Utilities::MPI::some_to_some(mpi_communicator, index_set_to_send));

    for (const auto &cell : dof_handler_fine.cell_iterators())
      {
        if (cell_basis_map.find(cell->id()) != cell_basis_map.end())
          {
            for (const auto &coarse_cell : dof_handler.active_cell_iterators())
              if (coarse_cell->id() == cell->id())
                {
                  if (coarse_cell->is_locally_owned())
                    {
                      typename std::map<CellId, ElaBasis<dim>>::iterator
                        it_basis = cell_basis_map.find(cell->id());

                      std::vector<Vector<double>> local_solution =
                        (it_basis->second).get_global_solution();

                      unsigned int i = 0;

                      for (const auto &fine_cell :
                           GridTools::get_active_child_cells<DoFHandler<dim>>(
                             cell))
                        {
                          cell_dof_indices.clear();
                          fine_cell->get_dof_indices(local_dof_indices);

                          for (unsigned k = 0; k < local_dof_indices.size();
                               ++k)
                            {
                              if (!dof_index_indicator[local_dof_indices[k]])
                                {
                                  bool is_on_other_proc = false;

                                  for (unsigned int m = 0;
                                       m < Utilities::MPI::this_mpi_process(
                                             mpi_communicator);
                                       ++m)
                                    {
                                      if ((locally_owned_dofs_map_other.at(m))
                                            .is_element(local_dof_indices[k]))
                                        {
                                          is_on_other_proc = true;
                                          break;
                                        }
                                    }
                                  if (!is_on_other_proc)
                                    {
                                      cell_dof_indices.push_back(
                                        local_dof_indices[k]);
                                    }

                                  dof_index_indicator[local_dof_indices[k]] =
                                    true;
                                }
                            }

                          if (cell_dof_indices.size() != 0)
                            {
                              Vector<double> local_cleaned_solution(
                                cell_dof_indices.size());

                              unsigned int l = 0;

                              for (unsigned k = 0; k < local_dof_indices.size();
                                   ++k)
                                {
                                  if (local_dof_indices[k] ==
                                      cell_dof_indices[l])
                                    {
                                      local_cleaned_solution(l) =
                                        local_solution[i](k);
                                      ++l;
                                    }
                                  if (l >= cell_dof_indices.size())
                                    break;
                                }

                              locally_owned_solution_fine.set(
                                cell_dof_indices, local_cleaned_solution);
                            }

                          ++i;
                        }
                    }
                }
          }
      }

    locally_owned_solution_fine.compress(VectorOperation::insert);

    Vector<double> fine_scale_ms_solution(locally_owned_solution_fine);

    //   output_fine_solution(dof_handler_fine, fine_scale_ms_solution);

    return fine_scale_ms_solution;
  }


  template <int dim>
  void
  ElaMs<dim>::output_fine_solution(
    parallel::shared::Triangulation<dim> &triangulation_fine,
    DoFHandler<dim>                      &dof_handler_fine,
    const Vector<double>                 &fine_solution)
  {
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler_fine);

    // add the displacement to the output
    std::vector<std::string> solution_name(dim, "displacement");
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      interpretation(dim,
                     DataComponentInterpretation::component_is_part_of_vector);

    data_out.add_data_vector(fine_solution,
                             solution_name,
                             DataOut<dim>::type_dof_data,
                             interpretation);

    // add the linearized strain tensor to the output
    StrainPostprocessor<dim> strain_postproc;
    data_out.add_data_vector(fine_solution, strain_postproc);

    // add the linearized stress tensor to the output
    StressPostprocessor<dim> stress_postproc(ela_parameters);
    data_out.add_data_vector(fine_solution, stress_postproc);

    data_out.build_patches();

    // write the output files
    std::string filename = "output/fine_ms_partitioned/fine_scale_ms_solution";
    filename +=
      Utilities::int_to_string(triangulation_fine.locally_owned_subdomain(), 4);
    filename += std::string(".vtu");
    std::ofstream output(filename);
    data_out.write_vtu(output);

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      {
        std::vector<std::string> filenames;
        std::string              filename;
        for (unsigned int i = 0;
             i < Utilities::MPI::n_mpi_processes(mpi_communicator);
             ++i)
          {
            filename = "fine_ms_partitioned/fine_scale_ms_solution";
            filename += Utilities::int_to_string(i, 4);
            filename += std::string(".vtu");
            filenames.push_back(filename);
          }

        std::string master_file("output/fine_scale_ms_solution");
        master_file += std::string(".pvtu");
        std::ofstream master_output(master_file);
        data_out.write_pvtu_record(master_output, filenames);
      }
  }


  template <int dim>
  void
  ElaMs<dim>::run()
  {
    pcout << "MsFEM: "
          << "Running with "
          << "Trilinos"
          << " on " << Utilities::MPI::n_mpi_processes(mpi_communicator)
          << " MPI rank(s)..." << std::endl;

    {
      const Point<dim> p1 = ela_parameters.init_p1, p2 = ela_parameters.init_p2;

      const std::vector<unsigned int> repetitions =
        MyTools::get_repetitions(p1, p2);

      GridGenerator::subdivided_hyper_rectangle(
        triangulation, repetitions, p1, p2, true);

      triangulation.refine_global(ela_parameters.coarse_refinements);
    }

    setup_system();

    pcout << "   Number of active cells:       "
          << triangulation.n_global_active_cells() << std::endl
          << "   Number of degrees of freedom: " << dof_handler.n_dofs()
          << std::endl;

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