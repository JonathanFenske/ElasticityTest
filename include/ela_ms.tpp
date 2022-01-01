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
        MyTools::create_data_directory("output/coarse_ms_partitioned/");
        MyTools::create_data_directory("output/other_partitioned/");
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

    pcout << "cell_basis_map.clear();" << std::endl;
    cell_basis_map.clear();

    pcout << "create iterators" << std::endl;
    typename Triangulation<dim>::active_cell_iterator cell = dof_handler
                                                               .begin_active(),
                                                      endc = dof_handler.end();

    pcout << "iterating..." << std::endl;

    for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
          {
            processor_is_used = true;
            break;
          }
      }

    if (!processor_is_used)
      {
        std::cout << Utilities::MPI::this_mpi_process(mpi_communicator)
                  << " not used" << std::endl;
        std::cout << "pairing" << std::endl;
        std::pair<CellId, bool> tmp_pair(dof_handler.begin_active(),
                                         processor_is_used);
        std::cout << "allocating" << std::endl;
        std::vector<std::pair<CellId, bool>> first_cells(
          Utilities::MPI::n_mpi_processes(mpi_communicator));
        std::cout << "before gathering" << std::endl;
        // first_cells = Utilities::MPI::all_gather(mpi_communicator, tmp_pair);
      }


    // pcout << "second iterating..." << std::endl;
    // for (auto &first_cell_it : first_cells)
    //   {
    //     if (first_cell_it.second)
    //       {
    //         first_cell_id = first_cell_it.first;
    //         break;
    //       }
    //   }
    if (!processor_is_used)
      {
        std::cout << "before first cell" << std::endl;
      }
    first_cell_id = cell->id();

    cell = dof_handler.begin_active();

    // pcout << "third iterating..." << std::endl;
    for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
          {
            // pcout << "initializing cell problem" << std::endl;
            ElaBasis<dim> current_cell_problem(
              cell,
              first_cell_id,
              triangulation.locally_owned_subdomain(),
              mpi_communicator,
              ela_parameters);

            // std::cout << "make pair" << std::endl;
            std::pair<typename std::map<CellId, ElaBasis<dim>>::iterator, bool>
              result;

            // std::cout << "insert" << std::endl;
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
    if (processor_is_used)
      {
        typename std::map<CellId, ElaBasis<dim>>::iterator
          it_basis    = cell_basis_map.begin(),
          it_endbasis = cell_basis_map.end();

        for (; it_basis != it_endbasis; ++it_basis)
          {
            (it_basis->second).run();
          }
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
          (std::string("coarse_ms_solution.") +
           Utilities::int_to_string(triangulation.locally_owned_subdomain(),
                                    4) +
           std::string(".vtu"));
        std::ofstream output("output/coarse_ms_partitioned/" + coarse_filename);
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
            coarse_filenames.push_back(
              std::string("coarse_ms_partitioned/coarse_ms_solution-") + "." +
              Utilities::int_to_string(i, 4) + std::string(".vtu"));

        std::ofstream master_output(std::string("output/coarse_ms_solution") +
                                    std::string(".pvtu"));
        data_out.write_pvtu_record(master_output, coarse_filenames);

        std::ofstream fine_master_output(
          std::string("output/fine_ms_solution") + (".pvtu"));
        data_out.write_pvtu_record(fine_master_output, ordered_basis_filenames);
      }
  }


  template <int dim>
  void
  ElaMs<dim>::compute_errors(
    Vector<double>                                         &coarse_solution,
    Vector<double>                                         &fine_solution,
    std::map<CellId, std::vector<types::global_dof_index>> &dof_map_coarse,
    std::map<CellId, std::vector<types::global_dof_index>> &dof_map_fine)
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

      double L2error_ms, H1error_ms, L2error_coarse, H1error_coarse;

      {
        pcout << "Assembling the fine scale MsFEM solution..." << std::endl;
        Vector<double> ms_solution = get_fine_solution(dof_handler_fine);
        output_fine_solution(triangulation_fine,
                             dof_handler_fine,
                             ms_solution,
                             "fine_ms_assembled");

        pcout << "\nComputing the L2-error of the MsFEM solution..."
              << std::endl;
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
        pcout << "\nReordering coarse scale standard FEM solution..."
              << std::endl;
        coarse_solution =
          map_to_local_mesh(dof_handler, coarse_solution, dof_map_coarse);
        pcout << "Reordering coarse fine standard FEM solution..." << std::endl;
        fine_solution =
          map_to_local_mesh(dof_handler_fine, fine_solution, dof_map_fine);
        pcout << "\nComputing the L2-error of the coarse scale standard FEM"
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

        pcout << "output solutions..." << std::endl;
        output_fine_solution(triangulation_fine,
                             dof_handler_fine,
                             fine_solution,
                             "fine_std");
        output_fine_solution(triangulation,
                             dof_handler,
                             coarse_solution,
                             "coarse_std");
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
  Vector<double>
  ElaMs<dim>::map_to_local_mesh(
    const DoFHandler<dim> &local_dof_handler,
    const Vector<double>  &solution_vector,
    const std::map<CellId, std::vector<types::global_dof_index>> &dof_map)
  {
    pcout << solution_vector.size() << std::endl;
    IndexSet current_locally_owned_dofs;
    current_locally_owned_dofs = local_dof_handler.locally_owned_dofs();
    TrilinosWrappers::MPI::Vector vector_parallel(current_locally_owned_dofs,
                                                  mpi_communicator);

    unsigned int dofs_per_cell = fe.n_dofs_per_cell();

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    std::vector<double>                  local_solution(dofs_per_cell);
    std::vector<types::global_dof_index> other_local_dof_indices;


    for (auto &cell : local_dof_handler.active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            other_local_dof_indices = dof_map.at(cell->id());
            solution_vector.extract_subvector_to(other_local_dof_indices,
                                                 local_solution);
            cell->get_dof_indices(local_dof_indices);

            vector_parallel.set(local_dof_indices, local_solution);
          }
      }

    vector_parallel.compress(VectorOperation::insert);

    return Vector<double>(vector_parallel);
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

    for (const auto &coarse_cell : dof_handler.active_cell_iterators())
      {
        if (coarse_cell->is_locally_owned())
          {
            for (const auto &cell : dof_handler_fine.cell_iterators())
              {
                if (coarse_cell->id() == cell->id())
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
                        fine_cell->get_dof_indices(local_dof_indices);

                        locally_owned_solution_fine.set(local_dof_indices,
                                                        local_solution[i]);

                        ++i;
                      }
                  }
              }
          }
      }

    locally_owned_solution_fine.compress(VectorOperation::insert);

    Vector<double> fine_scale_ms_solution(locally_owned_solution_fine);

    return fine_scale_ms_solution;
  }


  template <int dim>
  void
  ElaMs<dim>::output_fine_solution(
    parallel::shared::Triangulation<dim> &triangulation_fine,
    DoFHandler<dim>                      &dof_handler_fine,
    const Vector<double>                 &fine_solution,
    std::string                           name)
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
    std::string filename = "output/other_partitioned/";
    filename += name;
    filename += std::string("_solution");
    filename +=
      Utilities::int_to_string(triangulation_fine.locally_owned_subdomain(), 4);
    filename += std::string(".vtu");
    std::ofstream output(filename);
    data_out.write_vtu(output);

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      {
        std::vector<std::string> filenames;
        for (unsigned int i = 0;
             i < Utilities::MPI::n_mpi_processes(mpi_communicator);
             ++i)
          {
            filename = "other_partitioned/";
            filename += name;
            filename += std::string("_solution");
            filename += Utilities::int_to_string(i, 4);
            filename += std::string(".vtu");
            filenames.push_back(filename);
          }

        std::string master_file("output/");
        master_file += name;
        master_file += std::string("_solution");
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
    pcout << "assemble" << std::endl;
    assemble_system();
    pcout << "solve" << std::endl;
    solve();

    send_global_weights_to_cell();

    {
      TimerOutput::Scope t(computing_timer, "output");
      pcout << "really" << std::endl;
      output_results();
    }

    pcout << "done" << std::endl;

    computing_timer.print_summary();
    computing_timer.reset();
    pcout << std::endl;
  }
} // namespace Elasticity

#endif // _INCLUDE_ELA_MS_TPP_