// Deal.II
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/logstream.h>
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

// STL
#include <cmath>
#include <fstream>
#include <iostream>



namespace elasticity
{
  using namespace dealii;

  double E  = 210.e9;
  double nu = 0.3;

  template <int dim>
  class SurfaceForce : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;
  };

  template <int dim>
  double
  SurfaceForce<dim>::value(const Point<dim> &p, const unsigned int) const
  {
    // if(p(0) > 4)
    //   return -100;//-100 *fabs(std::sin(M_PI*p(0)/5));
    // else
    //   return 0.;
    return 0;
  }


  template <int dim>
  class lambda : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;
  };


  template <int dim>
  double
  lambda<dim>::value(const Point<dim> &p, const unsigned int) const
  {
    return E / (2 * (1 + nu));
    // return -(std::sin(M_PI*p(0)/15)+2);//(-0.1 * p(0) + 2.5) * 1e9;////*;
  }


  template <int dim>
  class mu : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;
  };


  template <int dim>
  double
  mu<dim>::value(const Point<dim> &p, const unsigned int) const
  {
    int fr = 10;
    return E * nu / ((1 + nu) * (1 - 2 * nu)) *
           (0.1 * std::sin(2 * fr * M_PI * p(0) / 20) + 1); //(-0.025*p(0)+0.75)
  }


  template <int dim>
  const Point<dim>
  roty(const Point<dim> &p)
  {
    Point<dim> q = p;
    if (dim == 3)
      {
        // Point<dim> axis = {0,1,0};
        // Tensor<2,dim> rot =
        // Physics::Transformations::Rotations::rotation_matrix_3d(axis,-M_PI);
        // return Point<dim>(rot*q);
        return Point<dim>(-q(2), q(1), q(0));
      }
    else
      {
        return q;
      }
  }


  template <int dim>
  const Point<dim>
  backroty(const Point<dim> &p)
  {
    Point<dim> q = p;
    if (dim == 3)
      {
        // Point<dim> axis = {0,1,0};
        // Tensor<2,dim> rot =
        // Physics::Transformations::Rotations::rotation_matrix_3d(axis,M_PI);
        // return Point<dim>(rot*q);
        return Point<dim>(q(2), q(1), -q(0));
      }
    else
      {
        return q;
      }
  }


  template <int dim>
  class ElaProblem
  {
  public:
    ElaProblem();
    void
    run();

  private:
    std::tuple<Point<dim>, Point<dim>>
    get_init_vert(std::vector<double> p) const;
    void
    setup_system();
    void
    assemble_system();
    void
    solve();
    void
    refine_grid();
    void
                                              output_results(const unsigned int cycle) const;
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
    ConditionalOStream                        pcout;
    TimerOutput                               computing_timer;
    const bool                                direct_solver;
  };


  template <int dim>
  ElaProblem<dim>::ElaProblem()
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
    , direct_solver(false)
  {}


  template <int dim>
  std::tuple<Point<dim>, Point<dim>>
  ElaProblem<dim>::get_init_vert(std::vector<double> p) const
  {
    Point<dim> p1_tmp, p2_tmp;
    for (int i = 0; i < dim; ++i)
      p1_tmp(i) = p[i], p2_tmp(i) = p[dim + i];
    return {p1_tmp, p2_tmp};
  }


  template <int dim>
  void
  ElaProblem<dim>::setup_system()
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

    // for (auto &face : triangulation.active_face_iterators())
    //   {
    //     if (face->at_boundary() && (face->center()[0] < 1) &&
    //     (face->center()[0] > -1)
    //           && (face->center()[1] < 1) && (face->center()[1] > -1) &&
    //           (face->boundary_id() == 1))
    //     {
    //       face->set_boundary_id(100);
    //     }
    //   }


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
  ElaProblem<dim>::assemble_system()
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
    std::vector<double>   lambda_values(n_q_points);
    std::vector<double>   mu_values(n_q_points);
    lambda<dim>           lambda;
    mu<dim>               mu;
    FullMatrix<double>    cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>        cell_rhs(dofs_per_cell);
    Vector<double>        cell_rhs_tmp(dofs_per_cell);
    const double          rhs_value = -9.81 * 7.85e3;
    SurfaceForce<dim>     surface_force;
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
                      ((component_i == dim - 1) ? rhs_value : 0.) *
                      fe_values.JxW(q_point);
                  }
                cell_rhs_tmp = cell_rhs;
              }
            for (unsigned int face_number = 0;
                 face_number < GeometryInfo<dim>::faces_per_cell;
                 ++face_number)
              if (cell->face(face_number)->at_boundary() &&
                  (cell->face(face_number)->boundary_id() == 4))
                {
                  fe_face_values.reinit(cell, face_number);
                  for (unsigned int q_point = 0; q_point < n_face_q_points;
                       ++q_point)
                    {
                      const double neumann_value = surface_force.value(
                        fe_face_values.quadrature_point(q_point)); // *
                      // fe_face_values.normal_vector(q_point));
                      for (unsigned int i = 0; i < dofs_per_cell; ++i)
                        cell_rhs(i) +=
                          (fe_face_values.shape_value(i,
                                                      q_point) * // phi_i(x_q)
                           neumann_value *                       // g(x_q)
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
  ElaProblem<dim>::solve()
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


  template <int dim>
  void
  ElaProblem<dim>::refine_grid()
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
  ElaProblem<dim>::output_results(const unsigned int cycle) const
  {
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    std::vector<std::string> solution_names;
    switch (dim)
      {
        case 1:
          solution_names.emplace_back("displacement");
          break;
        case 2:
          solution_names.emplace_back("x_displacement");
          solution_names.emplace_back("y_displacement");
          break;
        case 3:
          solution_names.emplace_back("x_displacement");
          solution_names.emplace_back("y_displacement");
          solution_names.emplace_back("z_displacement");
          break;
        default:
          Assert(false, ExcNotImplemented());
      }

    data_out.add_data_vector(locally_relevant_solution, solution_names);

    std::vector<std::string> solution_name(dim, "displacement");
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      interpretation(dim,
                     DataComponentInterpretation::component_is_part_of_vector);

    data_out.add_data_vector(locally_relevant_solution,
                             solution_name,
                             DataOut<dim>::type_dof_data,
                             interpretation);

    Vector<float> subdomain(triangulation.n_active_cells());
    for (unsigned int i = 0; i < subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();
    data_out.add_data_vector(subdomain, "subdomain");

    data_out.build_patches();

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
  ElaProblem<dim>::run()
  {
    pcout << "Running with "
          << "Trilinos"
          << " on " << Utilities::MPI::n_mpi_processes(mpi_communicator)
          << " MPI rank(s)..." << std::endl;

    const auto [p1, p2]         = get_init_vert({-10, 0, 0, 10, 1, 1});
    const unsigned int n_cycles = 1;
    for (unsigned int cycle = 0; cycle < n_cycles; ++cycle)
      {
        pcout << "Cycle " << cycle << ':' << std::endl;
        if (cycle == 0)
          {
            const std::vector<unsigned int> repetitions = {20, 1, 1};
            GridGenerator::subdivided_hyper_rectangle(
              triangulation, repetitions, p1, p2, true);

            // GridGenerator::cylinder(triangulation, 10., 0.1);

            triangulation.refine_global(2);
          }
        else
          {
            // GridTools::transform(backroty<dim>, triangulation);
            refine_grid();
          }

        // GridTools::transform(roty<dim>, triangulation);
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
} // namespace elasticity


int
main(int argc, char *argv[])
{
  try
    {
      using namespace dealii;
      using namespace elasticity;

      deallog.depth_console(2);

      Utilities::MPI::MPI_InitFinalize mpi_initialization(
        argc, argv, dealii::numbers::invalid_unsigned_int);

      ElaProblem<3> Ela_problem_3d;
      Ela_problem_3d.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  return 0;
}