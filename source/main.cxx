#include "ela_std.h"
#include "ela_ms.h"

int main (int argc, char *argv[])
{
  try
    {
      using namespace dealii;
      using namespace Elasticity;

      deallog.depth_console(2);

      // For debugging or if we use PetSc it may be nice to limit the threads on
      // each process.
#if DEBUG
      dealii::MultithreadInfo::set_thread_limit(1);
      dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(
        argc, argv, /* max_threads */ 1);

#else
      dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(
        argc, argv, dealii::numbers::invalid_unsigned_int);
#endif

      const bool say_hello_from_cluster =
#if DEBUG
        true;
#else
        false;
#endif

      if (say_hello_from_cluster)
        {
          char processor_name[MPI_MAX_PROCESSOR_NAME];
          int  name_len;
          MPI_Get_processor_name(processor_name, &name_len);

          std::string proc_name(processor_name, name_len);

          std::cout << "Hello from   " << proc_name << "   Rank:   "
                    << dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
                    << "   out of   "
                    << dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)
                    << "   | cores = " << dealii::MultithreadInfo::n_cores()
                    << "   | threads = " << dealii::MultithreadInfo::n_threads()
                    << std::endl;
        }
      const bool direct_solver = true;
      const bool neumann_bc = false;

      ElaStd<3> ela_std(direct_solver, neumann_bc);
      ela_std.run();

      ElaMs<3> ela_ms(direct_solver, neumann_bc);
      ela_ms.run();
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