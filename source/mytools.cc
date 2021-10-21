#include "mytools.h"

#include <sys/stat.h>
#include <sys/types.h>

#include <iostream>

#include "mytools.tpp"

namespace MyTools
{
  template void
  set_dirichlet_id(const Point<2>                          &p1,
                   const Point<2>                          &p2,
                   const unsigned int                       face_id,
                   const unsigned int                       id,
                   parallel::distributed::Triangulation<2> &triangulation);

  template void
  set_dirichlet_id(const Point<3>                          &p1,
                   const Point<3>                          &p2,
                   const unsigned int                       face_id,
                   const unsigned int                       id,
                   parallel::distributed::Triangulation<3> &triangulation);

  template const std::vector<unsigned int>
  get_repetitions(const Point<2> &p1, const Point<2> &p2);

  template const std::vector<unsigned int>
  get_repetitions(const Point<3> &p1, const Point<3> &p2);

  void
  create_data_directory(const char *dir_name)
  {
    struct stat info;

    if (stat(dir_name, &info) == -1)
      {
        const int dir_err =
          mkdir(dir_name, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (-1 == dir_err)
          {
            throw std::runtime_error(
              "Error creating directory! It might already exist or you do not have write permissions in this folder.");
          }
      }
  }

  template class Rotation<3>;
} // namespace MyTools