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


  // RandomNumberUInt::RandomNumberUInt(const unsigned int b,
  //                                    const bool same_on_all_ranks = true)
  //   : a(0)
  //   , b(b)
  //   , same_on_all_ranks(same_on_all_ranks)
  //   , timeSeed((
  //       same_on_all_ranks ?
  //         0.0 :
  //         std::chrono::high_resolution_clock::now().time_since_epoch().count()))
  //   , seed_sequence{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >>
  //   32)}
  // {
  //   rng.seed(seed_sequence);
  // }

  // RandomNumberUInt::RandomNumberUInt(const unsigned int a,
  //                                    const unsigned int b,
  //                                    const bool same_on_all_ranks = true)
  //   : a(a)
  //   , b(b)
  //   , same_on_all_ranks(same_on_all_ranks)
  //   , timeSeed((
  //       same_on_all_ranks ?
  //         0.0 :
  //         std::chrono::high_resolution_clock::now().time_since_epoch().count()))
  //   , seed_sequence{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >>
  //   32)}
  // {
  //   rng.seed(seed_sequence);
  // }

  // void
  // RandomNumberUInt::reinit()
  // {
  //   // re-initialize the random number generator with time-dependent seed
  //   timeSeed =
  //     (same_on_all_ranks ?
  //        0.0 :
  //        std::chrono::high_resolution_clock::now().time_since_epoch().count());
  //   std::seed_seq seed_sequence{uint32_t(timeSeed & 0xffffffff),
  //                               uint32_t(timeSeed >> 32)};
  //   rng.seed(seed_sequence);
  // }

  // double
  // RandomNumberUInt::generate()
  // {
  //   return rng();
  // }

  template class Rotation<3>;
} // namespace MyTools