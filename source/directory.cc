#include <directory.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <iostream>

namespace Tools
{
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

} // namespace Tools