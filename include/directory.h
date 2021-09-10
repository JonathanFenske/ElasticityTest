#ifndef INCLUDE_OHER_TOOLS_H_
#define INCLUDE_OHER_TOOLS_H_

#include <sys/stat.h>

#include <stdexcept>
#include <string>

/*!
 * @namespace Tools
 *
 * @brief Namespace containing all tools that do not fit other more specific namespaces.
 */
namespace Tools
{
  /*!
   * @brief Creates a directory with a name from a given string
   *
   * @param dir_name
   */
  void
  create_data_directory(const char *dir_name);
} // namespace Tools

#endif /* INCLUDE_OHER_TOOLS_H_ */