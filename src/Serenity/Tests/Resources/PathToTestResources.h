/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SERENITY_PATHTOTESTRESOURCES_H
#define SERENITY_PATHTOTESTRESOURCES_H

#include <boost/dll.hpp>
#include <boost/filesystem.hpp>
#include <string>

namespace Scine {
namespace Serenity {
namespace Tests {

std::string pathToTestResources() {
  boost::filesystem::path pathToProgram = boost::dll::program_location().parent_path();
  return pathToProgram.string() + "/Resources";
}

} /* namespace Tests */
} /* namespace Serenity */
} /* namespace Scine */

#endif // SERENITY_PATHTOTESTRESOURCES_H
