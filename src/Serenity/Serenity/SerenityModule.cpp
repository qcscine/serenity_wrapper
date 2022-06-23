/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
/* Internal Includes */
#include "Serenity/SerenityModule.h"
#include "Serenity/Calculators/CCCalculator.h"
#include "Serenity/Calculators/DFTCalculator.h"
#include "Serenity/Calculators/HFCalculator.h"
/* External Includes */
#include <Core/DerivedModule.h>
#include <Core/Exceptions.h>
#include <Utils/CalculatorBasics/StatesHandler.h>
#include <Utils/Settings.h>

namespace Scine {
namespace Serenity {

std::string SerenityModule::name() const noexcept {
  return "Serenity";
}

using InterfaceModelMap =
    boost::mpl::map<boost::mpl::pair<Scine::Core::Calculator, boost::mpl::vector<DFTCalculator, HFCalculator, CCCalculator>>>;

boost::any SerenityModule::get(const std::string& interface, const std::string& model) const {
  boost::any resolved = Scine::Core::DerivedModule::resolve<InterfaceModelMap>(interface, model);

  // Throw an exception if we could not match an interface or model
  if (resolved.empty()) {
    throw Scine::Core::ClassNotImplementedError();
  }

  return resolved;
}

bool SerenityModule::has(const std::string& interface, const std::string& model) const noexcept {
  return Scine::Core::DerivedModule::has<InterfaceModelMap>(interface, model);
}

std::vector<std::string> SerenityModule::announceInterfaces() const noexcept {
  return Scine::Core::DerivedModule::announceInterfaces<InterfaceModelMap>();
}

std::vector<std::string> SerenityModule::announceModels(const std::string& interface) const noexcept {
  return Scine::Core::DerivedModule::announceModels<InterfaceModelMap>(interface);
}

std::shared_ptr<Scine::Core::Module> SerenityModule::make() {
  return std::make_shared<SerenityModule>();
}

std::vector<std::shared_ptr<Scine::Core::Module>> moduleFactory() {
  return {Scine::Serenity::SerenityModule::make()};
}

} /* namespace Serenity */
} /* namespace Scine */

#ifdef __MINGW32__
/* MinGW builds are problematic. We build with default visibility, and adding
 * an attribute __dllexport__ specifically for this singular symbol leads to the
 * loss of all other weak symbols. Essentially, here we have just expanded the
 * BOOST_DLL_ALIAS macro in order to declare the type-erased const void*
 * 'moduleFactory' without any symbol visibility attribute additions that could
 * confuse the MinGW linker, which per Boost DLL documentation is unable to mix
 * weak attributes and __dllexport__ correctly.
 *
 * If ever the default visibility for this translation unit is changed, we
 * will have to revisit this bit of code for the MinGW platform again.
 *
 * Additionally, more recent Boost releases may have fixed this problem.
 * See the macro BOOST_DLL_FORCE_ALIAS_INSTANTIATIONS as used in the library's
 * example files.
 */
extern "C" {
const void* moduleFactory =
    reinterpret_cast<const void*>(reinterpret_cast<intptr_t>(&Serenity::ScineInterface::moduleFactory));
}
#else
BOOST_DLL_ALIAS(Scine::Serenity::moduleFactory, moduleFactory)
#endif
