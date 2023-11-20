/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SERENITY_SERENITYSTATE_H_
#define SERENITY_SERENITYSTATE_H_

/* Serenity Includes*/
#include <io/Filesystem.h>
#include <settings/Settings.h>
#include <system/SystemController.h>
/* Scine Includes */
#include <Core/BaseClasses/StateHandableObject.h>

namespace Sty = Serenity;
namespace Scine {
namespace Serenity {

/**
 * @brief The definition of a loadable serenity state.
 */
class SerenityState : public Scine::Core::State {
 public:
  SerenityState(std::shared_ptr<Sty::SystemController> s) : system(s){};
  ~SerenityState() {
    // TODO
    //      remove_all(system->getSettings().path);
  }
  std::shared_ptr<Sty::SystemController> system;
};

} /* namespace Serenity */
} /* namespace Scine */

#endif /* SERENITY_SERENITYSTATE_H_ */
