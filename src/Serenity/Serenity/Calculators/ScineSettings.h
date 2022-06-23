/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SERENITY_SCINESETTINGS_H_
#define SERENITY_SCINESETTINGS_H_

/* Scine Includes */
#include <Utils/Settings.h>

namespace Serenity {
class Settings;
}
namespace Sty = Serenity;

namespace Scine {
namespace Serenity {

/**
 * @brief A SCINE style settings object for the settings relevanet to DFT in Serenity.
 */
class ScineSettings : public Scine::Utils::Settings {
 public:
  /**
   * @brief Construct a new ScineSettings object.
   */
  ScineSettings();
  /**
   * @brief Applies these SCINE style settings to the Serenity style settings.
   * @param settings The Serenity style settings.
   */
  void applyTo(Sty::Settings& settings);
  /**
   * @brief If spin mode is set to 'any', this selects a fitting serenity spin mode depending on the spin multiplicity
   */
  void resolveSpinMode();
};

} /* namespace Serenity */
} /* namespace Scine */

#endif /* SERENITY_SCINESETTINGS_H_ */
