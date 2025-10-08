/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SERENITY_SERENITYEMBEDDINGCALCULATORSETTINGS_H
#define SERENITY_SERENITYEMBEDDINGCALCULATORSETTINGS_H

#include "Serenity/Calculators/ScineSettings.h"

namespace Sty = Serenity;

namespace Scine {
namespace Serenity {

/**
 * @class SerenityEmbeddingCalculatorSettings SerenityEmbeddingCalculatorSettings.h
 * @brief The settings for the Serenity embedding calculator.
 */
class SerenityEmbeddingCalculatorSettings : public ScineSettings {
 public:
  /**
   * @brief Constructor. Adds all setting members.
   */
  SerenityEmbeddingCalculatorSettings();
  ~SerenityEmbeddingCalculatorSettings() = default;
};
} /* namespace Serenity */
} /* namespace Scine */
#endif // SERENITY_SERENITYEMBEDDINGCALCULATORSETTINGS_H
