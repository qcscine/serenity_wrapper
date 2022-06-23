/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SERENITY_CCCALCULATOR_H_
#define SERENITY_CCCALCULATOR_H_

#include "Serenity/Calculators/CalculatorBase.h"
/* Serenity Includes */
#include <settings/Options.h>
/* Scine Includes */
#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/Settings.h>
#include <Utils/Technical/CloneInterface.h>
#include <string>

namespace Serenity {
class Geometry;
class SystemController;
} // namespace Serenity
namespace Sty = Serenity;

namespace Scine {
namespace Utils {}
namespace Serenity {

/**
 * @brief An implementation of the Scine::Core::Calculator for single system CC calculations.
 */
class CCCalculator : public Scine::Utils::CloneInterface<CCCalculator, CalculatorBase> {
 public:
  static constexpr const char* model = "CC";
  static constexpr const char* program = "Serenity";
  /// @brief Default Constructor
  CCCalculator() = default;
  /// @brief Default Destructor.
  ~CCCalculator() = default;
  /// @brief Copy Constructor.
  CCCalculator(const CCCalculator& other) = default;
  /**
   * @brief Getter for the name of the underlying method.
   * @returns Returns the name of the underlying method.
   */
  std::string name() const final;
  /**
   * @brief
   * @return Scine::Utils::PropertyList
   */
  Scine::Utils::PropertyList possibleProperties() const override;
  /**
   * @brief Check if the method family is supported by this calculator.
   * @param methodFamily The method family as all caps string.
   * @return true  If it is supported.
   * @return false If it is not supported.
   */
  bool supportsMethodFamily(const std::string& methodFamily) const final;

 protected:
  void applyFixedSettings(Sty::Settings& settings) const final;
  template<Sty::Options::SCF_MODES ScfMode>
  void calculateImpl();
  void calculateImplRestricted() final {
    this->calculateImpl<Sty::RESTRICTED>();
  }
  void calculateImplUnrestricted() final {
    this->calculateImpl<Sty::UNRESTRICTED>();
  }
  inline std::vector<std::string> availableSolvationModels() const final {
    return {};
  }
};

} /* namespace Serenity */
} /* namespace Scine */

#endif /* SERENITY_CCCALCULATOR_H_ */
