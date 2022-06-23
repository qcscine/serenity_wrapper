/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SERENITY_CALCULATORBASE_H_
#define SERENITY_CALCULATORBASE_H_

/* Serenity Includes */
#include "data/matrices/DensityMatrix.h"
#include "settings/Options.h"
/* Scine Includes */
#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/Settings.h>
#include <Utils/Technical/CloneInterface.h>
#include <string>

namespace Serenity {
class Geometry;
class SystemController;
template<Options::SCF_MODES ScfMode, class T, typename E>
class SpinPolarizedData;
class Settings;
} // namespace Serenity
// Namespace alias to avoid abiguity
namespace Sty = Serenity;

namespace Scine {
namespace Utils {}
namespace Serenity {

class ScineSettings;

/**
 * @brief A base class for all Scine::Calculators implemented in Serenity.
 *
 * This class comprises some of the common mapping features required in order to couple
 * Serenity and the Scine interfaces.
 */
class CalculatorBase : public Scine::Utils::CloneInterface<Scine::Utils::Abstract<CalculatorBase>, Scine::Core::Calculator> {
 public:
  /// @brief Default Constructor
  CalculatorBase();
  /// @brief Default Destructor.
  ~CalculatorBase() override;
  /// @brief Copy Constructor.
  CalculatorBase(const CalculatorBase& other);
  /**
   * @brief Sets new structure and initializes the underlying method with the parameter given in the settings.
   * @param structure The structure to be assigned.
   */
  void setStructure(const Scine::Utils::AtomCollection& structure) final;
  /**
   * @brief Getter for the underlying element types and positions.
   */
  std::unique_ptr<Scine::Utils::AtomCollection> getStructure() const final;
  /**
   * @brief Allows to modify the positions of the underlying Utils::AtomCollection
   * @param newPositions the new positions to be assigned to the underlying Utils::AtomCollection
   */
  void modifyPositions(Scine::Utils::PositionCollection newPositions) final;
  /**
   * @brief Getter for the coordinates of the underlying Utils::AtomCollection
   */
  const Scine::Utils::PositionCollection& getPositions() const final;
  /**
   * @brief Sets the properties to calculate.
   * @param requiredProperties a Scine::Utils::PropertyList object, contains an enum class that work as
   *        a bitset, switching on and off the bits representing a property.
   */
  void setRequiredProperties(const Scine::Utils::PropertyList& requiredProperties) final;
  /**
   * @brief Gets the properties to calculate.
   * @return requiredProperties a Scine::Utils::PropertyList object, contains an enum class that work as
   *        a bitset, switching on and off the bits representing a property.
   */
  virtual Scine::Utils::PropertyList getRequiredProperties() const final;
  /**
   * @brief Getter for the name of the underlying method.
   * @returns Returns the name of the underlying method.
   */
  virtual std::string name() const = 0;
  /**
   * @brief Returns a list of all properties this Calculator can produce.
   * @return Scine::Utils::PropertyList A list of all properties this Calculator can produce.
   */
  virtual Scine::Utils::PropertyList possibleProperties() const = 0;
  /**
   * @brief The main function running calculations.
   * @param dummy   A dummy parameter.
   * @return Scine::Utils::Results Return the result of the calculation.
   */
  const Scine::Utils::Results& calculate(std::string dummy) final;
  /**
   * @brief Accessor for the Settings used in this method wrapper.
   * @returns Scine::Utils::Settings& The Settings.
   */
  Scine::Utils::Settings& settings() final;
  /**
   * @brief Const accessor for the Settings used in this method wrapper.
   * @returns const Scine::Utils::Settings& The Settings.
   */
  const Scine::Utils::Settings& settings() const final;
  /**
   * @brief Accessor for the Results stored in this method wrapper.
   * @returns Scine::Utils::Results& The results of the previous calculation.
   */
  Scine::Utils::Results& results() final;
  /**
   * @brief Const accessor for the Results used in this method wrapper.
   * @returns const Scine::Utils::Results& The results of the previous calculation.
   */
  const Scine::Utils::Results& results() const final;
  /**
   * @brief Exchange the current state/system for a different one.
   * @param state The new state/system.
   */
  void loadState(std::shared_ptr<Scine::Core::State> state) final;
  /**
   * @brief Get a copy of current state/system.
   * @return std::shared_ptr<Scine::Core::State> The current state/system.
   */
  std::shared_ptr<Scine::Core::State> getState() const final;

 protected:
  std::unique_ptr<ScineSettings> _settings;
  std::unique_ptr<Scine::Utils::Results> _results;
  Scine::Utils::PropertyList _requiredProperties;
  std::shared_ptr<Sty::SystemController> _system;
  std::shared_ptr<Sty::Geometry> _geometry;
  std::unique_ptr<Scine::Utils::PositionCollection> _scinePositions;
  bool _moved;

  /**
   * @brief Apply all settings required to be a fixed value as determined by the Calculator type.
   * @param settings The Serenity::Settings to be modified.
   */
  virtual void applyFixedSettings(Sty::Settings& settings) const = 0;
  /**
   * @brief The actual implementation of the calculation(s).
   *
   * Any implementation is expected to use the protected members in this base class
   * to populate the results of the/this calculation, which are also a member variable.
   *
   * These functions could have been pure virtual templates if C++ would have allowed for it.
   */
  virtual void calculateImplRestricted() = 0;
  /**
   * @brief The actual implementation of the calculation(s).
   *
   * Any implementation is expected to use the protected members in this base class
   * to populate the results of the/this calculation, which are also a member variable.
   *
   * These functions could have been pure virtual templates if C++ would have allowed for it.
   */
  virtual void calculateImplUnrestricted() = 0;

  /**
   * @brief The available solvation models for each implementation
   */
  virtual std::vector<std::string> availableSolvationModels() const = 0;

  template<Sty::Options::SCF_MODES ScfMode>
  Scine::Utils::DensityMatrix convertDensityMatrix(Sty::DensityMatrix<ScfMode> dmat,
                                                   Sty::SpinPolarizedData<ScfMode, unsigned int, void> nEl) const;
  template<Sty::Options::SCF_MODES ScfMode>
  std::vector<double> getMullikenCharges() const;
  template<Sty::Options::SCF_MODES ScfMode>
  std::vector<double> getHirshfeldCharges() const;

 private:
  template<Sty::Options::SCF_MODES ScfMode>
  std::vector<double> populationToCharges(const Sty::SpinPolarizedData<ScfMode, Eigen::VectorXd>& populations) const;
  inline static double
  electronPopulationAtAtom(const Sty::SpinPolarizedData<Sty::Options::SCF_MODES::RESTRICTED, Eigen::VectorXd>& populations,
                           int i) {
    return populations[i];
  };
  inline static double
  electronPopulationAtAtom(const Sty::SpinPolarizedData<Sty::Options::SCF_MODES::UNRESTRICTED, Eigen::VectorXd>& populations,
                           int i) {
    return populations.alpha[i] + populations.beta[i];
  };
};

} /* namespace Serenity */
} /* namespace Scine */

#endif /* SERENITY_CALCULATORBASE_H_ */
