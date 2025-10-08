/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SERENITY_SERENITYEMBEDDINGCALCULATOR_H
#define SERENITY_SERENITYEMBEDDINGCALCULATOR_H

#include "Serenity/Calculators/SerenityEmbeddingCalculatorSettings.h"
/* Serenity Includes */
#include <postHF/LocalCorrelation/LocalCorrelationController.h>
#include <settings/EmbeddingSettings.h>
#include <settings/Options.h>
#include <tasks/BasisSetTruncationTask.h>
#include <tasks/LocalizationTask.h>
/* Scine Includes */
#include "Serenity/Calculators/CalculatorBase.h"
#include <Core/Interfaces/Calculator.h>
#include <Core/Interfaces/EmbeddingCalculator.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/Settings.h>
#include <Utils/Technical/CloneInterface.h>
#include <memory>
#include <string>
#include <vector>

namespace Sty = Serenity;
namespace Scine {
namespace Serenity {

/**
 * @class SerenityEmbeddingCalculator SerenityEmbeddingCalculator.h
 * @brief QM/QM embedding calculator for serenity.
 */
class SerenityEmbeddingCalculator
  : public Scine::Utils::CloneInterface<SerenityEmbeddingCalculator, Scine::Core::EmbeddingCalculator, Scine::Core::Calculator> {
 public:
  static constexpr const char* model = "QMQM";
  static constexpr const char* program = "Serenity";
  /**
   * @brief Constructor.
   */
  SerenityEmbeddingCalculator();
  /**
   * @brief Default destructor.
   */
  ~SerenityEmbeddingCalculator() = default;
  /**
   * @brief Copy constructor.
   * @param other The calculator to be copied.
   */
  SerenityEmbeddingCalculator(const SerenityEmbeddingCalculator& other);
  /**
   * @brief Sets the underlying calculators.
   */
  void setUnderlyingCalculators(std::vector<std::shared_ptr<Core::Calculator>> underlyingCalculators) override;
  /**
   * @brief Getter for the name of the underlying method.
   * @returns Returns the name of the underlying method.
   */
  std::string name() const final;
  /**
   * @brief Check if the method family is supported by this calculator.
   * @param methodFamily The method family as all caps string.
   * @return Returns true if it is supported and false otherwise.
   */
  bool supportsMethodFamily(const std::string& methodFamily) const final;
  /**
   * @brief Gets the underlying calculators.
   */
  std::vector<std::shared_ptr<Core::Calculator>> getUnderlyingCalculators() const override;
  /**
   * @brief Changes the molecular structure to calculate.
   * @param structure A new Utils::AtomCollection to save.
   */
  void setStructure(const Utils::AtomCollection& structure) override;
  /**
   * @brief Gets the molecular structure as a const Utils::AtomCollection&.
   * @return a const Utils::AtomCollection&.
   */
  std::unique_ptr<Utils::AtomCollection> getStructure() const override;
  /**
   * @brief Allows to modify the positions of the underlying Utils::AtomCollection.
   * @param newPositions The new positions to be assigned to the underlying Utils::AtomCollection.
   */
  void modifyPositions(Utils::PositionCollection newPositions) override;
  /**
   * @brief Getter for the coordinates of the underlying Utils::AtomCollection.
   */
  const Utils::PositionCollection& getPositions() const override;
  /**
   * @brief Sets the properties to calculate.
   * @param requiredProperties A Utils::PropertyList, a sequence of bits that represent the
   *        properties that must be calculated.
   */
  void setRequiredProperties(const Utils::PropertyList& requiredProperties) override;
  /**
   * @brief Getter for the properties to calculate.
   */
  Utils::PropertyList getRequiredProperties() const override;
  /**
   * @brief Returns the list of the possible properties to calculate.
   */
  Utils::PropertyList possibleProperties() const override;
  /**
   * @brief The main function running calculations.
   *
   * @param description   The calculation description.
   * @return Utils::Result Return the result of the calculation. The object contains the
   *                       properties that were given as requirement by the
   *                       Calculator::setRequiredProperties function.
   */
  const Utils::Results& calculate(std::string description) override;
  /**
   * @brief Accessor for the settings.
   * @return Utils::Settings& The settings.
   */
  Utils::Settings& settings() override;
  /**
   * @brief Constant accessor for the settings.
   * @return const Utils::Settings& The settings.
   */
  const Utils::Settings& settings() const override;
  /**
   * @brief Accessor for the saved instance of Utils::Results.
   * @return Utils::Results& The results of the previous calculation.
   */
  Utils::Results& results() override;
  /**
   * @brief Constant accessor for the Utils::Results.
   * @return const Utils::Results& The results of the previous calculation.
   */
  const Utils::Results& results() const override;
  /**
   * @brief Implements Core::StateHandableObject::getState().
   * @return std::shared_ptr<Core::State> The current state
   */
  std::shared_ptr<Core::State> getState() const final;
  /**
   * @brief Implements Core::StateHandableObject::loadState().
   * @param state The new state.
   */
  void loadState(std::shared_ptr<Core::State> /*state*/) final;
  /**
   * @brief Whether the calculator has no underlying Python code and can therefore
   * release the global interpreter lock in Python bindings
   */
  bool allowsPythonGILRelease() const override {
    const auto underlyingCalculators = getUnderlyingCalculators();
    return std::all_of(underlyingCalculators.begin(), underlyingCalculators.end(),
                       [](const auto& c) { return c->allowsPythonGILRelease(); });
  };
  /**
   * @brief Does nothing since all Serenity BaseCalculators have the same settings.
   */
  void addUnderlyingSettings() override;

  /**
   * @brief Calculate the core Hamiltonian and prepare the orbital coefficient,
   *        and energy objects for a CAS calculation using the embedding.
   */
  template<Sty::Options::SCF_MODES SCFMode>
  void prepareCAS();

 protected:
  std::unique_ptr<Utils::Settings> _settings;
  std::vector<std::shared_ptr<CalculatorBase>> _underlyingCalculators;
  std::unique_ptr<Scine::Utils::Results> _results;
  Scine::Utils::PropertyList _requiredProperties;
  std::shared_ptr<Sty::Geometry> _geometry;
  std::unique_ptr<Scine::Utils::PositionCollection> _scinePositions;
  std::unique_ptr<double> _finalEnergy;
  std::vector<std::shared_ptr<Eigen::MatrixXd>> _oneParticleIntegralsAlpha;
  std::vector<std::shared_ptr<Eigen::MatrixXd>> _oneParticleIntegralsBeta;
  std::vector<std::shared_ptr<Eigen::MatrixXd>> _coefficientMatrixAlpha;
  std::vector<std::shared_ptr<Eigen::MatrixXd>> _coefficientMatrixBeta;

 private:
  template<Sty::Options::SCF_MODES SCFMode>
  void runEmbeddingTask();
  static bool isLocalCorrelation(std::shared_ptr<CalculatorBase> calculatorBase);
  void updateEmbeddingSettings(Sty::EmbeddingSettings& embeddingSettings);
  void updateLocalizationSettings(Sty::LocalizationTaskSettings& localizationTaskSettings);
  static void updateLocalCorrelationSettings(Sty::LocalCorrelationSettings& localCorrelationSettings,
                                             const Utils::Settings& settings, bool isLocalCorrelationSystem = false);
  void updateBasisSetTruncationSettings(Sty::BasisSetTruncationTaskSettings& basisSetTruncationTaskSettings);
  void setIOOptions();
  template<Sty::Options::SCF_MODES SCFMode>
  double getCorrelationEnergy(std::shared_ptr<Sty::SystemController> activeSystemController,
                              const Sty::Options::PNO_METHOD& pnoMethod);
  template<Sty::Options::SCF_MODES SCFMode>
  double getReferenceEnergy(std::shared_ptr<Sty::SystemController> activeSystemController);
  void applySettingsToUnderlyingCalculators();
  void adjustCharges();
  bool subsystemsAreGeneratedFromIndices() const;
  void subsystemAtomAssignmentSanityCheck();
  void assignMethodToSubsystems();
  template<Sty::Options::SCF_MODES SCFMode>
  void runStaticEmbeddingForMultipleSubsystems();
  template<Sty::Options::SCF_MODES SCFMode>
  void runTopDownEmbeddingWithSingleSystemLocalCorrelation();
  template<Sty::Options::SCF_MODES SCFMode>
  void runTopDownEmbeddingOnlyMeanField();
  template<Sty::Options::SCF_MODES SCFMode>
  static Eigen::MatrixXd compactIntegrals(const Sty::SPMatrix<SCFMode>& spMatrix);
  template<Sty::Options::SCF_MODES SCFMode>
  static Sty::SPMatrix<SCFMode>
  extractValenceOrbitals(const Sty::SPMatrix<SCFMode>& coefficients,
                         const Sty::SpinPolarizedData<SCFMode, std::vector<unsigned int>>& valenceRange);
  template<Sty::Options::SCF_MODES SCFMode>
  static Sty::SpinPolarizedData<SCFMode, std::vector<unsigned int>>
  joinOrbitalRanges(const Sty::SpinPolarizedData<SCFMode, std::vector<unsigned int>>& rangeOne,
                    const Sty::SpinPolarizedData<SCFMode, std::vector<unsigned int>>& rangeTwo);
  template<Sty::Options::SCF_MODES SCFMode>
  std::shared_ptr<Sty::SystemController> loadSupersystemIfGiven();
};

} /* namespace Serenity */
} /* namespace Scine */

#endif // SERENITY_SERENITYEMBEDDINGCALCULATOR_H
