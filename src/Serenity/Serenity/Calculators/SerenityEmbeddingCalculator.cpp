/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Serenity/Calculators/SerenityEmbeddingCalculator.h"
#include "Serenity/Calculators/CCCalculator.h"
#include "Serenity/Calculators/SerenityEmbeddingCalculatorSettings.h"
#include "Serenity/Utilities/SerenityConversionFunctions.h"
/* Serenity includes */
#include "Utils/CalculatorBasics/CalculationRoutines.h"
#include <data/ElectronicStructure.h>
#include <data/OrbitalController.h>
#include <energies/EnergyComponentController.h>
#include <energies/EnergyContributions.h>
#include <geometry/Geometry.h>
#include <io/FCIDumpFileWriter.h>
#include <io/FormattedOutput.h>
#include <io/FormattedOutputStream.h>
#include <io/IOOptions.h>
#include <settings/Settings.h>
#include <system/SystemController.h>
#include <tasks/DFTEmbeddedLocalCorrelationTask.h>
#include <tasks/FCIDumpFileWriterTask.h>
#include <tasks/OrbitalsIOTask.h>
#include <tasks/SystemAdditionTask.h>
#include <tasks/TDEmbeddingTask.h>
#include <tasks/TopDownStaticEmbeddingTask.h>
#include <algorithm>
#include <boost/algorithm/string.hpp>

namespace Scine {
namespace Serenity {

SerenityEmbeddingCalculator::SerenityEmbeddingCalculator() {
  this->_settings = std::make_unique<SerenityEmbeddingCalculatorSettings>();
  this->_results = std::make_unique<Scine::Utils::Results>();
}

SerenityEmbeddingCalculator::SerenityEmbeddingCalculator(const SerenityEmbeddingCalculator& rhs)
  : CloneInterface(rhs), _requiredProperties(rhs._requiredProperties) {
  this->setLog(rhs.getLog());
  setUnderlyingCalculators(rhs.getUnderlyingCalculators());
  auto valueCollection = dynamic_cast<const Utils::UniversalSettings::ValueCollection&>(rhs.settings());
  this->_settings =
      std::make_unique<Utils::Settings>(Utils::Settings(valueCollection, rhs.settings().getDescriptorCollection()));
  _results = std::make_unique<Scine::Utils::Results>(*rhs._results);
  _geometry = rhs._geometry;
}
void SerenityEmbeddingCalculator::setUnderlyingCalculators(std::vector<std::shared_ptr<Core::Calculator>> underlyingCalculators) {
  _underlyingCalculators = {};
  for (const auto& calculator : underlyingCalculators) {
    auto casted = std::dynamic_pointer_cast<CalculatorBase>(calculator);
    if (casted == nullptr) {
      throw std::runtime_error("The serenity embedding calculator works only combined with Serenity calculators.");
    }
    _underlyingCalculators.push_back(casted);
  }
  if (_underlyingCalculators.size() < 2) {
    throw std::runtime_error("Embedding calculations are only possible for at least two systems.");
  }
}
std::string SerenityEmbeddingCalculator::name() const {
  return "QMQM";
}
bool SerenityEmbeddingCalculator::supportsMethodFamily(const std::string& methodFamily) const {
  return methodFamily == this->name();
}
std::vector<std::shared_ptr<Core::Calculator>> SerenityEmbeddingCalculator::getUnderlyingCalculators() const {
  std::vector<std::shared_ptr<Core::Calculator>> toReturn;
  for (const auto& calculator : _underlyingCalculators) {
    auto casted = std::dynamic_pointer_cast<Core::Calculator>(calculator);
    if (casted == nullptr) {
      throw std::runtime_error("Calculator base class cast failed! This is an implementation error.");
    }
    toReturn.push_back(casted);
  }
  return toReturn;
}
void SerenityEmbeddingCalculator::setStructure(const Utils::AtomCollection& structure) {
  if (!this->subsystemsAreGeneratedFromIndices()) {
    throw std::runtime_error("The structure cannot be directly manipulated for the QM/QM embedding in serenity\n"
                             "if no atom indices for the subsystems are provided. It must be set through the\n"
                             "underlying calculators.");
  }
  _geometry = std::make_shared<Sty::Geometry>(SerenityConversionFunctions::atomCollectionToGeometry(structure));
  _scinePositions = std::make_unique<Scine::Utils::PositionCollection>(structure.getPositions());
  _results = std::make_unique<Scine::Utils::Results>();

  // Update structures in the underlying calculators.
  this->subsystemAtomAssignmentSanityCheck();
  const auto& assignments = this->_settings->getIntListList("qmqm_atom_indices");
  for (unsigned int iCalculator = 0; iCalculator < _underlyingCalculators.size(); ++iCalculator) {
    std::vector<unsigned int> assignment;
    // This odd insert call is necessary to cast the ints to unsigned ints.
    assignment.insert(assignment.begin(), assignments[iCalculator].begin(), assignments[iCalculator].end());
    Utils::AtomCollection atomCollectionCopy(*this->getStructure());
    atomCollectionCopy.keepAtomsByIndices(assignment);
    _underlyingCalculators[iCalculator]->setStructure(atomCollectionCopy);
  }
}
std::unique_ptr<Utils::AtomCollection> SerenityEmbeddingCalculator::getStructure() const {
  if (this->subsystemsAreGeneratedFromIndices()) {
    if (!_geometry || !_scinePositions) {
      throw std::runtime_error("Missing geometry in a Serenity Calculator");
    }
    return std::make_unique<Scine::Utils::AtomCollection>(SerenityConversionFunctions::geometryToAtomCollection(*_geometry));
  }
  auto supersystem = std::make_unique<Utils::AtomCollection>();
  for (const auto& calculator : _underlyingCalculators) {
    *supersystem += *calculator->getStructure();
  }
  return supersystem;
}
void SerenityEmbeddingCalculator::modifyPositions(Utils::PositionCollection newPositions) {
  (void)newPositions;
  throw std::runtime_error("The structure cannot be directly manipulated for the QM/QM embedding in serenity."
                           " It must be set through the underlying calculators.");
}
const Utils::PositionCollection& SerenityEmbeddingCalculator::getPositions() const {
  throw std::runtime_error("The structure is not directly available from the Serenity QM/QM embedding calculator.");
}
void SerenityEmbeddingCalculator::setRequiredProperties(const Utils::PropertyList& requiredProperties) {
  for (auto& calculator : _underlyingCalculators) {
    calculator->setRequiredProperties(requiredProperties);
  }
  _requiredProperties = requiredProperties;
}
Utils::PropertyList SerenityEmbeddingCalculator::getRequiredProperties() const {
  return _requiredProperties;
}
Utils::PropertyList SerenityEmbeddingCalculator::possibleProperties() const {
  Utils::PropertyList properties{Utils::Property::Energy | Utils::Property::AtomicCharges |
                                 Utils::Property::OneElectronMatrix | Utils::Property::CoefficientMatrix |
                                 Utils::Property::PartialEnergies | Utils::Property::AOtoAtomMapping};
  return properties;
}
const Utils::Results& SerenityEmbeddingCalculator::calculate(std::string description) {
  this->setIOOptions();
  this->applySettingsToUnderlyingCalculators();
  const Sty::Options::SCF_MODES scfMode = _underlyingCalculators[0]->getSystemController()->getSCFMode();
  try {
    if (scfMode == Sty::RESTRICTED) {
      this->runEmbeddingTask<Sty::RESTRICTED>();
      this->prepareCAS<Sty::RESTRICTED>();
    }
    else {
      this->runEmbeddingTask<Sty::UNRESTRICTED>();
      this->prepareCAS<Sty::UNRESTRICTED>();
    }
  }
  catch (Sty::SerenityError& e) {
    throw Core::UnsuccessfulCalculationException(e.what());
  }
  if (!_finalEnergy) {
    throw Core::UnsuccessfulCalculationException("No energies available after QM/QM calculation with Serenity.");
  }
  _results = std::make_unique<Scine::Utils::Results>();
  _results->set<Utils::Property::Description>(std::move(description));
  _results->set<Scine::Utils::Property::ProgramName>("serenity");
  _results->set<Scine::Utils::Property::Energy>(*_finalEnergy);
  _results->set<Scine::Utils::Property::SuccessfulCalculation>(true);
  if (_requiredProperties.containsSubSet(Utils::Property::AOtoAtomMapping)) {
  }

  return *_results;
}
bool SerenityEmbeddingCalculator::isLocalCorrelation(std::shared_ptr<CalculatorBase> calculatorBase) {
  std::string method = calculatorBase->settings().getString("method");
  boost::algorithm::to_lower(method);
  if (method.find("dlpno") != std::string::npos) {
    return true;
  }
  if (calculatorBase->name() == CCCalculator::calculatorName) {
    throw std::runtime_error(
        "Canonical coupled cluster is not supported with QM/QM embedding through the Serenity wrapper.");
  }
  return false;
}
void SerenityEmbeddingCalculator::updateEmbeddingSettings(Sty::EmbeddingSettings& embeddingSettings) {
  auto value = this->_settings->getString("embedding_method");
  Sty::Options::resolve(value, embeddingSettings.embeddingMode);

  auto methodInput = Scine::Utils::CalculationRoutines::splitIntoMethodAndDispersion(
      this->_settings->getString("non_additive_xc_functional"));
  Sty::Options::resolve(methodInput.first, embeddingSettings.naddXCFunc);
  if (!methodInput.second.empty()) {
    Sty::Options::resolve(methodInput.second, embeddingSettings.dispersion);
  }
  embeddingSettings.fermiShift = this->_settings->getDouble("fermi_shift");
}
void SerenityEmbeddingCalculator::updateLocalizationSettings(Sty::LocalizationTaskSettings& localizationTaskSettings) {
  auto value = this->_settings->getString("loc_type");
  Sty::Options::resolve(value, localizationTaskSettings.locType);
  localizationTaskSettings.splitValenceAndCore = true;
  // If we want to save CAS integrals, we must also localize the virtual orbitals.
  localizationTaskSettings.localizeVirtuals = !this->_settings->getIntList("cas_systems").empty();
  localizationTaskSettings.useEnergyCutOff = false;
}
void SerenityEmbeddingCalculator::updateLocalCorrelationSettings(Sty::LocalCorrelationSettings& localCorrelationSettings,
                                                                 const Utils::Settings& settings,
                                                                 bool isLocalCorrelationSystem) {
  auto value = settings.getString("pno_settings");
  Sty::Options::resolve(value, localCorrelationSettings.pnoSettings);
  localCorrelationSettings.useFrozenCore = settings.getBool("use_frozen_core");
  std::string method = (isLocalCorrelationSystem) ? settings.getString("method") : "DLPNO-MP2";
  Sty::Options::resolve(method, localCorrelationSettings.method);
  localCorrelationSettings.ignoreMemoryConstraints = true;
}
void SerenityEmbeddingCalculator::updateBasisSetTruncationSettings(Sty::BasisSetTruncationTaskSettings& basisSetTruncationTaskSettings) {
  auto value = this->_settings->getString("trunc_algorithm");
  Sty::Options::resolve(value, basisSetTruncationTaskSettings.truncAlgorithm);
  basisSetTruncationTaskSettings.netThreshold = this->_settings->getDouble("net_threshold");
}
void SerenityEmbeddingCalculator::setIOOptions() {
  const bool printOutput = this->_settings->getBool("show_serenity_output");
  if (!printOutput) {
    Sty::GLOBAL_PRINT_LEVEL = Sty::Options::GLOBAL_PRINT_LEVELS::MINIMUM;
    Sty::iOOptions.printFinalOrbitalEnergies = false;
    Sty::iOOptions.printGeometry = false;
    Sty::iOOptions.printSCFCycleInfo = false;
    Sty::iOOptions.printSCFResults = false;
    Sty::iOOptions.printDebugInfos = false;
    Sty::iOOptions.printGridInfo = false;
    Sty::iOOptions.gridAccuracyCheck = false;
    Sty::iOOptions.timingsPrintLevel = 0;
  }
}
Utils::Settings& SerenityEmbeddingCalculator::settings() {
  return *this->_settings;
}
const Utils::Settings& SerenityEmbeddingCalculator::settings() const {
  return *this->_settings;
}
Utils::Results& SerenityEmbeddingCalculator::results() {
  return *this->_results;
}
const Utils::Results& SerenityEmbeddingCalculator::results() const {
  return *this->_results;
}
std::shared_ptr<Core::State> SerenityEmbeddingCalculator::getState() const {
  throw std::runtime_error("State handling is not supported for QM/QM embedding through the Serenity wrapper.");
}
void SerenityEmbeddingCalculator::loadState(std::shared_ptr<Core::State>) {
  throw std::runtime_error("State loading is not supported for QM/QM embedding through the Serenity wrapper.");
}
void SerenityEmbeddingCalculator::addUnderlyingSettings() {
  /*
   * Nothing to be done here. All Serenity base calculators have the same setting objects.
   */
}
void SerenityEmbeddingCalculator::applySettingsToUnderlyingCalculators() {
  // If we consider only top-down embedding, the basis set must also be forwarded.
  std::vector<std::string> settingsNotForwarded = {
      Utils::SettingsNames::spinMultiplicity,
      Utils::SettingsNames::molecularCharge,
      Utils::SettingsNames::method,
      "basis_auxJLabel",
      "basis_auxCLabel",
      "basis_makeSphericalBasis",
      "pno_settings",
      "use_frozen_core",
      "grid_accuracy",
      "grid_smallGridAccuracy",
      "grid_gridType",
      "basis_basisLibPath",
  };
  for (auto& calculator : this->_underlyingCalculators) {
    for (const auto& descriptor : calculator->settings().getDescriptorCollection()) {
      auto key = descriptor.first;
      if (std::find(settingsNotForwarded.begin(), settingsNotForwarded.end(), key) != settingsNotForwarded.end()) {
        continue;
      }
      if (this->_settings->valueExists(key)) {
        calculator->settings().modifyValue(key, _settings->getValue(key));
      }
    }
  }
  if (this->_settings->getBool("automated_charges")) {
    this->adjustCharges();
  }
  if (this->subsystemsAreGeneratedFromIndices()) {
    this->assignMethodToSubsystems();
  }
}
void SerenityEmbeddingCalculator::adjustCharges() {
  /*
   * We ensure that all environment systems are restricted. All alpha electron excess will
   * be formally assigned to the active system and the charge adjusted such that the total
   * charge of the supersystem is conserved. Note that the actual populations on the subsystems
   * will be determined during the embedding procedure if an automated partitioning scheme
   * is chosen.
   */
  const int totalCharge = this->_settings->getInt("molecular_charge");
  int assignedCharges = 0;
  for (unsigned int i = 1; i < _underlyingCalculators.size(); ++i) {
    auto calculator = _underlyingCalculators[i];
    auto atomCollection = calculator->getStructure();
    unsigned int totalNuclearCharge = 0;
    for (auto& element : atomCollection->getElements()) {
      totalNuclearCharge += Utils::ElementInfo::Z(element);
    }
    if (totalNuclearCharge % 2 != 0) {
      calculator->settings().modifyInt(Utils::SettingsNames::molecularCharge, 1);
      assignedCharges += 1;
    }
  }
  _underlyingCalculators[0]->settings().modifyInt(Utils::SettingsNames::molecularCharge, totalCharge - assignedCharges);
  unsigned int totalMultiplicity = this->_settings->getInt(Utils::SettingsNames::spinMultiplicity);
  _underlyingCalculators[0]->settings().modifyInt(Utils::SettingsNames::spinMultiplicity, totalMultiplicity);
}
bool SerenityEmbeddingCalculator::subsystemsAreGeneratedFromIndices() const {
  return !this->_settings->getIntListList("qmqm_atom_indices").empty();
}
void SerenityEmbeddingCalculator::subsystemAtomAssignmentSanityCheck() {
  if (!this->subsystemsAreGeneratedFromIndices()) {
    return;
  }
  if (this->_settings->getIntListList("qmqm_atom_indices").size() != _underlyingCalculators.size()) {
    throw std::runtime_error("The number of subsystem index lists and calculators does not match!");
  }
  std::vector<unsigned int> atomIndicesAssigned;
  for (const auto& assignment : this->_settings->getIntListList("qmqm_atom_indices")) {
    if (std::any_of(assignment.begin(), assignment.end(), [&](int index) { return index < 0; })) {
      throw std::runtime_error("Atom indices must be positive or zero.");
    }
    for (const auto& idx : assignment) {
      if (std::find(atomIndicesAssigned.begin(), atomIndicesAssigned.end(), idx) != atomIndicesAssigned.end()) {
        throw std::runtime_error("An atom is assigned to multiple subsystems during a QM/QM calculation."
                                 " The atom index is " +
                                 std::to_string(idx));
      }
    }
    atomIndicesAssigned.insert(atomIndicesAssigned.begin(), assignment.begin(), assignment.end());
  }
  if (atomIndicesAssigned.size() != _geometry->getNAtoms()) {
    throw std::runtime_error("Not all atoms in the QM region were assigned to a subsystem");
  }
}
void SerenityEmbeddingCalculator::assignMethodToSubsystems() {
  const std::string& fullMethod = this->_settings->getString(Utils::SettingsNames::method);
  const auto methods = Utils::CalculationRoutines::split(fullMethod, '/');
  // Maybe adjust the logic here for (QM/QM)/QM embedding, i.e., if we have more than two QM embedding layers.
  if (methods.size() != _underlyingCalculators.size()) {
    throw std::runtime_error(
        "The number of methods must match the number of underlying calculators in QM/QM embedding");
  }
  for (unsigned int iCalculator = 0; iCalculator < _underlyingCalculators.size(); ++iCalculator) {
    _underlyingCalculators[iCalculator]->settings().modifyString(Utils::SettingsNames::method, methods[iCalculator]);
  }
}
template<>
Sty::SpinPolarizedData<Sty::RESTRICTED, std::vector<unsigned int>> SerenityEmbeddingCalculator::joinOrbitalRanges(
    const Sty::SpinPolarizedData<Sty::RESTRICTED, std::vector<unsigned int>>& rangeOne,
    const Sty::SpinPolarizedData<Sty::RESTRICTED, std::vector<unsigned int>>& rangeTwo) {
  auto result = rangeOne;
  result.insert(result.end(), rangeTwo.begin(), rangeTwo.end());
  return result;
}
template<>
Sty::SpinPolarizedData<Sty::UNRESTRICTED, std::vector<unsigned int>> SerenityEmbeddingCalculator::joinOrbitalRanges(
    const Sty::SpinPolarizedData<Sty::UNRESTRICTED, std::vector<unsigned int>>& rangeOne,
    const Sty::SpinPolarizedData<Sty::UNRESTRICTED, std::vector<unsigned int>>& rangeTwo) {
  auto result = rangeOne;
  result.alpha.insert(result.alpha.end(), rangeTwo.alpha.begin(), rangeTwo.alpha.end());
  result.beta.insert(result.beta.end(), rangeTwo.beta.begin(), rangeTwo.beta.end());
  return result;
}
template<>
Sty::SPMatrix<Sty::RESTRICTED> SerenityEmbeddingCalculator::extractValenceOrbitals(
    const Sty::SPMatrix<Sty::RESTRICTED>& coefficients,
    const Sty::SpinPolarizedData<Sty::RESTRICTED, std::vector<unsigned int>>& valenceRange) {
  Sty::SPMatrix<Sty::RESTRICTED> valenceOrbitalCoefficients(Eigen::MatrixXd::Zero(coefficients.rows(), valenceRange.size()));
  unsigned int counter = 0;
  for (const auto& iValenceOrb : valenceRange) {
    valenceOrbitalCoefficients.col(counter) = coefficients.col(iValenceOrb).eval();
    ++counter;
  }
  return valenceOrbitalCoefficients;
}
template<>
Sty::SPMatrix<Sty::UNRESTRICTED> SerenityEmbeddingCalculator::extractValenceOrbitals(
    const Sty::SPMatrix<Sty::UNRESTRICTED>& coefficients,
    const Sty::SpinPolarizedData<Sty::UNRESTRICTED, std::vector<unsigned int>>& valenceRange) {
  Sty::SPMatrix<Sty::UNRESTRICTED> valenceOrbitalCoefficients;
  valenceOrbitalCoefficients.alpha.resize(coefficients.alpha.rows(), valenceRange.alpha.size());
  valenceOrbitalCoefficients.beta.resize(coefficients.beta.rows(), valenceRange.beta.size());
  valenceOrbitalCoefficients.alpha.setZero();
  valenceOrbitalCoefficients.beta.setZero();
  unsigned int counter = 0;
  for (const auto& iValenceOrb : valenceRange.alpha) {
    valenceOrbitalCoefficients.alpha.col(counter) = coefficients.alpha.col(iValenceOrb).eval();
    ++counter;
  }
  counter = 0;
  for (const auto& iValenceOrb : valenceRange.beta) {
    valenceOrbitalCoefficients.beta.col(counter) = coefficients.beta.col(iValenceOrb).eval();
    ++counter;
  }
  return valenceOrbitalCoefficients;
}
template<>
Eigen::MatrixXd SerenityEmbeddingCalculator::compactIntegrals(const Sty::SPMatrix<Sty::RESTRICTED>& spMatrix) {
  return Eigen::MatrixXd(spMatrix);
}
template<>
Eigen::MatrixXd SerenityEmbeddingCalculator::compactIntegrals(const Sty::SPMatrix<Sty::UNRESTRICTED>& spMatrix) {
  /*
   * In the case of an unrestricted reference, we store the data as blocks below each other:
   *
   * M_all = M_alpha
   *         M_beta
   */
  const unsigned int nRows = spMatrix.alpha.rows() + spMatrix.beta.rows();
  const unsigned int nCols = std::max(spMatrix.alpha.cols(), spMatrix.beta.cols());
  Eigen::MatrixXd compactedMatrix = Eigen::MatrixXd::Zero(nRows, nCols);
  compactedMatrix.topLeftCorner(spMatrix.alpha.rows(), spMatrix.alpha.cols()) = spMatrix.alpha;
  compactedMatrix.bottomLeftCorner(spMatrix.beta.rows(), spMatrix.beta.cols()) = spMatrix.beta;
  return compactedMatrix;
}
template<Sty::Options::SCF_MODES SCFMode>
void SerenityEmbeddingCalculator::prepareCAS() {
  if (this->_settings->getIntList("cas_systems").empty()) {
    return;
  }
  Sty::printSectionTitle("Prepare CAS Integrals");
  for (const auto& calculatorIndex : this->_settings->getIntList("cas_systems")) {
    auto activeCalculator = _underlyingCalculators[calculatorIndex];
    auto activeSystem = activeCalculator->getSystemController();
    std::vector<std::shared_ptr<Sty::SystemController>> environmentSystems;
    for (const auto& calculator : _underlyingCalculators) {
      if (calculator != activeCalculator) {
        environmentSystems.push_back(calculator->getSystemController());
      }
    }
    const auto fullValenceRange = activeSystem->template getActiveOrbitalController<SCFMode>()->getAllValenceOrbitalIndices();
    const auto coreOrbitals = activeSystem->template getActiveOrbitalController<SCFMode>()
                                  ->getValenceOrbitalIndices(activeSystem->template getNOccupiedOrbitals<SCFMode>())
                                  .second;
    const auto coreAndValenceOrbitalRanges = joinOrbitalRanges(coreOrbitals, fullValenceRange);

    Sty::FCIDumpFileWriterTaskSettings settings;
    this->updateEmbeddingSettings(settings.embedding);
    Sty::FCIDumpFileWriter<SCFMode> fciDumpFileWriter(activeSystem, environmentSystems, settings);
    const auto oneParticleIntegrals = fciDumpFileWriter.getOneParticleIntegrals(coreAndValenceOrbitalRanges);
    const auto& coefficientMatrix = activeSystem->template getActiveOrbitalController<SCFMode>()->getCoefficients();

    Sty::SPMatrix<SCFMode> valenceOrbitalCoefficients = extractValenceOrbitals(coefficientMatrix, coreAndValenceOrbitalRanges);
    auto& results = activeCalculator->results();
    results.set<Utils::Property::OneElectronMatrix>(compactIntegrals(oneParticleIntegrals));
    results.set<Utils::Property::CoefficientMatrix>(
        SerenityConversionFunctions::molecularOrbitalsFromCoefficients(coefficientMatrix));

    Sty::OrbitalsIOTask<SCFMode> orbitalsIoTask(activeSystem);
    orbitalsIoTask.settings.fileFormat = Sty::Options::ORBITAL_FILE_TYPES::MOLDEN;
    orbitalsIoTask.settings.write = true;
    orbitalsIoTask.run();

    std::unordered_map<std::string, double> partialEnergies;
    partialEnergies.insert(std::make_pair("total_uncorrelated_energy", fciDumpFileWriter.getTotalUncorrelatedEnergy()));
    partialEnergies.insert(std::make_pair("uncorrelated_active_space_energy",
                                          fciDumpFileWriter.getUncorrelatedActiveSpaceEnergy(coreAndValenceOrbitalRanges)));
    partialEnergies.insert(
        std::make_pair("total_core_energy", fciDumpFileWriter.getTotalCoreEnergy(coreAndValenceOrbitalRanges)));
    results.set<Utils::Property::PartialEnergies>(partialEnergies);

    activeCalculator->setRequiredProperties(Utils::Property::OrbitalFragmentPopulations |
                                            Utils::Property::NAlphaElectrons | Utils::Property::NBetaElectrons);
    activeCalculator->storeProperties<SCFMode>();
  }
}
template void SerenityEmbeddingCalculator::prepareCAS<Sty::RESTRICTED>();
template void SerenityEmbeddingCalculator::prepareCAS<Sty::UNRESTRICTED>();

template<Sty::Options::SCF_MODES SCFMode>
void SerenityEmbeddingCalculator::runStaticEmbeddingForMultipleSubsystems() {
  std::vector<std::shared_ptr<Sty::SystemController>> activeSystems;
  std::vector<std::shared_ptr<CalculatorBase>> activeCalculators;
  std::vector<std::shared_ptr<Sty::SystemController>> environmentSystems;
  unsigned int counter = 0;
  for (auto& calculator : _underlyingCalculators) {
    const auto& casSystems = this->_settings->getIntList("cas_systems");
    if (isLocalCorrelation(calculator) || std::find(casSystems.begin(), casSystems.end(), counter) != casSystems.end()) {
      activeSystems.push_back(calculator->getSystemController());
      activeCalculators.push_back(calculator);
    }
    else {
      environmentSystems.push_back(calculator->getSystemController());
    }
    ++counter;
  }
  Sty::TopDownStaticEmbeddingTask<SCFMode> tdStaticEmbeddingTask(activeSystems, environmentSystems,
                                                                 loadSupersystemIfGiven<SCFMode>());
  this->updateLocalizationSettings(tdStaticEmbeddingTask.settings.loc);
  this->updateEmbeddingSettings(tdStaticEmbeddingTask.settings.lcSettings.embeddingSettings);
  for (const auto& calculator : activeCalculators) {
    if (isLocalCorrelation(calculator)) {
      updateLocalCorrelationSettings(tdStaticEmbeddingTask.settings.lcSettings, calculator->settings(), true);
      break;
    }
  }
  tdStaticEmbeddingTask.settings.split.orbitalThreshold = this->_settings->getDouble("orbital_threshold");
  auto value = this->_settings->getString("system_partitioning");
  Sty::Options::resolve(value, tdStaticEmbeddingTask.settings.split.systemPartitioning);
  tdStaticEmbeddingTask.run();
  _finalEnergy = std::make_unique<double>(tdStaticEmbeddingTask.getFinalEnergy());
}
template<Sty::Options::SCF_MODES SCFMode>
void SerenityEmbeddingCalculator::runTopDownEmbeddingWithSingleSystemLocalCorrelation() {
  std::vector<std::shared_ptr<Sty::SystemController>> allOtherSystems;
  auto activeSystem = _underlyingCalculators[0]->getSystemController();
  for (unsigned int i = 1; i < _underlyingCalculators.size(); ++i) {
    allOtherSystems.insert(allOtherSystems.begin(), _underlyingCalculators[i]->getSystemController());
  }
  Sty::DFTEmbeddedLocalCorrelationTask lcTdTask(activeSystem, allOtherSystems, loadSupersystemIfGiven<SCFMode>());
  this->updateLocalizationSettings(lcTdTask.settings.loc);
  this->updateEmbeddingSettings(lcTdTask.settings.lcSettings.embeddingSettings);
  updateLocalCorrelationSettings(lcTdTask.settings.lcSettings, _underlyingCalculators[0]->settings(), true);
  this->updateBasisSetTruncationSettings(lcTdTask.settings.trunc);
  lcTdTask.settings.split.orbitalThreshold = this->_settings->getDouble("orbital_threshold");
  auto value = this->_settings->getString("system_partitioning");
  Sty::Options::resolve(value, lcTdTask.settings.split.systemPartitioning);
  lcTdTask.settings.add.addOccupiedOrbitals = false;
  lcTdTask.run();
  double correlationEnergy = this->template getCorrelationEnergy<SCFMode>(activeSystem, lcTdTask.settings.lcSettings.method);
  double referenceEnergy = this->template getReferenceEnergy<SCFMode>(activeSystem);
  _finalEnergy = std::make_unique<double>(correlationEnergy + referenceEnergy);
}
template<Sty::Options::SCF_MODES SCFMode>
void SerenityEmbeddingCalculator::runTopDownEmbeddingOnlyMeanField() {
  std::vector<std::shared_ptr<Sty::SystemController>> allOtherSystems;
  auto activeSystem = _underlyingCalculators[0]->getSystemController();
  for (unsigned int i = 1; i < _underlyingCalculators.size(); ++i) {
    allOtherSystems.insert(allOtherSystems.begin(), _underlyingCalculators[i]->getSystemController());
  }
  Sty::Settings environmentSettings = allOtherSystems[0]->getSettings();
  auto totalEnvironment = std::make_shared<Sty::SystemController>(std::make_shared<Sty::Geometry>(), environmentSettings);
  Sty::SystemAdditionTask<SCFMode> additionTask(totalEnvironment, allOtherSystems);
  additionTask.settings.addOccupiedOrbitals = false;
  additionTask.run();

  Sty::TDEmbeddingTask<SCFMode> tdTask(activeSystem, totalEnvironment);
  this->updateEmbeddingSettings(tdTask.settings.embedding);
  updateLocalCorrelationSettings(tdTask.settings.lcSettings, _underlyingCalculators[0]->settings());
  auto value = this->_settings->getString("trunc_algorithm");
  Sty::Options::resolve(value, tdTask.settings.truncAlgorithm);
  tdTask.settings.netThreshold = this->_settings->getDouble("net_threshold");
  tdTask.settings.load = this->_settings->getString("supersystem_path");
  Sty::LocalizationTaskSettings localizationTaskSettings;
  this->updateLocalizationSettings(localizationTaskSettings);
  tdTask.settings.locType = localizationTaskSettings.locType;
  tdTask.settings.splitValenceAndCore = localizationTaskSettings.splitValenceAndCore;
  tdTask.settings.useFermiLevel = false;
  value = this->_settings->getString("system_partitioning");
  Sty::Options::resolve(value, tdTask.settings.systemPartitioning);
  tdTask.settings.orbitalThreshold = this->_settings->getDouble("orbital_threshold");
  tdTask.run();
  _finalEnergy = std::make_unique<double>(this->template getReferenceEnergy<SCFMode>(activeSystem));
}
template<Sty::Options::SCF_MODES SCFMode>
double SerenityEmbeddingCalculator::getReferenceEnergy(std::shared_ptr<Sty::SystemController> activeSystemController) {
  auto energyComponentController =
      activeSystemController->template getElectronicStructure<SCFMode>()->getEnergyComponentController();
  double referenceEnergy = NAN;
  if (energyComponentController->checkEnergyComponentFromChildren(Sty::ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_WF_DFT)) {
    referenceEnergy = energyComponentController->getEnergyComponent(Sty::ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_WF_DFT);
  }
  else if (energyComponentController->checkEnergyComponentFromChildren(Sty::ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_DFT_DFT)) {
    referenceEnergy = energyComponentController->getEnergyComponent(Sty::ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_DFT_DFT);
  }
  else if (energyComponentController->checkEnergyComponentFromChildren(Sty::ENERGY_CONTRIBUTIONS::PBE_SUPERSYSTEM_ENERGY_WFT_DFT)) {
    referenceEnergy = energyComponentController->getEnergyComponent(Sty::ENERGY_CONTRIBUTIONS::PBE_SUPERSYSTEM_ENERGY_WFT_DFT);
  }
  else if (energyComponentController->checkEnergyComponentFromChildren(Sty::ENERGY_CONTRIBUTIONS::PBE_SUPERSYSTEM_ENERGY_DFT_DFT)) {
    referenceEnergy = energyComponentController->getEnergyComponent(Sty::ENERGY_CONTRIBUTIONS::PBE_SUPERSYSTEM_ENERGY_DFT_DFT);
  }
  else {
    throw std::runtime_error(
        "No embedded energy is available from Serenity after QM/QM calculation. Something went wrong!");
  }
  return referenceEnergy;
}
template<Sty::Options::SCF_MODES SCFMode>
double SerenityEmbeddingCalculator::getCorrelationEnergy(std::shared_ptr<Sty::SystemController> activeSystemController,
                                                         const Sty::Options::PNO_METHOD& pnoMethod) {
  double correlationEnergy = 0.0;
  auto electronicStructure = activeSystemController->template getElectronicStructure<SCFMode>();
  switch (pnoMethod) {
    case Sty::Options::PNO_METHOD::SC_MP2:
    case Sty::Options::PNO_METHOD::DLPNO_MP2: {
      correlationEnergy = electronicStructure->getEnergy(Sty::ENERGY_CONTRIBUTIONS::MP2_CORRECTION);
      break;
    }
    case Sty::Options::PNO_METHOD::DLPNO_CCSD: {
      correlationEnergy = electronicStructure->getEnergy(Sty::ENERGY_CONTRIBUTIONS::CCSD_CORRECTION);
      break;
    }
    case Sty::Options::PNO_METHOD::DLPNO_CCSD_T0: {
      correlationEnergy = electronicStructure->getEnergy(Sty::ENERGY_CONTRIBUTIONS::CCSD_CORRECTION);
      correlationEnergy += electronicStructure->getEnergy(Sty::ENERGY_CONTRIBUTIONS::TRIPLES_CORRECTION);
      break;
    }
    case Sty::Options::PNO_METHOD::NONE: {
      correlationEnergy = 0.0;
      break;
    }
  }
  return correlationEnergy;
}

template<Sty::Options::SCF_MODES SCFMode>
void SerenityEmbeddingCalculator::runEmbeddingTask() {
  if (this->_settings->getBool("static_embedding")) {
    this->runStaticEmbeddingForMultipleSubsystems<SCFMode>();
    return;
  }
  if (isLocalCorrelation(_underlyingCalculators[0])) {
    this->runTopDownEmbeddingWithSingleSystemLocalCorrelation<SCFMode>();
    return;
  }
  this->runTopDownEmbeddingOnlyMeanField<SCFMode>();
}
template void SerenityEmbeddingCalculator::runEmbeddingTask<Sty::RESTRICTED>();
template void SerenityEmbeddingCalculator::runEmbeddingTask<Sty::UNRESTRICTED>();

template<Sty::Options::SCF_MODES SCFMode>
std::shared_ptr<Sty::SystemController> SerenityEmbeddingCalculator::loadSupersystemIfGiven() {
  std::shared_ptr<Sty::SystemController> supersystem;
  if (!this->_settings->getString("supersystem_path").empty()) {
    Sty::Settings supersystemSettings;
    supersystemSettings.load = this->_settings->getString("supersystem_path");
    supersystem = std::make_shared<Sty::SystemController>(supersystemSettings);
    return supersystem;
  }
  if (!this->_settings->getString("external_supersystem_mo_file").empty()) {
    if (_underlyingCalculators.empty())
      throw std::runtime_error("ERROR: At least one system must be available in order to built a supersystem.\n");
    Sty::Settings supersystemSettings =
        _underlyingCalculators[_underlyingCalculators.size() - 1]->getSystemController()->getSettings();
    supersystemSettings.name = "TMP_Supersystem";
    supersystemSettings.path = supersystemSettings.path + supersystemSettings.name + "/";
    supersystemSettings.charge = 0;
    supersystemSettings.spin = 0;
    supersystem = std::make_shared<Sty::SystemController>(std::make_shared<Sty::Geometry>(), supersystemSettings);
    // Addition
    std::vector<std::shared_ptr<Sty::SystemController>> allSystems;
    for (auto& calculator : _underlyingCalculators) {
      allSystems.push_back(calculator->getSystemController());
    }
    Sty::SystemAdditionTask<SCFMode> additionTask(supersystem, allSystems);
    additionTask.settings.addOccupiedOrbitals = false;
    additionTask.settings.checkSuperCharge = false;
    additionTask.settings.checkSuperGeom = false;
    additionTask.run();
    if (this->subsystemsAreGeneratedFromIndices()) {
      // In the case that the supersystem was generated from subsystem indices, we have no idea of the atom ordering
      // in the underlying calculators. Therefore, we directly take the _geometry object still holding the atom order
      // of the original input system.
      Sty::Settings updatedSupersystemSettings = supersystem->getSettings();
      supersystem = std::make_shared<Sty::SystemController>(_geometry, updatedSupersystemSettings);
    }

    // Load orbitals
    Sty::OrbitalsIOTask<SCFMode> orbitalsIOTask(supersystem);
    orbitalsIOTask.settings.path = this->_settings->getString("external_supersystem_mo_file");
    auto value = this->_settings->getString("external_supersystem_mo_type");
    Sty::Options::resolve(value, orbitalsIOTask.settings.fileFormat);
    orbitalsIOTask.settings.write = false;
    orbitalsIOTask.run();
    return supersystem;
  }
  return supersystem;
}

} /* namespace Serenity */
} /* namespace Scine */
