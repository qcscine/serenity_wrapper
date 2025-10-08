/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
/* Wrapper Includes */
#include "Serenity/Calculators/CalculatorBase.h"
#include "Serenity/Calculators/ScineSettings.h"
#include "Serenity/Calculators/SerenityState.h"
#include "Serenity/Utilities/SerenityConversionFunctions.h"
/* Serenity Includes */
#include <analysis/populationAnalysis/HirshfeldPopulationCalculator.h>
#include <analysis/populationAnalysis/MullikenPopulationCalculator.h>
#include <basis/AtomCenteredBasisController.h>
#include <data/ElectronicStructure.h>
#include <data/OrbitalController.h>
#include <data/grid/BasisFunctionOnGridController.h>
#include <data/grid/BasisFunctionOnGridControllerFactory.h>
#include <data/grid/DensityMatrixDensityOnGridController.h>
#include <data/grid/DensityOnGridCalculator.h>
#include <data/matrices/DensityMatrix.h>
#include <geometry/Geometry.h>
#include <integrals/OneElectronIntegralController.h>
#include <integrals/wrappers/Libint.h>
#include <io/FormattedOutputStream.h>
#include <math/Matrix.h>
#include <potentials/HCorePotential.h>
#include <system/SystemController.h>
/* Scine Includes */
#include <Utils/Geometry.h>
#include <Utils/Solvation/ImplicitSolvation.h>
#include <Utils/Technical/UniqueIdentifier.h>
#include <Utils/Typenames.h>

using namespace Serenity;

namespace Scine {
namespace Serenity {

CalculatorBase::CalculatorBase()
  : _results(std::make_unique<Scine::Utils::Results>()),
    _system(nullptr),
    _geometry(nullptr),
    _scinePositions(nullptr),
    _moved(true) {
  this->_settings = std::make_unique<ScineSettings>();
  auto& libint = Libint::getInstance();
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 2);
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 3);
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 4);
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 1, 2);
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 1, 3);
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 1, 4);
}

CalculatorBase::~CalculatorBase() {
  auto& libint = Libint::getInstance();
  libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 2);
  libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 3);
  libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 4);
  libint.freeEngines(LIBINT_OPERATOR::coulomb, 1, 2);
  libint.freeEngines(LIBINT_OPERATOR::coulomb, 1, 3);
  libint.freeEngines(LIBINT_OPERATOR::coulomb, 1, 4);
}

CalculatorBase::CalculatorBase(const CalculatorBase& other) {
  _system = nullptr;
  _settings = std::make_unique<ScineSettings>(*other._settings);
  _moved = other._moved;
  if (other._geometry) {
    _geometry = std::make_shared<Geometry>(other._geometry->getAtomSymbols(), *other._scinePositions);
  }
  else {
    _geometry = nullptr;
  }
  if (other._scinePositions) {
    _scinePositions = std::make_unique<Scine::Utils::PositionCollection>(*other._scinePositions);
  }
  else {
    _scinePositions = nullptr;
  }
  if (other._results) {
    _results = std::make_unique<Scine::Utils::Results>(*other._results);
  }
  else {
    _results = std::make_unique<Scine::Utils::Results>();
  }
  auto& libint = Libint::getInstance();
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 2);
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 3);
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 4);
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 1, 2);
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 1, 3);
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 1, 4);
}

Scine::Utils::Settings& CalculatorBase::settings() {
  return *_settings;
}

const Scine::Utils::Settings& CalculatorBase::settings() const {
  return *_settings;
}

Scine::Utils::PropertyList CalculatorBase::getRequiredProperties() const {
  return _requiredProperties;
}

void CalculatorBase::setStructure(const Scine::Utils::AtomCollection& structure) {
  _geometry = std::make_shared<Geometry>(SerenityConversionFunctions::atomCollectionToGeometry(structure));
  _scinePositions = std::make_unique<Scine::Utils::PositionCollection>(structure.getPositions());
  // TODO
  //  if (_system != nullptr)
  //    remove_all(_system->getSettings().path);
  _system = nullptr;
  _results = std::make_unique<Scine::Utils::Results>();
}

std::unique_ptr<Scine::Utils::AtomCollection> CalculatorBase::getStructure() const {
  if (!_geometry || !_scinePositions) {
    throw std::runtime_error("Missing geometry in a Serenity Calculator");
  }
  return std::make_unique<Scine::Utils::AtomCollection>(SerenityConversionFunctions::geometryToAtomCollection(*_geometry));
}

void CalculatorBase::modifyPositions(Scine::Utils::PositionCollection newPositions) {
  if (!_geometry || !_scinePositions) {
    throw std::runtime_error("Missing geometry in a Serenity Calculator");
  };
  auto diff = ((*_scinePositions) - newPositions).rowwise().norm();
  if (_system && diff.maxCoeff() > 0.1) {
    _system->setElectronicStructure<RESTRICTED>(nullptr);
    _system->setElectronicStructure<UNRESTRICTED>(nullptr);
  }
  (*_scinePositions) = Scine::Utils::PositionCollection(newPositions);
  auto old = iOOptions.printGridInfo;
  iOOptions.printGridInfo = false;
  _geometry->setCoordinates(newPositions);

  iOOptions.printGridInfo = old;
  _results = std::make_unique<Scine::Utils::Results>();
  this->_moved = true;
}

const Scine::Utils::PositionCollection& CalculatorBase::getPositions() const {
  if (!_scinePositions) {
    throw std::runtime_error("Missing geometry in a Serenity Calculator");
  }
  return *_scinePositions;
}

void CalculatorBase::setRequiredProperties(const Scine::Utils::PropertyList& requiredProperties) {
  _requiredProperties = requiredProperties;
}

Scine::Utils::Results& CalculatorBase::results() {
  return *_results;
}

const Scine::Utils::Results& CalculatorBase::results() const {
  return *_results;
}

void CalculatorBase::loadState(std::shared_ptr<Scine::Core::State> state) {
  auto castState = std::dynamic_pointer_cast<SerenityState>(state);
  if (!castState)
    throw Scine::Core::StateCastingException();

  //  TODO
  //  // Remove old system
  //  if (_system != nullptr)
  //    remove_all(_system->getSettings().path);
  //  _system = nullptr;

  // Load state as new system
  auto settings = Settings();
  // throws error for wrong input and updates 'any' entries
  Utils::Solvation::ImplicitSolvation::solvationNeededAndPossible(availableSolvationModels(), *_settings);
  _settings->applyTo(settings);
  this->applyFixedSettings(settings);
  // Generate a unique name
  Scine::Utils::UniqueIdentifier uid;
  settings.name = uid.getStringRepresentation();

  _geometry = std::make_shared<Geometry>(castState->system->getGeometry()->getAtomSymbols(),
                                         castState->system->getGeometry()->getCoordinates());
  _system = std::make_shared<SystemController>(_geometry, settings);

  // Load data into the new system generated from the state
  if (castState->system->hasElectronicStructure<RESTRICTED>()) {
    CoefficientMatrix<RESTRICTED> coeff(_system->getBasisController());
    auto orig = castState->system->getElectronicStructure<RESTRICTED>()->getMolecularOrbitals()->getCoefficients();
    (Eigen::MatrixXd) coeff = (Eigen::MatrixXd)orig;
    auto eval = castState->system->getElectronicStructure<RESTRICTED>()->getMolecularOrbitals()->getEigenvalues();
    auto nCoreOrbs = castState->system->getElectronicStructure<RESTRICTED>()->getMolecularOrbitals()->getNCoreOrbitals();
    auto orbitals = std::make_shared<OrbitalController<RESTRICTED>>(_system->getBasisController(), nCoreOrbs);
    orbitals->updateOrbitals(coeff, eval);
    auto es = std::make_shared<ElectronicStructure<RESTRICTED>>(orbitals, _system->getOneElectronIntegralController(),
                                                                castState->system->getNOccupiedOrbitals<RESTRICTED>());
    _system->setElectronicStructure<RESTRICTED>(es);
  }
  if (castState->system->hasElectronicStructure<UNRESTRICTED>()) {
    CoefficientMatrix<UNRESTRICTED> coeff(_system->getBasisController());
    auto orig = castState->system->getElectronicStructure<UNRESTRICTED>()->getMolecularOrbitals()->getCoefficients();
    coeff.alpha = orig.alpha;
    coeff.beta = orig.beta;
    auto eval = castState->system->getElectronicStructure<UNRESTRICTED>()->getMolecularOrbitals()->getEigenvalues();
    auto nCoreOrbs = castState->system->getElectronicStructure<UNRESTRICTED>()->getMolecularOrbitals()->getNCoreOrbitals();
    auto orbitals = std::make_shared<OrbitalController<UNRESTRICTED>>(_system->getBasisController(), nCoreOrbs);
    orbitals->updateOrbitals(coeff, eval);
    auto es = std::make_shared<ElectronicStructure<UNRESTRICTED>>(orbitals, _system->getOneElectronIntegralController(),
                                                                  castState->system->getNOccupiedOrbitals<UNRESTRICTED>());
    _system->setElectronicStructure<UNRESTRICTED>(es);
  }
  _results = std::make_unique<Scine::Utils::Results>();
}

std::shared_ptr<Scine::Core::State> CalculatorBase::getState() const {
  if (!_geometry) {
    throw std::runtime_error("Missing geometry in Serenity DFT Calculator");
  };
  auto settings = Settings();
  // throws error for wrong input and updates 'any' entries
  Utils::Solvation::ImplicitSolvation::solvationNeededAndPossible(availableSolvationModels(), *_settings);
  _settings->applyTo(settings);
  this->applyFixedSettings(settings);
  // Generate a unique name
  Scine::Utils::UniqueIdentifier uid;
  settings.name = uid.getStringRepresentation();
  auto geometry = std::make_shared<Geometry>(_geometry->getAtomSymbols(), _geometry->getCoordinates());
  auto system = std::make_shared<SystemController>(geometry, settings);

  if (_system) {
    if (_system->hasElectronicStructure<RESTRICTED>()) {
      CoefficientMatrix<RESTRICTED> coeff(system->getBasisController());
      (Eigen::MatrixXd) coeff =
          (Eigen::MatrixXd)_system->getElectronicStructure<RESTRICTED>()->getMolecularOrbitals()->getCoefficients();
      auto eval = _system->getElectronicStructure<RESTRICTED>()->getMolecularOrbitals()->getEigenvalues();
      auto nCoreOrbs = _system->getElectronicStructure<RESTRICTED>()->getMolecularOrbitals()->getNCoreOrbitals();
      auto orbitals = std::make_shared<OrbitalController<RESTRICTED>>(system->getBasisController(), nCoreOrbs);
      orbitals->updateOrbitals(coeff, eval);
      auto es = std::make_shared<ElectronicStructure<RESTRICTED>>(orbitals, system->getOneElectronIntegralController(),
                                                                  _system->getNOccupiedOrbitals<RESTRICTED>());
      system->setElectronicStructure<RESTRICTED>(es);
    }
    if (_system->hasElectronicStructure<UNRESTRICTED>()) {
      CoefficientMatrix<UNRESTRICTED> coeff(system->getBasisController());
      auto orig = _system->getElectronicStructure<UNRESTRICTED>()->getMolecularOrbitals()->getCoefficients();
      coeff.alpha = orig.alpha;
      coeff.beta = orig.beta;
      auto eval = _system->getElectronicStructure<UNRESTRICTED>()->getMolecularOrbitals()->getEigenvalues();
      auto nCoreOrbs = _system->getElectronicStructure<UNRESTRICTED>()->getMolecularOrbitals()->getNCoreOrbitals();
      auto orbitals = std::make_shared<OrbitalController<UNRESTRICTED>>(system->getBasisController(), nCoreOrbs);
      orbitals->updateOrbitals(coeff, eval);
      auto es = std::make_shared<ElectronicStructure<UNRESTRICTED>>(orbitals, system->getOneElectronIntegralController(),
                                                                    _system->getNOccupiedOrbitals<UNRESTRICTED>());
      system->setElectronicStructure<UNRESTRICTED>(es);
    }
  }
  system->setDiskMode(true);
  return std::make_shared<SerenityState>(system);
}

const Scine::Utils::Results& CalculatorBase::calculate(std::string /*description*/) {
  if (!_geometry) {
    throw std::runtime_error("Missing geometry in Serenity Calculator");
  };

  if (!this->possibleProperties().containsSubSet(_requiredProperties)) {
    throw std::runtime_error("Unavailable Properties requested.");
  }

  // Modify output level
  bool showOutput = this->_settings->getBool("show_serenity_output");
  if (!showOutput) {
    GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::MINIMUM;
    iOOptions.printFinalOrbitalEnergies = false;
    iOOptions.printGeometry = false;
    iOOptions.printSCFCycleInfo = false;
    iOOptions.printSCFResults = false;
    iOOptions.printDebugInfos = false;
    iOOptions.printGridInfo = false;
    iOOptions.gridAccuracyCheck = false;
    iOOptions.timingsPrintLevel = 0;
  }
  _system = this->getSystemController();

  // Initialize the results
  _results = std::make_unique<Scine::Utils::Results>();

  // Run the actual calculation
  try {
    if (_system->getSettings().scfMode == RESTRICTED) {
      this->calculateImplRestricted();
    }
    else {
      this->calculateImplUnrestricted();
    }
  }
  catch (SerenityError& e) {
    throw Core::UnsuccessfulCalculationException(e.what());
  }

  // Reset output
  if (!showOutput) {
    iOOptions = IOOptions();
  }

  _results->set<Scine::Utils::Property::ProgramName>("serenity");

  return *_results;
}

template<>
Scine::Utils::DensityMatrix CalculatorBase::convertDensityMatrix(DensityMatrix<RESTRICTED> dmat,
                                                                 SpinPolarizedData<RESTRICTED, unsigned int, void> nEl) const {
  Scine::Utils::DensityMatrix ret;
  ret.setDensity(Eigen::MatrixXd(dmat), nEl);
  return ret;
}

template<>
Scine::Utils::DensityMatrix CalculatorBase::convertDensityMatrix(DensityMatrix<UNRESTRICTED> dmat,
                                                                 SpinPolarizedData<UNRESTRICTED, unsigned int, void> nEl) const {
  Scine::Utils::DensityMatrix ret;
  ret.setDensity(Eigen::MatrixXd(dmat.alpha), Eigen::MatrixXd(dmat.beta), nEl.alpha, nEl.beta);
  return ret;
}

template<Options::SCF_MODES ScfMode>
std::vector<double> CalculatorBase::getMullikenCharges() const {
  std::vector<double> charges;
  MullikenPopulationCalculator<ScfMode> calculator;
  auto populations = calculator.calculateMullikenPopulations(_system);
  return populationToCharges<ScfMode>(populations);
}

template<Options::SCF_MODES ScfMode>
std::vector<double> CalculatorBase::getHirshfeldCharges() const {
  auto basFuncOnGridController = BasisFunctionOnGridControllerFactory::produce(128, 0.0, 2, _system->getBasisController(),
                                                                               _system->getGridController());
  auto densOnGridCalc = std::make_shared<DensityOnGridCalculator<ScfMode>>(basFuncOnGridController, 0.0);
  auto densMatController = _system->getElectronicStructure<ScfMode>()->getDensityMatrixController();
  auto densOnGridController =
      std::make_shared<DensityMatrixDensityOnGridController<ScfMode>>(densOnGridCalc, densMatController);

  HirshfeldPopulationCalculator<ScfMode> hirshfeldCalculator(_system, densOnGridController);
  auto populations = hirshfeldCalculator.getAtomPopulations();
  return populationToCharges<ScfMode>(populations);
}

template<Options::SCF_MODES ScfMode>
std::vector<double> CalculatorBase::populationToCharges(const SpinPolarizedData<ScfMode, Eigen::VectorXd>& populations) const {
  std::vector<double> charges;
  const auto& serenityAtoms = _system->getAtoms();
  for (unsigned long i = 0; i < serenityAtoms.size(); ++i) {
    double pop = electronPopulationAtAtom(populations, i);
    charges.push_back(serenityAtoms[i]->getEffectiveCharge() - pop);
  }
  return charges;
}

template std::vector<double> CalculatorBase::populationToCharges<Options::SCF_MODES::RESTRICTED>(
    const SpinPolarizedData<Options::SCF_MODES::RESTRICTED, Eigen::VectorXd>&) const;
template std::vector<double> CalculatorBase::populationToCharges<Options::SCF_MODES::UNRESTRICTED>(
    const SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::VectorXd>&) const;
template std::vector<double> CalculatorBase::getMullikenCharges<Options::SCF_MODES::RESTRICTED>() const;
template std::vector<double> CalculatorBase::getMullikenCharges<Options::SCF_MODES::UNRESTRICTED>() const;
template std::vector<double> CalculatorBase::getHirshfeldCharges<Options::SCF_MODES::RESTRICTED>() const;
template std::vector<double> CalculatorBase::getHirshfeldCharges<Options::SCF_MODES::UNRESTRICTED>() const;

void CalculatorBase::storeCoreHamiltonian() {
  // reroute output
  std::ofstream out(_system->getSettings().path + "/hcore.cout.txt");
  std::streambuf* coutbuf = std::cout.rdbuf();
  std::cout.rdbuf(out.rdbuf());
  // calculate
  try {
    Sty::HCorePotential<RESTRICTED> hCorePotential(this->_system);
    const auto& hCoreMatrix = hCorePotential.getMatrix();
    _results->set<Scine::Utils::Property::OneElectronMatrix>(hCoreMatrix);
  }
  catch (Sty::SerenityError& e) {
    throw Core::UnsuccessfulCalculationException(e.what());
  }
  // reset output
  std::cout.rdbuf(coutbuf);
}

template<Sty::Options::SCF_MODES ScfMode>
void CalculatorBase::storeDensityMatrix() {
  auto dmat = _system->getElectronicStructure<ScfMode>()->getDensityMatrix();
  _results->set<Scine::Utils::Property::DensityMatrix>(this->convertDensityMatrix(dmat, _system->getNElectrons<ScfMode>()));
}
template void CalculatorBase::storeDensityMatrix<RESTRICTED>();
template void CalculatorBase::storeDensityMatrix<UNRESTRICTED>();

void CalculatorBase::storeAOToAtomMapping() {
  auto indices = _system->getAtomCenteredBasisController()->getBasisIndices();
  Scine::Utils::AtomsOrbitalsIndexes counts(indices.size());
  for (unsigned int i = 0; i < indices.size(); i++) {
    counts.addAtom(indices[i].second - indices[i].first);
  }
  _results->set<Scine::Utils::Property::AOtoAtomMapping>(counts);
}

template<Sty::Options::SCF_MODES ScfMode>
void CalculatorBase::storeAtomicCharges() {
  auto charges = getMullikenCharges<ScfMode>();
  _results->set<Scine::Utils::Property::AtomicCharges>(charges);
}
template void CalculatorBase::storeAtomicCharges<RESTRICTED>();
template void CalculatorBase::storeAtomicCharges<UNRESTRICTED>();

void CalculatorBase::storeOverlapMatrix() {
  _results->set<Scine::Utils::Property::OverlapMatrix>(_system->getOneElectronIntegralController()->getOverlapIntegrals());
}

template<Sty::Options::SCF_MODES ScfMode>
void CalculatorBase::storeElectronicOccupations() {
  auto occupation = Scine::Utils::LcaoUtils::ElectronicOccupation();
  if (ScfMode == Sty::RESTRICTED) {
    auto nElectrons = _system->getNElectrons<Sty::RESTRICTED>();
    occupation.fillLowestRestrictedOrbitalsWithElectrons(nElectrons);
  }
  else {
    auto nElectrons = _system->getNElectrons<Sty::UNRESTRICTED>();
    occupation.fillLowestUnrestrictedOrbitals(nElectrons.alpha, nElectrons.beta);
  }
  _results->set<Scine::Utils::Property::ElectronicOccupation>(occupation);
}
template void CalculatorBase::storeElectronicOccupations<RESTRICTED>();
template void CalculatorBase::storeElectronicOccupations<UNRESTRICTED>();

template<Sty::Options::SCF_MODES ScfMode>
void CalculatorBase::storeElectronicEnergy() {
  _results->set<Scine::Utils::Property::Energy>(_system->getElectronicStructure<ScfMode>()->getEnergy());
}
template void CalculatorBase::storeElectronicEnergy<RESTRICTED>();
template void CalculatorBase::storeElectronicEnergy<UNRESTRICTED>();

template<>
void CalculatorBase::storeNElectrons<RESTRICTED>() {
  auto nElectrons = _system->getNElectrons<RESTRICTED>();
  _results->set<Scine::Utils::Property::NAlphaElectrons>(nElectrons);
  _results->set<Scine::Utils::Property::NBetaElectrons>(nElectrons);
}
template<>
void CalculatorBase::storeNElectrons<UNRESTRICTED>() {
  auto nElectrons = _system->getNElectrons<UNRESTRICTED>();
  _results->set<Scine::Utils::Property::NAlphaElectrons>(nElectrons.alpha);
  _results->set<Scine::Utils::Property::NBetaElectrons>(nElectrons.beta);
}

template<Sty::Options::SCF_MODES ScfMode>
void CalculatorBase::storeProperties() {
  // Note: The electronic energy should be stored directly in the specific calculator, i.e., DFT/HF/CC calculators.
  if (_requiredProperties.containsSubSet(Utils::Property::AtomicCharges)) {
    this->storeAtomicCharges<ScfMode>();
  }
  if (_requiredProperties.containsSubSet(Utils::Property::NAlphaElectrons) ||
      _requiredProperties.containsSubSet(Utils::Property::NBetaElectrons)) {
    this->storeNElectrons<ScfMode>();
  }
  if (_requiredProperties.containsSubSet(Utils::Property::ElectronicOccupation)) {
    this->storeElectronicOccupations<ScfMode>();
  }
  if (_requiredProperties.containsSubSet(Utils::Property::OverlapMatrix)) {
    this->storeOverlapMatrix();
  }
  if (_requiredProperties.containsSubSet(Utils::Property::DensityMatrix)) {
    this->storeDensityMatrix<ScfMode>();
  }
  if (_requiredProperties.containsSubSet(Utils::Property::Gradients) ||
      _requiredProperties.containsSubSet(Utils::Property::PointChargesGradients)) {
    this->storeGradients(ScfMode);
  }
  if (_requiredProperties.containsSubSet(Utils::Property::OneElectronMatrix)) {
    this->storeCoreHamiltonian();
  }
  if (_requiredProperties.containsSubSet(Utils::Property::AOtoAtomMapping)) {
    this->storeAOToAtomMapping();
  }
  if (_requiredProperties.containsSubSet(Utils::Property::OrbitalFragmentPopulations)) {
    this->storeOrbitalFragmentPopulations<ScfMode>();
  }
  auto atomCollection = this->getStructure();
  Scine::Utils::ResultsAutoCompleter completer(*atomCollection);
  if (_requiredProperties.containsSubSet(Utils::Property::Energy)) {
    completer.setWantedProperties(Scine::Utils::Property::Energy);
  }
  if (_requiredProperties.containsSubSet(Scine::Utils::Property::BondOrderMatrix)) {
    completer.addOneWantedProperty(Scine::Utils::Property::BondOrderMatrix);
  }
  if (_requiredProperties.containsSubSet(Scine::Utils::Property::Hessian) or
      _requiredProperties.containsSubSet(Scine::Utils::Property::Thermochemistry)) {
    completer.addOneWantedProperty(Scine::Utils::Property::Thermochemistry);
    completer.setTemperature(_settings->getDouble(Scine::Utils::SettingsNames::temperature));
    completer.setPressure(_settings->getDouble(Scine::Utils::SettingsNames::pressure));
  }
  completer.generateProperties(*_results, *atomCollection);
  _results->set<Scine::Utils::Property::SuccessfulCalculation>(true);
}

template void CalculatorBase::storeProperties<RESTRICTED>();
template void CalculatorBase::storeProperties<UNRESTRICTED>();

bool CalculatorBase::propertyRequiresSCF() {
  std::vector<Scine::Utils::Property> scfProperties = {Scine::Utils::Property::Energy,
                                                       Scine::Utils::Property::Gradients,
                                                       Scine::Utils::Property::Hessian,
                                                       Scine::Utils::Property::BondOrderMatrix,
                                                       Scine::Utils::Property::Thermochemistry,
                                                       Scine::Utils::Property::AtomicCharges,
                                                       Scine::Utils::Property::DensityMatrix,
                                                       Scine::Utils::Property::ElectronicOccupation,
                                                       Scine::Utils::Property::PointChargesGradients};
  for (const auto& prop : scfProperties) {
    if (this->_requiredProperties.containsSubSet({prop})) {
      return true;
    }
  }
  return false;
}
std::shared_ptr<Sty::SystemController> CalculatorBase::getSystemController() {
  if (!_system) {
    // Parse current settings
    auto settings = Settings();
    // throws error for wrong input and updates 'any' entries
    Utils::Solvation::ImplicitSolvation::solvationNeededAndPossible(availableSolvationModels(), *_settings);
    // Apply user settings
    _settings->applyTo(settings);
    // Apply fixed settings and those that are specific to the Calculator implementation at hand.
    this->applyFixedSettings(settings);
    // Generate a unique name
    Scine::Utils::UniqueIdentifier uid;
    settings.name = uid.getStringRepresentation();
    // Generate the system
    _system = std::make_shared<SystemController>(_geometry, settings);
  }
  return _system;
}
template<>
Utils::SpinAdaptedMatrix
CalculatorBase::orbitalPopulationsToFragmentPopulations(const SPMatrix<Sty::RESTRICTED>& orbitalPopulations) {
  const unsigned int nAtoms = orbitalPopulations.rows();
  Eigen::VectorXd ghostAtomZeros = Eigen::VectorXd::Zero(nAtoms);
  const auto& atoms = this->getSystemController()->getGeometry()->getAtoms();
  for (unsigned int iAtom = 0; iAtom < nAtoms; ++iAtom) {
    const auto& atom = atoms[iAtom];
    if (atom->isDummy()) {
      continue;
    }
    ghostAtomZeros(iAtom) = 1.0;
  }
  // Eliminate the populations in the matrix on ghost atoms by multiplying with the ghostAtomZeros vector.
  // Then, we only have to take the column-wise sum to get the populations.
  Eigen::VectorXd orbitalFragmentPopulations =
      Eigen::VectorXd((orbitalPopulations.array().colwise() * ghostAtomZeros.array()).colwise().sum().transpose());
  return Utils::SpinAdaptedMatrix::createRestricted(orbitalFragmentPopulations.transpose());
}

template<>
Utils::SpinAdaptedMatrix
CalculatorBase::orbitalPopulationsToFragmentPopulations(const SPMatrix<Sty::UNRESTRICTED>& orbitalPopulations) {
  const auto alpha = this->orbitalPopulationsToFragmentPopulations<RESTRICTED>(orbitalPopulations.alpha);
  const auto beta = this->orbitalPopulationsToFragmentPopulations<RESTRICTED>(orbitalPopulations.beta);
  return Utils::SpinAdaptedMatrix::createUnrestricted(alpha.restrictedMatrix(), beta.restrictedMatrix());
}

template<Sty::Options::SCF_MODES ScfMode>
void CalculatorBase::storeOrbitalFragmentPopulations() {
  auto systemController = this->getSystemController();
  auto orbitalWisePopulations = MullikenPopulationCalculator<ScfMode>::calculateAtomwiseOrbitalPopulations(
      systemController->getActiveOrbitalController<ScfMode>()->getCoefficients(),
      systemController->getOneElectronIntegralController()->getOverlapIntegrals(),
      systemController->getAtomCenteredBasisController()->getBasisIndices());

  auto matrix = orbitalPopulationsToFragmentPopulations<ScfMode>(orbitalWisePopulations);
  _results->set<Scine::Utils::Property::OrbitalFragmentPopulations>(matrix);
}
template void CalculatorBase::storeOrbitalFragmentPopulations<RESTRICTED>();
template void CalculatorBase::storeOrbitalFragmentPopulations<UNRESTRICTED>();
} /* namespace Serenity */
} /* namespace Scine */
