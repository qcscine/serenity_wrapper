/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
/* Wrapper Includes */
#include "Serenity/Calculators/CalculatorBase.h"
#include "Serenity/Calculators/ScineSettings.h"
#include "Serenity/Calculators/SerenityState.h"
/* Serenity Includes */
#include <analysis/populationAnalysis/HirshfeldPopulationCalculator.h>
#include <analysis/populationAnalysis/MullikenPopulationCalculator.h>
#include <data/ElectronicStructure.h>
#include <data/OrbitalController.h>
#include <data/grid/BasisFunctionOnGridController.h>
#include <data/grid/BasisFunctionOnGridControllerFactory.h>
#include <data/grid/DensityMatrixDensityOnGridController.h>
#include <data/grid/DensityOnGridCalculator.h>
#include <data/matrices/DensityMatrix.h>
#include <geometry/Geometry.h>
#include <grid/GridControllerFactory.h>
#include <integrals/wrappers/Libint.h>
#include <io/FormattedOutputStream.h>
#include <math/Matrix.h>
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
  auto scine_elements = structure.getElements();
  std::vector<std::string> symbols;
  for (auto& e : scine_elements) {
    symbols.push_back(Scine::Utils::ElementInfo::symbol(e));
  }
  _geometry = std::make_shared<Geometry>(symbols, Eigen::MatrixXd(structure.getPositions()));
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
  std::vector<Scine::Utils::ElementType> scine_elements;
  for (auto& s : _geometry->getAtomSymbols()) {
    if (s.size() == 2) {
      s[1] = tolower(s[1]);
    }
    scine_elements.push_back(Scine::Utils::ElementInfo::elementTypeForSymbol(s));
  };
  Scine::Utils::PositionCollection scine(*_scinePositions);
  return std::make_unique<Scine::Utils::AtomCollection>(scine_elements, scine);
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
  //  auto old = iOOptions.printSystemInfoOnCreation;
  //  iOOptions.printSystemInfoOnCreation = false;
  _system = std::make_shared<SystemController>(_geometry, settings);
  //  iOOptions.printSystemInfoOnCreation = old;

  // Load data into the new system generated from the state
  if (castState->system->hasElectronicStructure<RESTRICTED>()) {
    CoefficientMatrix<RESTRICTED> coeff(_system->getBasisController());
    auto orig = castState->system->getElectronicStructure<RESTRICTED>()->getMolecularOrbitals()->getCoefficients();
    (Eigen::MatrixXd) coeff = (Eigen::MatrixXd)orig;
    auto eval = castState->system->getElectronicStructure<RESTRICTED>()->getMolecularOrbitals()->getEigenvalues();
    auto orbitals = std::make_shared<OrbitalController<RESTRICTED>>(_system->getBasisController());
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
    auto orbitals = std::make_shared<OrbitalController<UNRESTRICTED>>(_system->getBasisController());
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
  //  auto old = iOOptions.printSystemInfoOnCreation;
  //  iOOptions.printSystemInfoOnCreation = false;
  auto geometry = std::make_shared<Geometry>(_geometry->getAtomSymbols(), _geometry->getCoordinates());
  auto system = std::make_shared<SystemController>(geometry, settings);
  //  iOOptions.printSystemInfoOnCreation = old;

  if (_system) {
    if (_system->hasElectronicStructure<RESTRICTED>()) {
      CoefficientMatrix<RESTRICTED> coeff(system->getBasisController());
      (Eigen::MatrixXd) coeff =
          (Eigen::MatrixXd)_system->getElectronicStructure<RESTRICTED>()->getMolecularOrbitals()->getCoefficients();
      auto eval = _system->getElectronicStructure<RESTRICTED>()->getMolecularOrbitals()->getEigenvalues();
      auto orbitals = std::make_shared<OrbitalController<RESTRICTED>>(system->getBasisController());
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
      auto orbitals = std::make_shared<OrbitalController<UNRESTRICTED>>(system->getBasisController());
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
    throw std::runtime_error("Missing geometry in Serenity DFT Calculator");
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

  // System Initializations
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
    std::cout << std::endl;
    std::cout << "    Generated new Serenity system with UID:" << std::endl;
    std::cout << "        " << settings.name << std::endl;
    std::cout << std::endl;
    // Generate the system
    _system = std::make_shared<SystemController>(_geometry, settings);
  }

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

} /* namespace Serenity */
} /* namespace Scine */
