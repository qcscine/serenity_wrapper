/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Serenity/Calculators/DFTCalculator.h"
#include "Serenity/Calculators/ScineSettings.h"
/* Serenity Includes */
#include <basis/AtomCenteredBasisController.h>
#include <data/ElectronicStructure.h>
#include <data/matrices/DensityMatrix.h>
#include <dft/dispersionCorrection/DispersionCorrectionCalculator.h>
#include <geometry/Geometry.h>
#include <geometry/gradients/NumericalHessianCalc.h>
#include <integrals/OneElectronIntegralController.h>
#include <misc/SerenityError.h> //Errors.
#include <potentials/bundles/PotentialBundle.h>
#include <settings/Settings.h>
#include <system/SystemController.h>
#include <tasks/ScfTask.h>
/* Scine Includes */
#include <Utils/CalculatorBasics.h>
#include <Utils/DataStructures/AtomsOrbitalsIndexes.h>
#include <Utils/Geometry.h>
#include <Utils/Scf/LcaoUtils/ElectronicOccupation.h>
#include <Utils/Technical/UniqueIdentifier.h>
#include <Utils/Typenames.h>
#include <Utils/UniversalSettings/SettingsNames.h>

namespace Sty = Serenity;

namespace Scine {
namespace Serenity {

std::string DFTCalculator::name() const {
  return "SerenityDFTCalculator";
}

Scine::Utils::PropertyList DFTCalculator::possibleProperties() const {
  return Scine::Utils::Property::Energy | Scine::Utils::Property::Gradients | Scine::Utils::Property::Hessian |
         Scine::Utils::Property::BondOrderMatrix | Scine::Utils::Property::Thermochemistry |
         Scine::Utils::Property::AtomicCharges | Scine::Utils::Property::AOtoAtomMapping |
         Scine::Utils::Property::DensityMatrix | Scine::Utils::Property::OverlapMatrix |
         Scine::Utils::Property::ElectronicOccupation;
}

void DFTCalculator::applyFixedSettings(Sty::Settings& settings) const {
  settings.method = Sty::Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  auto methodInput = Scine::Utils::CalculationRoutines::splitIntoMethodAndDispersion(
      this->_settings->getString(Scine::Utils::SettingsNames::method));
  Sty::Options::resolve(methodInput.first, settings.dft.functional);
  if (!methodInput.second.empty()) {
    Sty::Options::resolve(methodInput.second, settings.dft.dispersion);
  }
}

template<Sty::Options::SCF_MODES ScfMode>
void DFTCalculator::calculateImpl() {
  // Calculate energy and electronic structure
  if (this->_moved) {
    Sty::ScfTask<ScfMode> scf(_system);
    scf.run();
    this->_moved = false;
  }
  auto es = _system->getElectronicStructure<ScfMode>();
  _results->set<Scine::Utils::Property::Energy>(es->getEnergy());

  // Calculate gradients
  if (_requiredProperties.containsSubSet(Scine::Utils::Property::Gradients)) {
    Eigen::MatrixXd gradients;
    auto potBundle = _system->getElectronicStructure<ScfMode>()->getPotentialBundle();
    gradients = potBundle->getGradients().eval();
    if (_system->getSettings().dft.dispersion != Sty::Options::DFT_DISPERSION_CORRECTIONS::NONE) {
      // Dispersion Correction components
      gradients += Sty::DispersionCorrectionCalculator::calcDispersionGradientCorrection(
          _system->getSettings().dft.dispersion, _system->getGeometry(), _system->getSettings().dft.functional);
    }
    _system->getGeometry()->setGradients(gradients);
    _results->set<Scine::Utils::Property::Gradients>(gradients);
  }

  // Calculate Hessian
  if (_requiredProperties.containsSubSet(Scine::Utils::Property::Hessian) or
      _requiredProperties.containsSubSet(Scine::Utils::Property::Thermochemistry)) {
    // reroute output
    std::ofstream out(_system->getSettings().path + "/hessian.cout.txt");
    std::streambuf* coutbuf = std::cout.rdbuf();
    std::cout.rdbuf(out.rdbuf());
    Sty::NumericalHessianCalc<ScfMode> hessianCalc(0.0e0, 0.001, true);
    // calculate
    try {
      auto hessian = hessianCalc.calcHessian(_system);
      _results->set<Scine::Utils::Property::Hessian>(hessian);
    }
    catch (Sty::SerenityError& e) {
      throw Core::UnsuccessfulCalculationException(e.what());
    }
    // reset output
    std::cout.rdbuf(coutbuf);
  }
  _results->set<Scine::Utils::Property::SuccessfulCalculation>(true);

  if (_requiredProperties.containsSubSet(Scine::Utils::Property::AtomicCharges)) {
    auto charges = getMullikenCharges<ScfMode>();
    _results->set<Scine::Utils::Property::AtomicCharges>(charges);
  }

  /*
   * Autocompletion part of the results
   */
  auto atomCollection = this->getStructure();
  Scine::Utils::ResultsAutoCompleter completer(*atomCollection);
  // Fill results with some basic data
  //  - AO to Atom Mapping
  auto indices = _system->getAtomCenteredBasisController()->getBasisIndices();
  Scine::Utils::AtomsOrbitalsIndexes counts(indices.size());
  for (unsigned int i = 0; i < indices.size(); i++) {
    counts.addAtom(indices[i].second - indices[i].first);
  }
  _results->set<Scine::Utils::Property::AOtoAtomMapping>(counts);
  //  - AO Overlap
  _results->set<Scine::Utils::Property::OverlapMatrix>(_system->getOneElectronIntegralController()->getOverlapIntegrals());
  //  - AO Density Matrix
  auto dmat = es->getDensityMatrix();
  _results->set<Scine::Utils::Property::DensityMatrix>(this->convertDensityMatrix(dmat, _system->getNElectrons<ScfMode>()));
  //  - Occupations
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
  // Autocomplete bond orders
  completer.setWantedProperties(Scine::Utils::Property::Energy);
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
}

bool DFTCalculator::supportsMethodFamily(const std::string& methodFamily) const {
  return methodFamily == "DFT";
}

} /* namespace Serenity */
} /* namespace Scine */
