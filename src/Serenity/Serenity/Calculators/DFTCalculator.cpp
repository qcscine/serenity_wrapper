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
         Scine::Utils::Property::ElectronicOccupation | Scine::Utils::Property::OneElectronMatrix |
         Scine::Utils::Property::PointChargesGradients | Scine::Utils::Property::OrbitalFragmentPopulations;
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

void DFTCalculator::storeGradients(Sty::Options::SCF_MODES ScfMode) {
  Eigen::MatrixXd gradients;
  std::shared_ptr<Eigen::MatrixXd> pointChargeGradients = nullptr;
  if (ScfMode == Sty::Options::SCF_MODES::RESTRICTED) {
    auto potBundle = _system->getElectronicStructure<Sty::Options::SCF_MODES::RESTRICTED>()->getPotentialBundle();
    gradients = potBundle->getGradients().eval();
    if (_system->hasExternalCharges()) {
      pointChargeGradients = std::make_shared<Eigen::MatrixXd>(potBundle->getPointChargeGradients());
    }
  }
  else {
    auto potBundle = _system->getElectronicStructure<Sty::Options::SCF_MODES::UNRESTRICTED>()->getPotentialBundle();
    gradients = potBundle->getGradients().eval();
    if (_system->hasExternalCharges()) {
      pointChargeGradients = std::make_shared<Eigen::MatrixXd>(potBundle->getPointChargeGradients());
    }
  }
  if (_system->getSettings().dft.dispersion != Sty::Options::DFT_DISPERSION_CORRECTIONS::NONE) {
    // Dispersion Correction components
    gradients += Sty::DispersionCorrectionCalculator::calcDispersionGradientCorrection(
        _system->getSettings().dft.dispersion, _system->getGeometry(), _system->getSettings().dft.functional);
  }
  _system->getGeometry()->setGradients(gradients);
  _results->set<Scine::Utils::Property::Gradients>(gradients);
  if (pointChargeGradients != nullptr) {
    _system->setPointChargeGradients(*pointChargeGradients);
    _results->set<Scine::Utils::Property::PointChargesGradients>(*pointChargeGradients);
  }
}

template<Sty::Options::SCF_MODES ScfMode>
void DFTCalculator::calculateImpl() {
  // Calculate energy and electronic structure
  if (this->_moved && this->propertyRequiresSCF()) {
    Sty::ScfTask<ScfMode> scf(_system);
    scf.run();
    this->_moved = false;
    this->storeElectronicEnergy<ScfMode>();
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
  this->storeProperties<ScfMode>();
}

bool DFTCalculator::supportsMethodFamily(const std::string& methodFamily) const {
  return methodFamily == "DFT";
}

} /* namespace Serenity */
} /* namespace Scine */
