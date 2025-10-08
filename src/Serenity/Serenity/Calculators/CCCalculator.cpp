/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Serenity/Calculators/CCCalculator.h"
#include "Serenity/Calculators/ScineSettings.h"
/* Serenity Includes */
#include <basis/AtomCenteredBasisController.h>
#include <data/ElectronicStructure.h>
#include <data/matrices/DensityMatrix.h>
#include <energies/EnergyContributions.h>
#include <geometry/Geometry.h>
#include <geometry/gradients/NumericalHessianCalc.h>
#include <integrals/OneElectronIntegralController.h>
#include <potentials/bundles/PotentialBundle.h>
#include <settings/Settings.h>
#include <system/SystemController.h>
#include <tasks/CoupledClusterTask.h>
#include <tasks/LocalizationTask.h>
#include <tasks/ScfTask.h>
/* Scine Includes */
#include <Utils/CalculatorBasics/ResultsAutoCompleter.h>
#include <Utils/DataStructures/AtomsOrbitalsIndexes.h>
#include <Utils/Geometry.h>
#include <Utils/Scf/LcaoUtils/ElectronicOccupation.h>
#include <Utils/Technical/UniqueIdentifier.h>
#include <Utils/Typenames.h>

namespace Sty = Serenity;

namespace Scine {
namespace Serenity {

std::string CCCalculator::name() const {
  return calculatorName;
}

Scine::Utils::PropertyList CCCalculator::possibleProperties() const {
  return Scine::Utils::Property::Energy | Scine::Utils::Property::AtomicCharges | Scine::Utils::Property::OverlapMatrix |
         Scine::Utils::Property::AOtoAtomMapping | Scine::Utils::Property::OneElectronMatrix;
}

void CCCalculator::applyFixedSettings(Sty::Settings& settings) const {
  settings.method = Sty::Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
}

void CCCalculator::storeGradients(Sty::Options::SCF_MODES ScfMode) {
  (void)ScfMode;
  throw std::runtime_error("Gradients of the energy with respect to the nuclear coordinates are not available from "
                           "Serenity for coupled cluster type methods.");
}

template<Sty::Options::SCF_MODES ScfMode>
void CCCalculator::calculateImpl() {
  if (ScfMode == Sty::Options::SCF_MODES::UNRESTRICTED)
    throw std::runtime_error("Unrestricted Coupled Cluster calculations are not yet supported in Serenity.");
  // Calculate energy and electronic structure
  Sty::Options::CC_LEVEL level = Sty::Options::CC_LEVEL::DLPNO_CCSD_T0;
  auto method = this->_settings->getString("method");
  Sty::Options::resolve(method, level);
  if (this->_moved && this->propertyRequiresSCF()) {
    Sty::ScfTask<ScfMode> scf(_system);
    scf.run();
    if (level == Sty::Options::CC_LEVEL::DLPNO_CCSD_T0 || level == Sty::Options::CC_LEVEL::CCSD_T) {
      Sty::LocalizationTask loc(_system);
      loc.settings.locType = Sty::Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO;
      loc.settings.splitValenceAndCore = true;
      loc.run();
    }
    Sty::CoupledClusterTask cc(_system);
    cc.settings.level = level;
    cc.run();
    this->_moved = false;
    _results->set<Scine::Utils::Property::SuccessfulCalculation>(true);
    auto es = _system->getElectronicStructure<ScfMode>();
    const double hf = es->getEnergy(Sty::ENERGY_CONTRIBUTIONS::HF_ENERGY);
    const double sd = es->getEnergy(Sty::ENERGY_CONTRIBUTIONS::CCSD_CORRECTION);
    double t = 0.0;
    if (level == Sty::Options::CC_LEVEL::DLPNO_CCSD_T0 || level == Sty::Options::CC_LEVEL::CCSD_T) {
      t = es->getEnergy(Sty::ENERGY_CONTRIBUTIONS::TRIPLES_CORRECTION);
    }
    const double total = hf + sd + t;
    _results->set<Scine::Utils::Property::Energy>(total);
  }
  this->storeProperties<ScfMode>();
}

bool CCCalculator::supportsMethodFamily(const std::string& methodFamily) const {
  return methodFamily == "CC";
}

} /* namespace Serenity */
} /* namespace Scine */
