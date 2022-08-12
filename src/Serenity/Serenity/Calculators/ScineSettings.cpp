/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Serenity/Calculators/ScineSettings.h"
/* Serenity Includes */
#include <settings/Settings.h>
/* Scine Includes */
#include <Utils/UniversalSettings/SettingsNames.h>
#include <algorithm>

namespace Sty = Serenity;

namespace Scine {
namespace Serenity {

using namespace Scine::Utils;

template<class T>
std::string toString(T& field) {
  std::string ret;
  Sty::Options::resolve(ret, field);
  return ret;
}

ScineSettings::ScineSettings() : Settings("SerenityDFTSettings") {
  auto defaults = Sty::Settings();

  using namespace Scine::Utils::UniversalSettings;

  // Technical
  BoolDescriptor show_serenity_output("Switch: turns Serenity text output on and off.");
  show_serenity_output.setDefaultValue(false);
  this->_fields.push_back("show_serenity_output", show_serenity_output);

  // Generalized duplicates (higher in hierarchy than the Serenity settings)
  IntDescriptor spin_multiplicity("The multiplicity.");
  spin_multiplicity.setDefaultValue(abs(defaults.spin) + 1);
  this->_fields.push_back(SettingsNames::spinMultiplicity, spin_multiplicity);

  IntDescriptor molecular_charge("The molecular charge.");
  molecular_charge.setDefaultValue(defaults.charge);
  this->_fields.push_back("molecular_charge", molecular_charge);

  IntDescriptor scf_max_iterations("The maximum number of SCF iterations.");
  scf_max_iterations.setDefaultValue(defaults.scf.maxCycles);
  this->_fields.push_back(SettingsNames::maxScfIterations, scf_max_iterations);

  DoubleDescriptor scf_energy_threshold("The energy convergence threshold.");
  scf_energy_threshold.setDefaultValue(defaults.scf.energyThreshold);
  this->_fields.push_back(SettingsNames::selfConsistenceCriterion, scf_energy_threshold);

  OptionListDescriptor spin_mode("The description of the spin mode such as restricted or unrestricted.");
  spin_mode.addOption("any");
  spin_mode.addOption("restricted");
  spin_mode.addOption("unrestricted");
  spin_mode.setDefaultOption("any");
  this->_fields.push_back(SettingsNames::spinMode, spin_mode);

  StringDescriptor method("The actual method used.");
  method.setDefaultValue("PBE");
  this->_fields.push_back(SettingsNames::method, method);

  StringDescriptor basis_set("The label of the basis set.");
  basis_set.setDefaultValue(defaults.basis.label);
  this->_fields.push_back(SettingsNames::basisSet, basis_set);

  DoubleDescriptor temperature("The temperature.");
  temperature.setDefaultValue(300.0);
  temperature.setMinimum(0.0);
  this->_fields.push_back(SettingsNames::temperature, temperature);

  DoubleDescriptor electronic_temperature("The electronic temperature.");
  electronic_temperature.setDefaultValue(0.0);
  electronic_temperature.setMinimum(0.0);
  electronic_temperature.setMaximum(0.0);
  this->_fields.push_back(SettingsNames::electronicTemperature, electronic_temperature);

  StringDescriptor solvation("The solvation method/description.");
  solvation.setDefaultValue("none");
  this->_fields.push_back(SettingsNames::solvation, solvation);

  StringDescriptor solvent("The solvent.");
  solvent.setDefaultValue("none");
  this->_fields.push_back(SettingsNames::solvent, solvent);

  // Serenity
  // - Basis - Block
  StringDescriptor basis_auxJLabel("Basis set label for the auxiliary basis for Coulomb integrals.");
  basis_auxJLabel.setDefaultValue(defaults.basis.auxJLabel);
  this->_fields.push_back("basis_auxJLabel", basis_auxJLabel);

  StringDescriptor basis_auxCLabel("Basis set label for the auxiliary basis for correlation treatments.");
  basis_auxCLabel.setDefaultValue(defaults.basis.auxCLabel);
  this->_fields.push_back("basis_auxCLabel", basis_auxCLabel);

  BoolDescriptor basis_makeSphericalBasis("Switch: use a spherical basis (or a cartesian one).");
  basis_makeSphericalBasis.setDefaultValue(defaults.basis.makeSphericalBasis);
  this->_fields.push_back("basis_makeSphericalBasis", basis_makeSphericalBasis);

  DoubleDescriptor basis_integralThreshold("The threshold for prescreening in integral evaluations.");
  basis_integralThreshold.setDefaultValue(defaults.basis.integralThreshold);
  this->_fields.push_back("basis_integralThreshold", basis_integralThreshold);

  StringDescriptor basis_basisLibPath("The path to the basis set files.");
  basis_basisLibPath.setDefaultValue(defaults.basis.basisLibPath);
  this->_fields.push_back("basis_basisLibPath", basis_basisLibPath);

  IntDescriptor basis_firstECP("The nuclear charge number of the first atom in the PSE to receive ECPs.");
  basis_firstECP.setDefaultValue(defaults.basis.firstECP);
  this->_fields.push_back("basis_firstECP", basis_firstECP);

  // - Grid - Block
  StringDescriptor grid_gridType("The identifier for the type of grid.");
  grid_gridType.setDefaultValue(toString(defaults.grid.gridType));
  this->_fields.push_back("grid_gridType", grid_gridType);

  IntDescriptor grid_smallGridAccuracy("The accuracy of the smaller integration grid used in temporary steps.");
  // grid_smallGridAccuracy.setDefaultValue(defaults.grid.smallGridAccuracy);
  grid_smallGridAccuracy.setDefaultValue(3);
  grid_smallGridAccuracy.setMinimum(1);
  grid_smallGridAccuracy.setMaximum(7);
  this->_fields.push_back("grid_smallGridAccuracy", grid_smallGridAccuracy);

  IntDescriptor grid_accuracy("The accuracy of the integration grid.");
  // grid_accuracy.setDefaultValue(defaults.grid.accuracy);
  grid_accuracy.setDefaultValue(5);
  grid_accuracy.setMinimum(1);
  grid_accuracy.setMaximum(7);
  this->_fields.push_back("grid_accuracy", grid_accuracy);

  // - SCF - Block
  StringDescriptor scf_initialguess("The initial guess to be used.");
  scf_initialguess.setDefaultValue(toString(defaults.scf.initialguess));
  this->_fields.push_back("scf_initialguess", scf_initialguess);

  IntDescriptor scf_seriesDampingInitialSteps("The number of initial dampening steps.");
  scf_seriesDampingInitialSteps.setDefaultValue(5);
  this->_fields.push_back("scf_seriesDampingInitialSteps", scf_seriesDampingInitialSteps);

  // - PCM - Block
  IntDescriptor pcm_alpha("The sharpness parameter for the molecular surface model function for DELLEY-type surfaces.");
  pcm_alpha.setDefaultValue(50);
  pcm_alpha.setMinimum(0);
  this->_fields.push_back("pcm_alpha", pcm_alpha);

  BoolDescriptor pcm_scaling("If true, the atom-radii used for the cavity construction are scaled by a factor of 1.2.");
  pcm_scaling.setDefaultValue(false);
  this->_fields.push_back("pcm_scaling", pcm_scaling);

  StringDescriptor pcm_radiiType("The atomic radii-set to be used in the cavity construction.");
  pcm_radiiType.setDefaultValue("uff");
  this->_fields.push_back("pcm_radiiType", pcm_radiiType);

  this->resetToDefaults();
}

void ScineSettings::applyTo(Sty::Settings& settings) {
  if (!this->valid()) {
    this->throwIncorrectSettings();
  }
  // tmp variable
  std::string value;

  // Mandatory
  // TODO
  settings.path = settings.path + "serenity_tmp/";

  // Serenity
  // - Basis - Block
  settings.basis.auxJLabel = this->getString("basis_auxJLabel");
  std::transform(settings.basis.auxJLabel.begin(), settings.basis.auxJLabel.end(), settings.basis.auxJLabel.begin(), ::toupper);
  settings.basis.auxCLabel = this->getString("basis_auxCLabel");
  std::transform(settings.basis.auxCLabel.begin(), settings.basis.auxCLabel.end(), settings.basis.auxCLabel.begin(), ::toupper);
  settings.basis.makeSphericalBasis = this->getBool("basis_makeSphericalBasis");
  settings.basis.integralThreshold = this->getDouble("basis_integralThreshold");
  settings.basis.basisLibPath = this->getString("basis_basisLibPath");
  settings.basis.firstECP = this->getInt("basis_firstECP");
  // - Grid - Block
  value = this->getString("grid_gridType");
  Sty::Options::resolve(value, settings.grid.gridType);
  settings.grid.smallGridAccuracy = this->getInt("grid_smallGridAccuracy");
  settings.grid.accuracy = this->getInt("grid_accuracy");
  // - SCF - Block
  value = this->getString("scf_initialguess");
  Sty::Options::resolve(value, settings.scf.initialguess);
  settings.scf.seriesDampingInitialSteps = this->getInt("scf_seriesDampingInitialSteps");
  // Spin mode
  this->resolveSpinMode();
  value = this->getString(SettingsNames::spinMode);
  Sty::Options::resolve(value, settings.scfMode);

  // - PCM - Block
  std::string solvent = this->getString(SettingsNames::solvent);
  std::string solvation = this->getString(SettingsNames::solvation);
  std::string pcm_radiiType = this->getString("pcm_radiiType");
  if (solvation.empty() || solvation == "none") {
    settings.pcm.use = false;
  }
  else {
    settings.pcm.use = true;
    Sty::Options::resolve(solvation, settings.pcm.solverType);
    Sty::Options::resolve(solvent, settings.pcm.solvent);
    Sty::Options::resolve(pcm_radiiType, settings.pcm.radiiType);
  }
  settings.pcm.alpha = this->getInt("pcm_alpha");
  settings.pcm.scaling = this->getBool("pcm_scaling");

  // Generalized duplicates
  settings.spin = (this->getInt(SettingsNames::spinMultiplicity) - 1);
  settings.charge = this->getInt("molecular_charge");
  value = this->getString(SettingsNames::spinMode);
  Sty::Options::resolve(value, settings.scfMode);
  settings.scf.energyThreshold = this->getDouble(SettingsNames::selfConsistenceCriterion);
  settings.scf.maxCycles = this->getInt(SettingsNames::maxScfIterations);
  settings.basis.label = this->getString(SettingsNames::basisSet);
  std::transform(settings.basis.label.begin(), settings.basis.label.end(), settings.basis.label.begin(), ::toupper);
}

void ScineSettings::resolveSpinMode() {
  if (this->getString(SettingsNames::spinMode) == "any") {
    int multiplicity = this->getInt(SettingsNames::spinMultiplicity);
    std::string spinMode = (multiplicity == 1) ? "restricted" : "unrestricted";
    this->modifyString(SettingsNames::spinMode, spinMode);
  }
}

} /* namespace Serenity */
} /* namespace Scine */
