/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "SerenityEmbeddingCalculatorSettings.h"
#include "Utils/UniversalSettings/SettingsNames.h"

namespace Scine {
namespace Serenity {

SerenityEmbeddingCalculatorSettings::SerenityEmbeddingCalculatorSettings() : ScineSettings() {
  Utils::UniversalSettings::StringDescriptor embeddingMethod("The embedding methodology.");
  embeddingMethod.setDefaultValue("FERMI");
  this->_fields.push_back("embedding_method", embeddingMethod);

  Utils::UniversalSettings::StringDescriptor naddXCFunctional("The non-additive exchange correlation functional.");
  naddXCFunctional.setDefaultValue("PBE-D3BJ");
  this->_fields.push_back(Utils::SettingsNames::nonAdditiveXCFunctional, naddXCFunctional);

  Utils::UniversalSettings::DoubleDescriptor fermiShift("The shift for Huzinaga type embedding.");
  fermiShift.setDefaultValue(1.0);
  this->_fields.push_back("fermi_shift", fermiShift);

  Utils::UniversalSettings::StringDescriptor localizationType("The orbital localization approach.");
  localizationType.setDefaultValue("IBO");
  this->_fields.push_back("loc_type", localizationType);

  Utils::UniversalSettings::StringDescriptor systemPartitioning("The system partitioning algorithm.");
  systemPartitioning.setDefaultValue("SPADE");
  this->_fields.push_back("system_partitioning", systemPartitioning);

  Utils::UniversalSettings::DoubleDescriptor orbitalPopulationThreshold(
      "The minimum orbital population used in the system partitioning.");
  orbitalPopulationThreshold.setDefaultValue(0.6);
  this->_fields.push_back("orbital_threshold", orbitalPopulationThreshold);

  Utils::UniversalSettings::StringDescriptor truncAlgorithm("The basis set truncation algorithm.");
  truncAlgorithm.setDefaultValue("None");
  this->_fields.push_back("trunc_algorithm", truncAlgorithm);

  Utils::UniversalSettings::DoubleDescriptor netThreshold(
      "The Mulliken net-population threshold for the basis set truncation.");
  netThreshold.setDefaultValue(1e-4);
  this->_fields.push_back("net_threshold", netThreshold);

  Utils::UniversalSettings::BoolDescriptor automatedChargesAndSpin(
      "If true, charges and spin multiplicities for the underlying"
      " calculators are determined automatically based on the supersystem charge and multiplicty."
      " Note that this will only select charges that could in principle be possible but may be unphysically"
      " for the subsystems. Therefore, this settings should only be used if the subsystem charges are determined"
      " automatically by the algorithm at a later stage. By default, true.");
  automatedChargesAndSpin.setDefaultValue(true);
  this->_fields.push_back("automated_charges", automatedChargesAndSpin);

  Utils::UniversalSettings::IntListListDescriptor qmqmAtomIndices(
      "The atom indices for each subsystem as a list of integer lists.\n"
      "If this list is empty, the subsystem partitioning must be set directly through\n"
      "the underlying calculators.");
  qmqmAtomIndices.setDefaultValue({});
  this->_fields.push_back(Utils::SettingsNames::qmqmAtomIndices, qmqmAtomIndices);

  Utils::UniversalSettings::BoolDescriptor staticEmbedding(
      "If true, no additional SCF is run for the subsystems after partitioning.");
  staticEmbedding.setDefaultValue(false);
  this->_fields.push_back("static_embedding", staticEmbedding);

  Utils::UniversalSettings::IntListDescriptor casSystems(
      "Calculator/system indices for which CAS integrals, coefficients, and energies are generated.\n"
      "See prepare_cas for more information. This means that the AO one particle\n"
      "integrals will correspond to the one particle integrals in the CAS. Furthermore,\n"
      "the partial energies of the uncorrelated CAS, the inactive orbitals, and the total\n"
      "uncorrelated energy will be written to the property partial_energies.");
  casSystems.setDefaultValue({});
  this->_fields.push_back("cas_systems", casSystems);

  Utils::UniversalSettings::StringDescriptor supersystemPath(
      "If given this path is used to load the supersystem. The path"
      "must point to a valid serenity system directory.");
  supersystemPath.setDefaultValue("");
  this->_fields.push_back("supersystem_path", supersystemPath);

  Utils::UniversalSettings::StringDescriptor externalSupersystemMOFile(
      "If given this path is used to load the supersystem MOs.");
  externalSupersystemMOFile.setDefaultValue("");
  this->_fields.push_back("external_supersystem_mo_file", externalSupersystemMOFile);

  Utils::UniversalSettings::OptionListDescriptor externalSupersystemMOType("The type of the external MOs to load.");
  externalSupersystemMOType.addOption("serenity");
  externalSupersystemMOType.addOption("turbomole");
  externalSupersystemMOType.addOption("molcas");
  externalSupersystemMOFile.setDefaultValue("serenity");
  this->_fields.push_back("external_supersystem_mo_type", externalSupersystemMOType);

  this->resetToDefaults();
}

} /* namespace Serenity */
} /* namespace Scine */
