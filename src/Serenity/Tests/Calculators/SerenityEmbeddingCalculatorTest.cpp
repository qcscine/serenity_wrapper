/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Serenity/Calculators/SerenityEmbeddingCalculator.h"
#include "Tests/Resources/PathToTestResources.h"
#include "Utils/CalculatorBasics/CalculationRoutines.h"
#include "gmock/gmock.h"
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <iostream>

using namespace testing;

namespace Scine {
namespace Serenity {
namespace Tests {

class SerenityEmbeddingCalculatorTest : public Test {};

TEST_F(SerenityEmbeddingCalculatorTest, ExactDFTinDFT) {
  auto waterA = Utils::ChemicalFileHandler::read(pathToTestResources() + "/waterA.xyz").first;
  auto waterB = Utils::ChemicalFileHandler::read(pathToTestResources() + "/waterB.xyz").first;
  auto waterDimer = waterA + waterB;
  auto superSystemCalculator = Utils::CalculationRoutines::getCalculator("dft", "serenity");
  superSystemCalculator->setStructure(waterDimer);
  superSystemCalculator->settings().modifyString(Utils::SettingsNames::method, "PBE-D3BJ");
  superSystemCalculator->settings().modifyString(Utils::SettingsNames::basisSet, "def2-svp");
  superSystemCalculator->settings().modifyBool("show_serenity_output", true);
  superSystemCalculator->setRequiredProperties(Utils::PropertyList(Utils::Property::Energy));
  auto superSystemResults = superSystemCalculator->calculate("Super system single point");

  for (const auto& staticEmbedding : {true, false}) {
    SerenityEmbeddingCalculator calculator;
    auto activeCalculator = Utils::CalculationRoutines::getCalculator("dft", "serenity");
    auto environmentCalculator = Utils::CalculationRoutines::getCalculator("dft", "serenity");
    activeCalculator->settings().modifyString(Utils::SettingsNames::method, "PBE-D3BJ");
    activeCalculator->settings().modifyString(Utils::SettingsNames::basisSet, "def2-svp");
    environmentCalculator->settings().modifyString(Utils::SettingsNames::basisSet, "def2-svp");
    environmentCalculator->settings().modifyString(Utils::SettingsNames::method, "PBE-D3BJ");
    calculator.settings().modifyBool("show_serenity_output", true);

    if (!calculator.settings().valid()) {
      calculator.settings().throwIncorrectSettings();
    }
    calculator.settings().modifyString("non_additive_xc_functional", "PBE-D3BJ");
    calculator.settings().modifyString(Utils::SettingsNames::basisSet, "def2-svp");
    activeCalculator->setStructure(waterA);
    environmentCalculator->setStructure(waterB);
    calculator.setUnderlyingCalculators({activeCalculator, environmentCalculator});
    calculator.setRequiredProperties(Utils::PropertyList(Utils::Property::Energy));
    calculator.settings().modifyBool("static_embedding", staticEmbedding);
    auto embeddingResults = calculator.calculate("DFT-in-DFT single point");
    EXPECT_NEAR(embeddingResults.get<Utils::Property::Energy>(), superSystemResults.get<Utils::Property::Energy>(), 1e-5);
  }
}

TEST_F(SerenityEmbeddingCalculatorTest, DLPNOCCSDTinDFT) {
  SerenityEmbeddingCalculator calculator;
  auto activeCalculator = Utils::CalculationRoutines::getCalculator("cc", "serenity");
  auto environmentCalculator = Utils::CalculationRoutines::getCalculator("dft", "serenity");
  activeCalculator->settings().modifyString(Utils::SettingsNames::method, "DLPNO-CCSD(T0)");
  activeCalculator->settings().modifyString(Utils::SettingsNames::basisSet, "def2-svp");
  environmentCalculator->settings().modifyString(Utils::SettingsNames::basisSet, "def2-svp");
  environmentCalculator->settings().modifyString(Utils::SettingsNames::method, "PBE-D3BJ");
  calculator.settings().modifyBool("show_serenity_output", true);

  if (!calculator.settings().valid()) {
    calculator.settings().throwIncorrectSettings();
  }
  calculator.settings().modifyString("non_additive_xc_functional", "PBE-D3BJ");
  calculator.settings().modifyString(Utils::SettingsNames::basisSet, "def2-svp");
  auto waterA = Utils::ChemicalFileHandler::read(pathToTestResources() + "/waterA.xyz").first;
  auto waterB = Utils::ChemicalFileHandler::read(pathToTestResources() + "/waterB.xyz").first;
  activeCalculator->setStructure(waterA);
  environmentCalculator->setStructure(waterB);
  calculator.setUnderlyingCalculators({activeCalculator, environmentCalculator});
  calculator.setRequiredProperties(Utils::PropertyList(Utils::Property::Energy));
  auto embeddingResults = calculator.calculate("CC-in-DFT single point");
  EXPECT_NEAR(-152.4562739471433, embeddingResults.get<Utils::Property::Energy>(), 1e-5);
}

TEST_F(SerenityEmbeddingCalculatorTest, TightDLPNOCCSDTinDFT) {
  SerenityEmbeddingCalculator calculator;
  auto activeCalculator = Utils::CalculationRoutines::getCalculator("cc", "serenity");
  auto environmentCalculator = Utils::CalculationRoutines::getCalculator("dft", "serenity");
  activeCalculator->settings().modifyString(Utils::SettingsNames::method, "DLPNO-CCSD(T0)");
  activeCalculator->settings().modifyString("pno_settings", "tight");
  activeCalculator->settings().modifyString(Utils::SettingsNames::basisSet, "def2-svp");
  environmentCalculator->settings().modifyString(Utils::SettingsNames::basisSet, "def2-svp");
  environmentCalculator->settings().modifyString(Utils::SettingsNames::method, "PBE-D3BJ");
  calculator.settings().modifyBool("show_serenity_output", true);

  if (!calculator.settings().valid()) {
    calculator.settings().throwIncorrectSettings();
  }
  calculator.settings().modifyString("non_additive_xc_functional", "PBE-D3BJ");
  calculator.settings().modifyString(Utils::SettingsNames::basisSet, "def2-svp");
  auto waterA = Utils::ChemicalFileHandler::read(pathToTestResources() + "/waterA.xyz").first;
  auto waterB = Utils::ChemicalFileHandler::read(pathToTestResources() + "/waterB.xyz").first;
  activeCalculator->setStructure(waterA);
  environmentCalculator->setStructure(waterB);
  calculator.setUnderlyingCalculators({activeCalculator, environmentCalculator});
  calculator.setRequiredProperties(Utils::PropertyList(Utils::Property::Energy));
  auto embeddingResults = calculator.calculate("CC-in-DFT single point");
  EXPECT_NEAR(-152.45640457455195, embeddingResults.get<Utils::Property::Energy>(), 1e-5);
}

TEST_F(SerenityEmbeddingCalculatorTest, TightDLPNOMP2TinDFT) {
  SerenityEmbeddingCalculator calculator;
  auto activeCalculator = Utils::CalculationRoutines::getCalculator("cc", "serenity");
  auto environmentCalculator = Utils::CalculationRoutines::getCalculator("dft", "serenity");
  activeCalculator->settings().modifyString(Utils::SettingsNames::method, "DLPNO-MP2");
  activeCalculator->settings().modifyString("pno_settings", "tight");
  activeCalculator->settings().modifyString(Utils::SettingsNames::basisSet, "def2-svp");
  environmentCalculator->settings().modifyString(Utils::SettingsNames::basisSet, "def2-svp");
  environmentCalculator->settings().modifyString(Utils::SettingsNames::method, "PBE-D3BJ");
  calculator.settings().modifyBool("show_serenity_output", true);
  calculator.settings().modifyString(Utils::SettingsNames::basisSet, "def2-svp");

  if (!calculator.settings().valid()) {
    calculator.settings().throwIncorrectSettings();
  }
  calculator.settings().modifyString("non_additive_xc_functional", "PBE-D3BJ");
  auto waterA = Utils::ChemicalFileHandler::read(pathToTestResources() + "/waterA.xyz").first;
  auto waterB = Utils::ChemicalFileHandler::read(pathToTestResources() + "/waterB.xyz").first;
  activeCalculator->setStructure(waterA);
  environmentCalculator->setStructure(waterB);
  calculator.setUnderlyingCalculators({activeCalculator, environmentCalculator});
  calculator.setRequiredProperties(Utils::PropertyList(Utils::Property::Energy));
  auto embeddingResults = calculator.calculate("CC-in-DFT single point");
  EXPECT_NEAR(-152.44472450387781, embeddingResults.get<Utils::Property::Energy>(), 1e-5);
}

TEST_F(SerenityEmbeddingCalculatorTest, ExactDFTinDFTPartitioning) {
  SerenityEmbeddingCalculator calculator;
  auto activeCalculator = Utils::CalculationRoutines::getCalculator("dft", "serenity");
  auto environmentCalculator = Utils::CalculationRoutines::getCalculator("dft", "serenity");
  activeCalculator->settings().modifyString(Utils::SettingsNames::method, "PBE-D3BJ");
  activeCalculator->settings().modifyString(Utils::SettingsNames::basisSet, "def2-svp");
  environmentCalculator->settings().modifyString(Utils::SettingsNames::basisSet, "def2-svp");
  environmentCalculator->settings().modifyString(Utils::SettingsNames::method, "PBE-D3BJ");
  calculator.settings().modifyBool("show_serenity_output", true);
  calculator.settings().modifyString("system_partitioning", "BEST_MATCH");
  calculator.settings().modifyString(Utils::SettingsNames::basisSet, "def2-svp");
  calculator.settings().modifyString("non_additive_xc_functional", "PBE-D3BJ");
  auto waterA = Utils::ChemicalFileHandler::read(pathToTestResources() + "/waterA.xyz").first;
  auto waterB = Utils::ChemicalFileHandler::read(pathToTestResources() + "/waterB.xyz").first;
  activeCalculator->setStructure(waterA);
  environmentCalculator->setStructure(waterB);
  calculator.setUnderlyingCalculators({activeCalculator, environmentCalculator});
  calculator.setRequiredProperties(Utils::PropertyList(Utils::Property::Energy));
  auto embeddingResults = calculator.calculate("DFT-in-DFT single point");

  auto waterDimer = waterA + waterB;
  auto superSystemCalculator = Utils::CalculationRoutines::getCalculator("dft", "serenity");
  superSystemCalculator->setStructure(waterDimer);
  superSystemCalculator->settings().modifyString(Utils::SettingsNames::method, "PBE-D3BJ");
  superSystemCalculator->settings().modifyString(Utils::SettingsNames::basisSet, "def2-svp");
  superSystemCalculator->setRequiredProperties(Utils::PropertyList(Utils::Property::Energy));
  auto superSystemResults = superSystemCalculator->calculate("Super system single point");

  EXPECT_NEAR(embeddingResults.get<Utils::Property::Energy>(), superSystemResults.get<Utils::Property::Energy>(), 1e-5);
}

TEST_F(SerenityEmbeddingCalculatorTest, automatedCalculatorConstruction) {
  EXPECT_NO_THROW(Utils::CalculationRoutines::getCalculator("dft/dft", "serenity/serenity"));
  auto calculator = Utils::CalculationRoutines::getCalculator("dft/dft", "serenity/serenity");
  EXPECT_EQ(calculator->name(), "QMQM");
  auto underlyingCalculators = std::static_pointer_cast<SerenityEmbeddingCalculator>(calculator)->getUnderlyingCalculators();
  EXPECT_EQ(underlyingCalculators.size(), 2);
  EXPECT_EQ(underlyingCalculators[0]->name(), "SerenityDFTCalculator");
  EXPECT_EQ(underlyingCalculators[1]->name(), "SerenityDFTCalculator");
  EXPECT_NO_THROW(Utils::CalculationRoutines::getCalculator("cc/dft", "serenity/serenity"));
  EXPECT_NO_THROW(Utils::CalculationRoutines::getCalculator("cc|cc/dft", "serenity/serenity"));
  EXPECT_NO_THROW(Utils::CalculationRoutines::getCalculator("cc|cc/dft|dft", "serenity/serenity"));
}

TEST_F(SerenityEmbeddingCalculatorTest, QMQMexternalCharges) {
  /*
   * This test runs LMP2-in-PBE-D3(BJ) embedding for a protein fragment with 37 atoms and external charges representing
   * the remaining protein/solvent environment.
   * Because of the size of the fragment, the calculation is a bit slow and takes around 290 s on my laptop. However, it
   * is a robust check if QM/QM/MM works through the serenity wrapper.
   */
  EXPECT_NO_THROW(Utils::CalculationRoutines::getCalculator("cc/dft", "serenity/serenity"));
  auto calculator = Utils::CalculationRoutines::getCalculator("cc/dft", "serenity/serenity");
  auto underlyingCalculators = std::static_pointer_cast<SerenityEmbeddingCalculator>(calculator)->getUnderlyingCalculators();
  auto proteinFragmentA = Utils::ChemicalFileHandler::read(pathToTestResources() + "/protein-ligand-fragment-A.xyz").first;
  auto proteinFragmentB = Utils::ChemicalFileHandler::read(pathToTestResources() + "/protein-ligand-fragment-B.xyz").first;
  auto activeCalculator = underlyingCalculators[0];
  auto environmentCalculator = underlyingCalculators[1];
  activeCalculator->settings().modifyString(Utils::SettingsNames::method, "DLPNO-MP2");
  environmentCalculator->settings().modifyString(Utils::SettingsNames::method, "PBE-D3BJ");
  activeCalculator->setStructure(proteinFragmentA);
  environmentCalculator->setStructure(proteinFragmentB);
  const std::string pointChargeFilePath = pathToTestResources() + "/protein-point-charges.pc";
  calculator->settings().modifyString(Utils::SettingsNames::basisSet, "def2-svp");
  calculator->settings().modifyString("point_charges_file", pointChargeFilePath);
  calculator->settings().modifyBool("show_serenity_output", true);
  calculator->settings().modifyString("system_partitioning", "BEST_MATCH");
  calculator->setRequiredProperties(Utils::PropertyList(Utils::Property::Energy));
  calculator->settings().modifyBool("automated_charges", true);
  auto embeddingResults = calculator->calculate("LMP2-in-DFT single point");
  EXPECT_NEAR(embeddingResults.get<Utils::Property::Energy>(), -1269.339759, 1e-5);
}

TEST_F(SerenityEmbeddingCalculatorTest, StructureAssignmentsThroughIndices) {
  /*
   * Same test as above but the underlying calculators are populated by the embedding calculator.
   */
  auto calculator = Utils::CalculationRoutines::getCalculator("cc/dft", "serenity/serenity");
  auto proteinFragment = Utils::ChemicalFileHandler::read(pathToTestResources() + "/protein-ligand-fragment.xyz").first;
  calculator->settings().modifyString(Utils::SettingsNames::basisSet, "def2-svp");
  calculator->settings().modifyString(Utils::SettingsNames::method, "DLPNO-MP2/PBE-D3BJ");
  const std::string pointChargeFilePath = pathToTestResources() + "/protein-point-charges.pc";
  calculator->settings().modifyString("point_charges_file", pointChargeFilePath);
  calculator->settings().modifyBool("show_serenity_output", true);
  calculator->settings().modifyString("system_partitioning", "BEST_MATCH");
  calculator->setRequiredProperties(Utils::PropertyList(Utils::Property::Energy));
  calculator->settings().modifyBool("automated_charges", true);
  std::vector<int> environmentSystemIndices = {19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33};
  std::vector<int> activeSystemIndices = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 34, 35, 36};
  calculator->settings().modifyIntListList("qmqm_atom_indices", {activeSystemIndices, environmentSystemIndices});
  calculator->setStructure(proteinFragment);
  auto embeddingResults = calculator->calculate("LMP2-in-DFT single point");
  EXPECT_NEAR(embeddingResults.get<Utils::Property::Energy>(), -1269.339759, 1e-5);
}

TEST_F(SerenityEmbeddingCalculatorTest, StaticEmbeddingMultipleSubsystems) {
  /*
   * This test uses a static embedding approach (no embedded SCF) and multiple quantum cores. These cores are denoted
   * as MP2|MP2/DFT, where '|' separates the method for the quantum cores.
   */
  auto calculator = Utils::CalculationRoutines::getCalculator("cc|cc/dft", "serenity/serenity");
  auto proteinFragment = Utils::ChemicalFileHandler::read(pathToTestResources() + "/protein-ligand-fragment.xyz").first;
  calculator->settings().modifyString(Utils::SettingsNames::basisSet, "def2-svp");
  calculator->settings().modifyString(Utils::SettingsNames::method, "DLPNO-MP2/DLPNO-MP2/PBE-D3BJ");
  const std::string pointChargeFilePath = pathToTestResources() + "/protein-point-charges.pc";
  calculator->settings().modifyString("point_charges_file", pointChargeFilePath);
  calculator->settings().modifyBool("show_serenity_output", true);
  calculator->settings().modifyString("system_partitioning", "BEST_MATCH");
  calculator->setRequiredProperties(Utils::PropertyList(Utils::Property::Energy));
  calculator->settings().modifyBool("automated_charges", true);
  calculator->settings().modifyBool("static_embedding", true);
  std::vector<int> environmentSystemIndices = {19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33};
  std::vector<int> activeSystemIndices1 = {4, 5, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 36};
  std::vector<int> activeSystemIndices2 = {0, 1, 2, 3, 6, 7, 34, 35};
  calculator->settings().modifyIntListList("qmqm_atom_indices",
                                           {activeSystemIndices1, activeSystemIndices2, environmentSystemIndices});
  calculator->setStructure(proteinFragment);
  auto embeddingResults = calculator->calculate("LMP2-in-DFT single point");
  EXPECT_NEAR(embeddingResults.get<Utils::Property::Energy>(), -1269.315704, 1e-5);
}

TEST_F(SerenityEmbeddingCalculatorTest, StaticEmbeddingLoadSupersystemMOs) {
  /*
   * This test loads orbitals for the supersystem from an external source.
   */
  auto calculator = Utils::CalculationRoutines::getCalculator("cc|cc/dft", "serenity/serenity");
  auto proteinFragment = Utils::ChemicalFileHandler::read(pathToTestResources() + "/protein-ligand-fragment.xyz").first;
  calculator->settings().modifyString(Utils::SettingsNames::basisSet, "def2-svp");
  calculator->settings().modifyString(Utils::SettingsNames::method, "DLPNO-MP2/DLPNO-MP2/PBE-D3BJ");
  const std::string pointChargeFilePath = pathToTestResources() + "/protein-point-charges.pc";
  calculator->settings().modifyString("point_charges_file", pointChargeFilePath);
  calculator->settings().modifyBool("show_serenity_output", true);
  calculator->settings().modifyString("system_partitioning", "BEST_MATCH");
  calculator->setRequiredProperties(Utils::PropertyList(Utils::Property::Energy));
  calculator->settings().modifyBool("automated_charges", true);
  calculator->settings().modifyBool("static_embedding", true);
  calculator->settings().modifyString("external_supersystem_mo_file", pathToTestResources());
  std::vector<int> environmentSystemIndices = {19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33};
  std::vector<int> activeSystemIndices1 = {4, 5, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 36};
  std::vector<int> activeSystemIndices2 = {0, 1, 2, 3, 6, 7, 34, 35};
  calculator->settings().modifyIntListList("qmqm_atom_indices",
                                           {activeSystemIndices1, activeSystemIndices2, environmentSystemIndices});
  calculator->setStructure(proteinFragment);
  auto embeddingResults = calculator->calculate("LMP2-in-DFT single point");
  EXPECT_NEAR(embeddingResults.get<Utils::Property::Energy>(), -1269.315704, 1e-5);
}

TEST_F(SerenityEmbeddingCalculatorTest, PrepareCAS) {
  auto calculator = Utils::CalculationRoutines::getCalculator("hf|hf/dft", "serenity/serenity");
  auto proteinFragment = Utils::ChemicalFileHandler::read(pathToTestResources() + "/protein-ligand-fragment.xyz").first;
  calculator->settings().modifyString(Utils::SettingsNames::basisSet, "def2-svp");
  calculator->settings().modifyString(Utils::SettingsNames::method, "HF/HF/PBE-D3BJ");
  const std::string pointChargeFilePath = pathToTestResources() + "/protein-point-charges.pc";
  calculator->settings().modifyString("point_charges_file", pointChargeFilePath);
  calculator->settings().modifyBool("show_serenity_output", true);
  calculator->settings().modifyString("system_partitioning", "BEST_MATCH");
  calculator->setRequiredProperties(Utils::PropertyList(Utils::Property::Energy));
  calculator->settings().modifyBool("automated_charges", true);
  calculator->settings().modifyBool("static_embedding", true);
  std::vector<int> environmentSystemIndices = {19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33};
  std::vector<int> activeSystemIndices1 = {4, 5, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 36};
  std::vector<int> activeSystemIndices2 = {0, 1, 2, 3, 6, 7, 34, 35};
  calculator->settings().modifyIntListList("qmqm_atom_indices",
                                           {activeSystemIndices1, activeSystemIndices2, environmentSystemIndices});
  calculator->setStructure(proteinFragment);
  calculator->settings().modifyIntList("cas_systems", {0, 1});
  auto embeddingResults = calculator->calculate("CAS preparation HF-in-DFT");

  auto qmQMCalculator = std::dynamic_pointer_cast<Scine::Core::EmbeddingCalculator>(calculator);
  auto qmCalculators = qmQMCalculator->getUnderlyingCalculators();

  std::vector<std::vector<double>> partialEnergyReferences = {{-1267.4771080647229, -1195.455652629103, -72.021455435619828},
                                                              {-1267.4771080647229, -314.35906054729764, -953.11804751742523}};
  const unsigned int nBasisFunctions = 360;
  unsigned int counter = 0;
  for (const auto& activeCalculator : qmCalculators) {
    if (activeCalculator->name() != "SerenityHFCalculator") {
      continue;
    }
    const auto& results = activeCalculator->results();
    EXPECT_TRUE(results.has<Utils::Property::OneElectronMatrix>());
    EXPECT_TRUE(results.has<Utils::Property::CoefficientMatrix>());
    EXPECT_TRUE(results.has<Utils::Property::PartialEnergies>());

    auto hCore = results.get<Utils::Property::OneElectronMatrix>();
    auto valenceCoefficients = results.get<Utils::Property::CoefficientMatrix>().restrictedMatrix();
    auto partialEnergies = results.get<Utils::Property::PartialEnergies>();
    EXPECT_EQ(hCore.rows(), hCore.cols());
    EXPECT_EQ(hCore.rows(), nBasisFunctions);
    EXPECT_EQ(valenceCoefficients.rows(), nBasisFunctions);
    EXPECT_EQ(valenceCoefficients.cols(), nBasisFunctions);

    EXPECT_TRUE(partialEnergies.find("total_uncorrelated_energy") != partialEnergies.end());
    EXPECT_TRUE(partialEnergies.find("uncorrelated_active_space_energy") != partialEnergies.end());
    EXPECT_TRUE(partialEnergies.find("total_core_energy") != partialEnergies.end());
    EXPECT_NEAR(partialEnergies["total_uncorrelated_energy"], partialEnergyReferences[counter][0], 1e-2);
    EXPECT_NEAR(partialEnergies["uncorrelated_active_space_energy"], partialEnergyReferences[counter][1], 1e-3);
    EXPECT_NEAR(partialEnergies["total_core_energy"], partialEnergyReferences[counter][2], 1e-2);
    counter++;
  }
}

} /* namespace Tests */
} /* namespace Serenity */
} /* namespace Scine */
