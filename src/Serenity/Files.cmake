cmake_minimum_required(VERSION 3.9)
set(SERENITY_MODULE_FILES
  "Serenity/Calculators/CalculatorBase.cpp"
  "Serenity/Calculators/CalculatorBase.h"
  "Serenity/Calculators/CCCalculator.cpp"
  "Serenity/Calculators/CCCalculator.h"
  "Serenity/Calculators/DFTCalculator.cpp"
  "Serenity/Calculators/DFTCalculator.h"
  "Serenity/Calculators/HFCalculator.cpp"
  "Serenity/Calculators/HFCalculator.h"
  "Serenity/Calculators/ScineSettings.cpp"
  "Serenity/Calculators/ScineSettings.h"
  "Serenity/Calculators/SerenityState.h"
  "Serenity/Calculators/SerenityEmbeddingCalculator.cpp"
  "Serenity/Calculators/SerenityEmbeddingCalculator.h"
  "Serenity/Calculators/SerenityEmbeddingCalculatorSettings.cpp"
  "Serenity/Calculators/SerenityEmbeddingCalculatorSettings.h"
  "Serenity/SerenityModule.cpp"
  "Serenity/SerenityModule.h"
  "Serenity/Utilities/SerenityConversionFunctions.cpp"
  "Serenity/Utilities/SerenityConversionFunctions.h"
)

set(SERENITY_TEST_FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/Tests/Resources/PathToTestResources.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/Tests/Calculators/SerenityEmbeddingCalculatorTest.cpp"
)
