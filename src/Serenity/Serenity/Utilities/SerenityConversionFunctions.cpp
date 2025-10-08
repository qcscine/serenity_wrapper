/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Serenity/Utilities/SerenityConversionFunctions.h"
#include "Core/BaseClasses/ObjectWithStructure.h"
#include "Utils/DataStructures/MolecularOrbitals.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Geometry/ElementInfo.h"
#include <math/Matrix.h>

namespace Scine {
namespace Serenity {

Sty::Geometry SerenityConversionFunctions::atomCollectionToGeometry(const Utils::AtomCollection& atomCollection) {
  auto scine_elements = atomCollection.getElements();
  std::vector<std::string> symbols;
  for (auto& e : scine_elements) {
    symbols.push_back(Scine::Utils::ElementInfo::symbol(e));
  }
  return Sty::Geometry(symbols, Eigen::MatrixXd(atomCollection.getPositions()));
}
Utils::AtomCollection SerenityConversionFunctions::geometryToAtomCollection(const Sty::Geometry& geometry) {
  std::vector<Scine::Utils::ElementType> scine_elements;
  for (auto& s : geometry.getAtomSymbols()) {
    if (s.size() == 2) {
      s[1] = tolower(s[1]);
    }
    scine_elements.push_back(Scine::Utils::ElementInfo::elementTypeForSymbol(s));
  }
  Scine::Utils::PositionCollection positionCollection(geometry.getCoordinates());
  return {scine_elements, positionCollection};
}
template<>
Utils::SpinAdaptedMatrix SerenityConversionFunctions::spMatrixToSpinAdaptedMatrix(const Sty::SPMatrix<Sty::RESTRICTED>& matrix) {
  return Utils::SpinAdaptedMatrix::createRestricted(Eigen::MatrixXd(matrix));
}
template<>
Utils::SpinAdaptedMatrix SerenityConversionFunctions::spMatrixToSpinAdaptedMatrix(const Sty::SPMatrix<Sty::UNRESTRICTED>& matrix) {
  return Utils::SpinAdaptedMatrix::createUnrestricted(matrix.alpha, matrix.beta);
}
template<>
Utils::MolecularOrbitals
SerenityConversionFunctions::molecularOrbitalsFromCoefficients(const Sty::SPMatrix<Sty::RESTRICTED>& coefficients) {
  return Utils::MolecularOrbitals::createFromRestrictedCoefficients(coefficients);
}
template<>
Utils::MolecularOrbitals
SerenityConversionFunctions::molecularOrbitalsFromCoefficients(const Sty::SPMatrix<Sty::UNRESTRICTED>& coefficients) {
  return Utils::MolecularOrbitals::createFromUnrestrictedCoefficients(coefficients.alpha, coefficients.beta);
}

} /* namespace Serenity */
} /* namespace Scine */
