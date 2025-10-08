/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SERENITY_CONVERSIONFUNCTIONS_H
#define SERENITY_CONVERSIONFUNCTIONS_H

#include <data/matrices/SPMatrix.h>
#include <geometry/Geometry.h>

namespace Sty = Serenity;
namespace Scine {
namespace Utils {
class AtomCollection;
class MolecularOrbitals;
class SpinAdaptedMatrix;
} // namespace Utils

namespace Serenity {

/**
 * @brief Conversion functions between Serenity and SCINE objects.
 */
class SerenityConversionFunctions {
 public:
  SerenityConversionFunctions() = delete;
  ~SerenityConversionFunctions() = delete;
  /**
   * @brief Convert an atom collection to a serenity geometry.
   * @param atomCollection The atom collection.
   * @return The Serenity geometry.
   */
  static Sty::Geometry atomCollectionToGeometry(const Utils::AtomCollection& atomCollection);
  /**
   * @brief Convert a serenity geometry to an atom collection.
   * @param geometry The geometry.
   * @return The atom collection.
   */
  static Utils::AtomCollection geometryToAtomCollection(const Sty::Geometry& geometry);

  /**
   * @brief Convert a serenity coefficient matrix to a MolecularOrbitals object.
   * @param coefficients The coefficients.
   * @return The molecular orbitals.
   */
  template<Sty::Options::SCF_MODES SCFMode>
  static Utils::MolecularOrbitals molecularOrbitalsFromCoefficients(const Sty::SPMatrix<SCFMode>& coefficients);
  /**
   * @brief Convert a serenity spin polarized matrix to its utils analogue.
   * @param matrix The spin polarized matrix.
   * @return The utils spin adapted matrix.
   */
  template<Sty::Options::SCF_MODES SCFMode>
  static Utils::SpinAdaptedMatrix spMatrixToSpinAdaptedMatrix(const Sty::SPMatrix<SCFMode>& matrix);
};

} /* namespace Serenity */
} /* namespace Scine */

#endif // SERENITY_CONVERSIONFUNCTIONS_H
