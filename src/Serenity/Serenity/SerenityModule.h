/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SERENITY_SERENITYMODULE_H_
#define SERENITY_SERENITYMODULE_H_

/* External Includes */
#include <Core/Module.h>
#include <boost/dll/alias.hpp>
#include <memory>

namespace Scine {
namespace Serenity {
/**
 * @brief The SCINE Module implementation for Serenity.
 */
class SerenityModule : public Scine::Core::Module {
 public:
  std::string name() const noexcept final;

  boost::any get(const std::string& interface, const std::string& model) const final;

  bool has(const std::string& interface, const std::string& model) const noexcept final;

  std::vector<std::string> announceInterfaces() const noexcept final;

  std::vector<std::string> announceModels(const std::string& concept) const noexcept final;

  static std::shared_ptr<Module> make();
};

} /* namespace Serenity */
} /* namespace Scine */

#endif /* SERENITY_SERENITYMODULE_H_ */
