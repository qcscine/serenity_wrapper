__copyright__ = """This file is part of SCINE Utilities.
This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

from conans import ConanFile


class TestPackageConan(ConanFile):
    def build(self):
        pass

    def test(self):
        if self.options["scine_serenity_wrapper"].python:
            self.output.info("Trying to import 'scine_serenity_wrapper'")
            import scine_serenity_wrapper
            self.output.info("Import worked")

            from scine_serenity_wrapper.Tests.test_calculators import run_all_tests
            run_all_tests()

