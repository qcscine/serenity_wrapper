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

            import scine_utilities as utils

            h2 = utils.AtomCollection(
                [utils.ElementType.H, utils.ElementType.H],
                [[-0.7, 0.0, 0.0], [0.7, 0.0, 0.0]]
            )

            # Test R-KS-DFT
            module_manager = utils.core.ModuleManager.get_instance()
            calculator = module_manager.get('calculator', 'dft')
            assert calculator.name() == 'SerenityDFTCalculator'
            calculator.structure = h2
            calculator.settings['method'] = 'pbe-d3bj'
            calculator.settings['basis_set'] = 'def2-tzvp'
            results = calculator.calculate()
            assert results.successful_calculation
            assert abs(results.energy - -1.166043) < 1e-6
            self.output.info("R-KS-DFT worked")

            # Test U-KS-DFT
            module_manager = utils.core.ModuleManager.get_instance()
            calculator = module_manager.get('calculator', 'dft')
            assert calculator.name() == 'SerenityDFTCalculator'
            calculator.structure = h2
            calculator.settings['spin_mode'] = 'unrestricted'
            calculator.settings['method'] = 'pbe-d3bj'
            calculator.settings['basis_set'] = 'def2-tzvp'
            results = calculator.calculate()
            assert results.successful_calculation
            assert abs(results.energy - -1.166043) < 1e-6
            self.output.info("U-KS-DFT worked")

            # Test R-HF
            module_manager = utils.core.ModuleManager.get_instance()
            calculator = module_manager.get('calculator', 'hf')
            assert calculator.name() == 'SerenityHFCalculator'
            calculator.structure = h2
            calculator.settings['method'] = 'hf'
            calculator.settings['basis_set'] = 'def2-tzvp'
            results = calculator.calculate()
            assert results.successful_calculation
            assert abs(results.energy - -1.132535) < 1e-6
            self.output.info("U-KS-DFT worked")

            # Test U-HF
            module_manager = utils.core.ModuleManager.get_instance()
            calculator = module_manager.get('calculator', 'hf')
            assert calculator.name() == 'SerenityHFCalculator'
            calculator.structure = h2
            calculator.settings['spin_mode'] = 'unrestricted'
            calculator.settings['method'] = 'hf'
            calculator.settings['basis_set'] = 'def2-tzvp'
            results = calculator.calculate()
            assert results.successful_calculation
            assert abs(results.energy - -1.132535) < 1e-6
            self.output.info("U-KS-DFT worked")

            # Test R-CCSD(T)
            module_manager = utils.core.ModuleManager.get_instance()
            calculator = module_manager.get('calculator', 'cc')
            assert calculator.name() == 'SerenityCCCalculator'
            calculator.structure = h2
            calculator.settings['method'] = 'ccsd(t)'
            calculator.settings['basis_set'] = 'def2-tzvp'
            results = calculator.calculate()
            assert results.successful_calculation
            assert abs(results.energy - -1.168261) < 1e-6
            self.output.info("R-CCSD(T) worked")

            # Test R-DLPNO-CCSD(T0)
            module_manager = utils.core.ModuleManager.get_instance()
            calculator = module_manager.get('calculator', 'cc')
            assert calculator.name() == 'SerenityCCCalculator'
            calculator.structure = h2
            calculator.settings['method'] = 'dlpno-ccsd(t0)'
            calculator.settings['basis_set'] = 'def2-tzvp'
            results = calculator.calculate()
            assert results.successful_calculation
            assert abs(results.energy - -1.168311) < 1e-6
            self.output.info("R-DLPNO-CCSD(T) worked")

