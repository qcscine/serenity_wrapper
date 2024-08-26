__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import pytest
import scine_utilities as utils

def create_h2() -> utils.AtomCollection:
    return utils.AtomCollection(
        [utils.ElementType.H, utils.ElementType.H],
        [[-0.7, 0.0, 0.0], [0.7, 0.0, 0.0]]
    )

def test_dft_restricted() -> None:
    h2 = create_h2()
    module_manager = utils.core.ModuleManager.get_instance()
    calculator = module_manager.get('calculator', 'dft')
    assert calculator.name() == 'SerenityDFTCalculator'
    calculator.structure = h2
    calculator.settings['method'] = 'pbe-d3bj'
    calculator.settings['basis_set'] = 'def2-tzvp'
    calculator.set_required_properties([utils.Property.Energy])
    results = calculator.calculate()
    assert results.successful_calculation
    assert results.energy
    assert abs(results.energy - -1.166043) < 1e-6

def test_dft_restricted_other_properties() -> None:
    h2 = create_h2()
    module_manager = utils.core.ModuleManager.get_instance()
    calculator = module_manager.get('calculator', 'dft')
    assert calculator.name() == 'SerenityDFTCalculator'
    calculator.structure = h2
    calculator.settings['method'] = 'pbe-d3bj'
    calculator.settings['basis_set'] = 'def2-tzvp'
    calculator.set_required_properties([utils.Property.AOtoAtomMapping,
                                        utils.Property.AtomicCharges,
                                        # utils.Property.OneElectronMatrix,
                                        utils.Property.OverlapMatrix,
                                        utils.Property.Thermochemistry,
                                        utils.Property.Gradients])
    results = calculator.calculate()
    assert results.successful_calculation
    assert results.energy is not None
    assert abs(results.energy - -1.166043) < 1e-6
    assert results.gradients is not None
    assert results.ao_to_atom_mapping is not None
    assert results.atomic_charges is not None
    assert results.overlap_matrix is not None
    # assert results.one_electron_matrix is not None

# TODO This should be included again as soon as the serenity wrapper avoids running SCF calculations for everything.
# def test_dft_restricted_non_scf_properties() -> None:
#     h2 = create_h2()
#     module_manager = utils.core.ModuleManager.get_instance()
#     calculator = module_manager.get('calculator', 'dft')
#     assert calculator.name() == 'SerenityDFTCalculator'
#     calculator.structure = h2
#     calculator.settings['method'] = 'pbe-d3bj'
#     calculator.settings['basis_set'] = 'def2-tzvp'
#     calculator.set_required_properties([utils.Property.AOtoAtomMapping,
#                                         utils.Property.OneElectronMatrix,
#                                         utils.Property.OverlapMatrix,
#                                         ])
#     results = calculator.calculate()
#     assert results.successful_calculation
#     assert results.ao_to_atom_mapping is not None
#     assert results.overlap_matrix is not None
#     assert results.one_electron_matrix is not None
#     assert results.energy is None

def test_dft_unrestricted() -> None:
    h2 = create_h2()
    module_manager = utils.core.ModuleManager.get_instance()
    calculator = module_manager.get('calculator', 'dft')
    assert calculator.name() == 'SerenityDFTCalculator'
    calculator.structure = h2
    calculator.settings['spin_mode'] = 'unrestricted'
    calculator.settings['method'] = 'pbe-d3bj'
    calculator.settings['basis_set'] = 'def2-tzvp'
    calculator.set_required_properties([utils.Property.Energy])
    results = calculator.calculate()
    assert results.successful_calculation
    assert results.energy
    assert abs(results.energy - -1.166043) < 1e-6

def test_hf_restricted() -> None:
    h2 = create_h2()
    module_manager = utils.core.ModuleManager.get_instance()
    calculator = module_manager.get('calculator', 'hf')
    assert calculator.name() == 'SerenityHFCalculator'
    calculator.structure = h2
    calculator.settings['method'] = 'hf'
    calculator.settings['basis_set'] = 'def2-tzvp'
    calculator.set_required_properties([utils.Property.Energy])
    results = calculator.calculate()
    assert results.successful_calculation
    assert results.energy
    assert abs(results.energy - -1.132535) < 1e-6

def test_hf_unrestricted() -> None:
    h2 = create_h2()
    module_manager = utils.core.ModuleManager.get_instance()
    calculator = module_manager.get('calculator', 'hf')
    assert calculator.name() == 'SerenityHFCalculator'
    calculator.structure = h2
    calculator.settings['spin_mode'] = 'unrestricted'
    calculator.settings['method'] = 'hf'
    calculator.settings['basis_set'] = 'def2-tzvp'
    calculator.set_required_properties([utils.Property.Energy])
    results = calculator.calculate()
    assert results.successful_calculation
    assert results.energy
    assert abs(results.energy - -1.132535) < 1e-6

def test_ccsd_t_restricted() -> None:
    h2 = create_h2()
    module_manager = utils.core.ModuleManager.get_instance()
    calculator = module_manager.get('calculator', 'cc')
    assert calculator.name() == 'SerenityCCCalculator'
    calculator.structure = h2
    calculator.settings['method'] = 'ccsd(t)'
    calculator.settings['basis_set'] = 'def2-tzvp'
    calculator.set_required_properties([utils.Property.Energy])
    results = calculator.calculate()
    assert results.successful_calculation
    assert results.energy
    assert abs(results.energy - -1.168261) < 1e-6

def test_dlpno_ccsd_t0_restricted() -> None:
    h2 = create_h2()
    module_manager = utils.core.ModuleManager.get_instance()
    calculator = module_manager.get('calculator', 'cc')
    assert calculator.name() == 'SerenityCCCalculator'
    calculator.structure = h2
    calculator.settings['method'] = 'dlpno-ccsd(t0)'
    calculator.settings['basis_set'] = 'def2-tzvp'
    calculator.set_required_properties([utils.Property.Energy])
    results = calculator.calculate()
    assert results.successful_calculation
    assert results.energy
    assert abs(results.energy - -1.168311) < 1e-6

def run_all_tests() -> None:
    test_dft_restricted()
    test_dft_unrestricted()
    test_dft_restricted_other_properties()
    test_hf_restricted()
    test_hf_unrestricted()
    test_ccsd_t_restricted()
    test_dlpno_ccsd_t0_restricted()