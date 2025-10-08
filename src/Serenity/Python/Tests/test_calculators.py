__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import pytest
import os
import scine_utilities as utils
import numpy as np

def create_h2() -> utils.AtomCollection:
    return utils.AtomCollection(
        [utils.ElementType.H, utils.ElementType.H],
        [[-0.7, 0.0, 0.0], [0.7, 0.0, 0.0]]
    )

def test_dft_restricted() -> None:
    print("Running test:", "test_dft_restricted")
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
    print("Running test:", "test_dft_restricted_other_properties")
    h2 = create_h2()
    module_manager = utils.core.ModuleManager.get_instance()
    calculator = module_manager.get('calculator', 'dft')
    assert calculator.name() == 'SerenityDFTCalculator'
    calculator.structure = h2
    calculator.settings['method'] = 'pbe-d3bj'
    calculator.settings['basis_set'] = 'def2-tzvp'
    calculator.settings['show_serenity_output'] = True
    calculator.set_required_properties([utils.Property.AOtoAtomMapping,
                                        utils.Property.AtomicCharges,
                                        utils.Property.OneElectronMatrix,
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
    assert results.one_electron_matrix is not None

def test_dft_restricted_non_scf_properties() -> None:
    print("Running test:", "test_dft_restricted_non_scf_properties")
    h2 = create_h2()
    module_manager = utils.core.ModuleManager.get_instance()
    calculator = module_manager.get('calculator', 'dft')
    assert calculator.name() == 'SerenityDFTCalculator'
    calculator.structure = h2
    calculator.settings['method'] = 'pbe-d3bj'
    calculator.settings['basis_set'] = 'def2-tzvp'
    calculator.set_required_properties([utils.Property.AOtoAtomMapping,
                                        utils.Property.OneElectronMatrix,
                                        utils.Property.OverlapMatrix,
                                        ])
    results = calculator.calculate()
    assert results.successful_calculation
    assert results.ao_to_atom_mapping is not None
    assert results.overlap_matrix is not None
    assert results.one_electron_matrix is not None
    assert results.energy is None

def test_dft_unrestricted() -> None:
    print("Running test:", "test_dft_unrestricted")
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
    print("Running test:", "test_hf_restricted")
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
    print("Running test:", "test_hf_unrestricted")
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
    print("Running test:", "test_ccsd_t_restricted")
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
    print("Running test:", "test_dlpno_ccsd_t0_restricted")
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

def test_electrostatic_embedding() -> None:
    print("Running test:", "test_electrostatic_embedding")
    from scine_serenity_wrapper.Tests.resources import resource_path
    atom_collection, _ = utils.io.read(os.path.join(resource_path(), "protein-ligand-fragment.xyz"))
    module_manager = utils.core.ModuleManager.get_instance()
    calculator = module_manager.get('calculator', 'hf')
    assert calculator.name() == 'SerenityHFCalculator'
    calculator.structure = atom_collection
    calculator.settings['method'] = 'hf'
    calculator.settings['basis_set'] = 'def2-svp'
    calculator.settings['point_charges_file'] = os.path.join(resource_path(), "protein-point-charges.pc")
    calculator.set_required_properties([utils.Property.OneElectronMatrix])
    results = calculator.calculate()
    assert results.energy is None  # Serenity should not have run an SCF calculation.
    assert results.one_electron_matrix is not None
    some_reference_values = [-4.481078953931e+01,   1.481350968171e+01,   6.863725168208e+00]
    for i in range(3):
        ref = some_reference_values[i]
        assert abs(results.one_electron_matrix[0, i] - ref) < 1e-6
        assert abs(results.one_electron_matrix[i, 0] - ref) < 1e-6

def test_electrostatic_embedding_open_shell() -> None:
    print("Running test:", "test_electrostatic_embedding_open_shell")
    from scine_serenity_wrapper.Tests.resources import resource_path
    atom_collection, _ = utils.io.read(os.path.join(resource_path(), "protein-ligand-fragment.xyz"))
    module_manager = utils.core.ModuleManager.get_instance()
    calculator = module_manager.get('calculator', 'hf')
    assert calculator.name() == 'SerenityHFCalculator'
    calculator.structure = atom_collection
    calculator.settings['method'] = 'hf'
    calculator.settings['spin_mode'] = 'unrestricted'
    calculator.settings['molecular_charge'] = +1
    calculator.settings['spin_multiplicity'] = 2
    calculator.settings['basis_set'] = 'def2-svp'
    calculator.settings['point_charges_file'] = os.path.join(resource_path(), "protein-point-charges.pc")
    calculator.set_required_properties([utils.Property.OneElectronMatrix])
    results = calculator.calculate()
    assert results.energy is None  # Serenity should not have run an SCF calculation.
    assert results.one_electron_matrix is not None

def test_point_charge_gradients() -> None:
    print("Running test:", "test_point_charge_gradients")
    from scine_serenity_wrapper.Tests.resources import resource_path
    atom_collection, _ = utils.io.read(os.path.join(resource_path(), "etoh.xyz"))
    module_manager = utils.core.ModuleManager.get_instance()
    calculator = module_manager.get('calculator', 'hf')
    assert calculator.name() == 'SerenityHFCalculator'
    calculator.structure = atom_collection
    calculator.settings['method'] = 'hf'
    calculator.settings['basis_set'] = 'def2-svp'
    calculator.settings['point_charges_file'] = os.path.join(resource_path(), "etoh-water-charges.pc")
    calculator.set_required_properties([utils.Property.Energy, utils.Property.Gradients, utils.Property.PointChargesGradients])

    reference_gradients = np.asarray([
        [0.00918164, -0.00930644, 0.00432299],
        [-0.0104084, 0.0270658, -0.00271828],
        [-0.00519274, -0.00452777, 0.00452071],
        [-0.00225914, -0.00382822, -0.000416061],
        [-0.00231768, 0.00308402, -0.00529075],
        [0.00163825, -0.00353335, 0.0278773],
        [0.00215748, -0.000452374, -0.00555908],
        [-0.00920128, -0.00455667, 0.00422912],
        [0.0164212, -0.00409334, -0.027123]
    ])
    reference_point_charge_gradients = np.asarray([
      [ 1.54667e-06,  2.28381e-06,  3.23415e-06],
      [-7.55103e-07 ,-9.68785e-07, -1.51259e-06],
      [-8.73348e-07 ,-1.13626e-06, -1.76019e-06],
      [-6.80872e-06 ,  1.2159e-08, -6.83198e-07],
      [ 3.37468e-06 , 1.05529e-07,  2.11097e-07],
      [3.10415e-06  ,-3.4387e-08 , 3.34131e-07],
      [-2.39119e-06 , 3.21388e-06, -2.39693e-07],
      [1.292e-06, -1.59313e-06,  1.38987e-07],
      [1.09579e-06, -1.47915e-06,  1.03665e-07],
      [-9.37566e-06,    2.414e-05, -1.02289e-05],
      [5.54925e-06 ,-1.19456e-05 , 4.75587e-06],
      [3.69042e-06 ,-1.17051e-05 , 3.87317e-06]
    ])

    results = calculator.calculate()
    assert results.energy is not None
    assert results.gradients is not None
    assert results.point_charges_gradients is not None
    assert np.amax(np.abs(results.gradients - reference_gradients)) < 5e-6
    assert np.amax(np.abs(results.point_charges_gradients[:12,:] - reference_point_charge_gradients)) < 1e-8


def run_all_tests() -> None:
    test_dft_restricted()
    test_dft_unrestricted()
    test_dft_restricted_non_scf_properties()
    test_dft_restricted_other_properties()
    test_hf_restricted()
    test_hf_unrestricted()
    test_ccsd_t_restricted()
    test_dlpno_ccsd_t0_restricted()
    test_point_charge_gradients()