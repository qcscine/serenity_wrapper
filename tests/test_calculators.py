#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
See LICENSE.txt for details.
"""

import pytest
import scine_utilities as utils
import scine_serenity_wrapper

def create_h2() -> utils.AtomCollection:
    return utils.AtomCollection(
        [utils.Elements.H, utils.Elements.H],
        [[0.0, 0.0, 0.0], [0.7, 0.0, 0.0]]
    )

def test_dft_restricted() -> None:
    h2 = create_h2()
    module_manager = utils.core.ModuleManager.get_instance()
    calculator = module_manager.get('calculator', 'dft')
    assert calculator.name() == 'SerenityDFTCalculator'
    calculator.structure = h2
    calculator.settings['method'] = 'pbe-d3bj'
    calculator.settings['basis_set'] = 'def2-tzvp'
    results = calculator.calculate()
    assert results.successful_calculation
    assert abs(results.energy - -1.166043) < 1e-6

def test_dft_unrestricted() -> None:
    h2 = create_h2()
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

def test_hf_restricted() -> None:
    h2 = create_h2()
    module_manager = utils.core.ModuleManager.get_instance()
    calculator = module_manager.get('calculator', 'hf')
    assert calculator.name() == 'SerenityHFCalculator'
    calculator.structure = h2
    calculator.settings['method'] = 'hf'
    calculator.settings['basis_set'] = 'def2-tzvp'
    results = calculator.calculate()
    assert results.successful_calculation
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
    results = calculator.calculate()
    assert results.successful_calculation
    assert abs(results.energy - -1.132535) < 1e-6

def test_ccsd_t_restricted() -> None:
    h2 = create_h2()
    module_manager = utils.core.ModuleManager.get_instance()
    calculator = module_manager.get('calculator', 'cc')
    assert calculator.name() == 'SerenityCCCalculator'
    calculator.structure = h2
    calculator.settings['method'] = 'ccsd(t)'
    calculator.settings['basis_set'] = 'def2-tzvp'
    results = calculator.calculate()
    assert results.successful_calculation
    assert abs(results.energy - -1.168261) < 1e-6

def test_dlpno_ccsd_t0_restricted() -> None:
    h2 = create_h2()
    module_manager = utils.core.ModuleManager.get_instance()
    calculator = module_manager.get('calculator', 'cc')
    assert calculator.name() == 'SerenityCCCalculator'
    calculator.structure = h2
    calculator.settings['method'] = 'dlpno-ccsd(t0)'
    calculator.settings['basis_set'] = 'def2-tzvp'
    results = calculator.calculate()
    assert results.successful_calculation
    assert abs(results.energy - -1.168311) < 1e-6

