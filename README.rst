SCINE - Serenity Wrapper -- A Wrapper for Serenity
==================================================

Introduction
------------

SCINE: Serenity Wrapper is a wrapper around
`Serenity <https://github.com/qcserenity/serenity>`_; it exports the following
methods:

- Density Functional Theory (DFT)
- Hartree-Fock (HF)
- (Local) Coupled Cluster (CC, DLPNO-CC)

into the SCINE tool chain.
Each method is represented by its own ``Calculator`` and the entire wrapper
constitutes a SCINE module that can be loaded dynamically at runtime.
For more information on these concepts see the ``Scine::Core``
`repository <https://github.com/qcscine/core>`_.

License and Copyright Information
---------------------------------

The SCINE Serenity wrapper is distributed under the BSD 3-clause "New" or
"Revised" License. For more license and copyright information, see the file
``LICENSE.txt`` in the repository.

Note: This license does not cover the original Serenity source code.
For the copyright information of the Serenity code please follow the linked
git submodule to the developers repository.

Installation and Usage
----------------------

The wrapper can be built and installed using the following commands::

    git submodule update --init
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release -DSCINE_BUILD_PYTHON_BINDINGS=ON -DCMAKE_INSTALL_PREFIX=<desired path> ..
    make -j<number of cores to use>
    make install

This will generate and install both the main Serenity code and the wrapper in
the form of the file ``serenity.module.so`` that can be used in SCINE.

In order to make Serenity available to SCINE the following two environment
variables need to be set::

    export SERENITY_RESOURCES=<desired path>/share/serenity/data/
    export SCINE_MODULE_PATH=$SCINE_MODULE_PATH:<desired path>/lib

Afterwards, SCINE programs such as `ReaDuct <https://github.com/qcscine/readuct>`_
will pick up Serenity's existence and it will be possible to request the
implemented methods.

The SCINE Serenity wrapper is also available via Python.
The underlying SCINE module can be loaded and its implemented calculators
accessed using the standard ``scine_utilities`` Python bindings.
A minimal workflow could look like this::

    import scine_utilities as su
    import scine_serenity_wrapper
    
    # Generate structure
    structure = su.AtomCollection()
    structure.elements = [su.ElementType.H, su.ElementType.H]
    structure.positions = [[-0.7, 0, 0], [0.7, 0, 0]]
    
    # Get calculator
    manager = su.core.ModuleManager.get_instance()
    calculator = manager.get('calculator', 'DFT')
    
    # Configure calculator
    calculator.structure = structure
    calculator.set_required_properties([su.Property.Energy, su.Property.Gradients])
    
    # Calculate
    results = calculator.calculate()
    print(results.energy)
    print(results.gradients)

How to Cite
-----------

When publishing results obtained with the SCINE Serenity wrapper, please cite the corresponding
release as archived on `Zenodo <https://doi.org/10.5281/zenodo.6695038>`_ (DOI
10.5281/zenodo.6695038; please use the DOI of the respective release).

This wrapper should also not be mistaken for the actual Serenity code it wraps.
For the latter code and its citations, we refer you to the original
Serenity repository. There you will find the references of the actual methods
used. They are listed in the `README.md <https://github.com/qcserenity/serenity/blob/master/README.md>`_.

Support and Contact
-------------------

In case you should encounter problems or bugs with the wrapper, please write a
short message to scine@phys.chem.ethz.ch.
