__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import os
import setuptools

def package_files(directory):
    paths = []
    os.chdir('scine_serenity_wrapper')
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join(path, filename))
    os.chdir('..')
    return paths

# Read README.rst for the long description
with open("README.rst", "r", encoding="utf-8") as fh:
    long_description = fh.read()


class EmptyListWithLength(list):
    """ Makes the wheel a binary distribution and platlib compliant. """

    def __len__(self):
        return 1


# Define the setup
setuptools.setup(
    name="scine_serenity_wrapper",
    version="@Serenity_VERSION@",
    author="ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group",
    author_email="scine@phys.chem.ethz.ch",
    description="A wrapper for Serenity",
    long_description=long_description,
    url="https://www.scine.ethz.ch",
    packages=["scine_serenity_wrapper"],
    include_package_data=True,
    package_data={"scine_serenity_wrapper": ['*.txt' @serenity_PY_DEPS@, "Tests/*"]},
    install_requires=["scine_utilities"],
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: C++",
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Chemistry"
    ],
    zip_safe=False,
    test_suite='pytest',
    tests_require=['pytest'],
    ext_modules=EmptyListWithLength()
)
