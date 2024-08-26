from dev.conan.base import ScineConan
import sys

class ScineSerenityConan(ScineConan):
    name = "scine_serenity_wrapper"
    version = "3.1.0"
    url = "https://github.com/qcscine/serenity_wrapper"
    description = """A wrapper around Serenity (https://github.com/qcserenity/serenity),
it exports DFT, HF and other quantum chemistry capabilities into the SCINE tool
chain."""
    options = {
        "shared": [True, False],
        "python": [True, False],
        "tests": [True, False],
        "microarch": ["detect", "none"],
        "python_version": "ANY"
    }
    python_version_string = str(sys.version_info.major) + \
        "." + str(sys.version_info.minor)
    default_options = {
        "shared": True,
        "python": False,
        "tests": False,
        "microarch": "none",
        "python_version": python_version_string,
    }
    exports = "dev/conan/*.py"
    exports_sources = [
        "dev/cmake/*", "src/*", "CMakeLists.txt", "README.rst", "LICENSE.txt",
        "dev/conan/hook.cmake", "dev/conan/glue/*"
    ]
    requires = ["scine_utilities/[=10.0.0]", "serenity/1.5.2"]
    cmake_name = "Serenity"
