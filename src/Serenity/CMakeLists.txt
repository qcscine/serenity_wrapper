# Load lists of source files
include(Files.cmake)

# Find/Import dependecies
include(ImportUtilsOS)
import_utils_os()
include(ImportCore)
import_core()

add_library(Serenity SHARED ${SERENITY_MODULE_FILES})
set_target_properties(Serenity PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(Serenity PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)

if(WIN32)
  set_target_properties(Serenity PROPERTIES
    PREFIX ""
    IMPORT_PREFIX ""
    WINDOWS_EXPORT_ALL_SYMBOLS ON
  )
endif()

set_target_properties(Serenity PROPERTIES
  OUTPUT_NAME serenity.module
  PREFIX ""
)

target_link_libraries(Serenity
  PRIVATE
    Scine::UtilsOS
    serenity
  PUBLIC
    Scine::CoreHeaders
)
target_compile_options(Serenity PUBLIC
  $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>
)

# Add namespaced aliases
add_library(Scine::Serenity ALIAS Serenity)
add_library(Scine::SerenityModule ALIAS Serenity)

# Install
install(
  TARGETS Serenity
  EXPORT scineSerenityTargets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

# Python Bindings
if(SCINE_BUILD_PYTHON_BINDINGS)
  include(FindPythonInterpreter)
  find_python_interpreter()

  add_custom_command(TARGET Serenity POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:Serenity> ${CMAKE_CURRENT_BINARY_DIR}/scine_serenity_wrapper
    COMMENT "Copying 'serenity.module.so' to 'scine_serenity'"
  )

  install(CODE
  "execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install --prefix=${CMAKE_INSTALL_PREFIX} --upgrade --no-deps ${CMAKE_CURRENT_BINARY_DIR}
                   RESULT_VARIABLE retcode)
   if(NOT \$\{retcode\} EQUAL 0)
     message(FATAL_ERROR \"Fatal error when installing Python module using PIP.\")
   endif()"
  )
  # Copy the setup file and README.rst into the binary dir so you can install
  # using pip or make a distribution with the compiled binary
  file(
    COPY ${CMAKE_CURRENT_SOURCE_DIR}/Python/__init__.py
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/scine_serenity_wrapper
  )

  # Figure out the location of serenity's data directory dependending on
  # whether serenity is imported or in-tree
  get_target_property(_serenity_imported serenity IMPORTED)
  if(_serenity_imported)
    # From lib/cmake/serenity to share/serenity/data
    get_filename_component(_serenity_data_dir ${serenity_DIR}/../../../share/serenity/data ABSOLUTE)
  else()
    set(_serenity_data_dir ${serenity_SOURCE_DIR}/data)
  endif()
  unset(_serenity_imported)
  add_custom_command(TARGET Serenity POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${_serenity_data_dir}/initialGuess ${CMAKE_CURRENT_BINARY_DIR}/scine_serenity_wrapper/data/initialGuess
    COMMENT "Copying Serenity initial guess files to 'scine_serenity_wrapper'"
  )
  add_custom_command(TARGET Serenity POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${_serenity_data_dir}/basis ${CMAKE_CURRENT_BINARY_DIR}/scine_serenity_wrapper/data/basis
    COMMENT "Copying Serenity basis files to 'scine_serenity_wrapper'"
  )
  unset(_serenity_data_dir)

  # Figure out python dependencies
  include(TargetLibName)
  set(_module_name "serenity.module${CMAKE_SHARED_LIBRARY_SUFFIX}")
  set(serenity_PY_DEPS ", \"${_module_name}\", \"data\"")
  unset(_module_name)
  target_lib_type(Scine::UtilsOS _utils_libtype)
  if(_utils_libtype STREQUAL "SHARED_LIBRARY")
    if(APPLE)
      set_target_properties(Serenity PROPERTIES
        BUILD_WITH_INSTALL_RPATH ON
        INSTALL_RPATH "@loader_path;@loader_path/../lib"
      )
    elseif(UNIX)
      set_target_properties(Serenity PROPERTIES
        BUILD_WITH_INSTALL_RPATH ON
        INSTALL_RPATH "\$ORIGIN;\$ORIGIN/../lib"
      )
    endif()
    add_custom_command(TARGET Serenity POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:Scine::UtilsOS> ${CMAKE_CURRENT_BINARY_DIR}/scine_serenity_wrapper
      COMMENT "Copying dependent shared libraries into python package directory"
    )
    target_lib_filename(Scine::UtilsOS _utils_name)
    string(APPEND serenity_PY_DEPS ", \"${_utils_name}\"")
    unset(_utils_name)
  endif()
  unset(_utils_libtype)

  # Add setuptools file
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Python/setup.py
    ${CMAKE_CURRENT_BINARY_DIR}/setup.py
    @ONLY
  )
  unset(serenity_PY_DEPS)
  file(
    COPY ${PROJECT_SOURCE_DIR}/README.rst
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
  )

  # Copy license files
  file(
    COPY ${PROJECT_SOURCE_DIR}/LICENSE.txt ${EXTERNAL_LICENSE_FILES}
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/scine_serenity_wrapper
  )
endif()

# Component-level cmake package files
scine_install_component_cmake_files(
  COMPONENT ScineSerenity
  EXPORT_NAME scineSerenityTargets
)
