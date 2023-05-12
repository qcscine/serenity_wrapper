cmake_minimum_required(VERSION 3.9)

# Setup the Scine Project
project(Serenity
  VERSION 2.0.0
  DESCRIPTION "SCINE Serenity Wrapper."
)

# Set the module path for universal cmake files inclusion
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/dev/cmake)

# Find Serenity
find_package(serenity QUIET)
if(NOT TARGET serenity)
  include(DownloadProject)
  download_project(
    PROJ serenity
    GIT_REPOSITORY https://github.com/qcserenity/serenity.git
    GIT_TAG 1.4.0
    QUIET
    UPDATE_DISCONNECTED 1
  )
  set(_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
  set(BUILD_SHARED_LIBS OFF)

  # harmonize wrapper/scine settings with serenity settings
  set(SERENITY_MARCH ${SCINE_MARCH} CACHE STRING "" FORCE)
  set(SERENITY_ENABLE_TESTS ${SCINE_BUILD_TESTS} CACHE BOOL "Activate generation of unittests" FORCE)
  set(SERENITY_USE_BLAS ${SCINE_USE_BLAS} CACHE BOOL "" FORCE)
  set(SERENITY_USE_LAPACK ${SCINE_USE_LAPACKE} CACHE BOOL "" FORCE)
  set(SERENITY_USE_INTEL_MKL ${SCINE_USE_INTEL_MKL} CACHE BOOL "" FORCE)

  add_subdirectory(${serenity_SOURCE_DIR} ${serenity_BINARY_DIR})
  if(TARGET serenity)
    message(STATUS "Serenity was not found in your PATH, so it was downloaded.")
  else()
    string(CONCAT error_msg
      "Serenity was not found in your PATH and could not be established through "
      "a download. Try specifying serenity_DIR or altering CMAKE_PREFIX_PATH to "
      "point to a candidate Serenity installation base directory."
    )
    message(FATAL_ERROR ${error_msg})
  endif()
  set(BUILD_SHARED_LIBS ${_BUILD_SHARED_LIBS})
  unset(_BUILD_SHARED_LIBS)
endif()


# Component setup
include(ComponentSetup)
scine_setup_component()

# Enable testing
if(SCINE_BUILD_TESTS)
  enable_testing()
endif()

# Subdirectories
add_subdirectory(src)

# Documentation
scine_component_documentation(UtilsOSDocumentation)
