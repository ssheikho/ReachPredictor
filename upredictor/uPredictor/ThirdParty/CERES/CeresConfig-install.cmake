# Ceres Solver - A fast non-linear least squares minimizer
# Copyright 2015 Google Inc. All rights reserved.
# http://ceres-solver.org/
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# * Neither the name of Google Inc. nor the names of its contributors may be
#   used to endorse or promote products derived from this software without
#   specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# Authors: pablo.speciale@gmail.com (Pablo Speciale)
#          alexs.mac@gmail.com (Alex Stewart)
#

# Config file for Ceres Solver - Find Ceres & dependencies.
#
# This file is used by CMake when find_package(Ceres) is invoked and either
# the directory containing this file either is present in CMAKE_MODULE_PATH
# (if Ceres was installed), or exists in the local CMake package registry if
# the Ceres build directory was exported.
#
# This module defines the following variables:
#
# Ceres_FOUND / CERES_FOUND: True if Ceres has been successfully
#                            found. Both variables are set as although
#                            FindPackage() only references Ceres_FOUND
#                            in Config mode, given the conventions for
#                            <package>_FOUND when FindPackage() is
#                            called in Module mode, users could
#                            reasonably expect to use CERES_FOUND
#                            instead.
#
# CERES_VERSION: Version of Ceres found.
#
# CERES_INCLUDE_DIRS: Include directories for Ceres and the
#                     dependencies which appear in the Ceres public
#                     API and are thus required to use Ceres.
#
# CERES_LIBRARIES: Libraries for Ceres and all
#                  dependencies against which Ceres was
#                  compiled. This will not include any optional
#                  dependencies that were disabled when Ceres was
#                  compiled.
#
# The following variables are also defined for legacy compatibility
# only.  Any new code should not use them as they do not conform to
# the standard CMake FindPackage naming conventions.
#
# CERES_INCLUDES = ${CERES_INCLUDE_DIRS}.

# Called if we failed to find Ceres or any of its required dependencies,
# unsets all public (designed to be used externally) variables and reports
# error message at priority depending upon [REQUIRED/QUIET/<NONE>] argument.
macro(CERES_REPORT_NOT_FOUND REASON_MSG)
  # FindPackage() only references Ceres_FOUND, and requires it to be
  # explicitly set FALSE to denote not found (not merely undefined).
  set(Ceres_FOUND FALSE)
  set(CERES_FOUND FALSE)
  unset(CERES_INCLUDE_DIRS)
  unset(CERES_LIBRARIES)

  # Reset the CMake module path to its state when this script was called.
  set(CMAKE_MODULE_PATH ${CALLERS_CMAKE_MODULE_PATH})

  # Note <package>_FIND_[REQUIRED/QUIETLY] variables defined by
  # FindPackage() use the camelcase library name, not uppercase.
  if (Ceres_FIND_QUIETLY)
    message(STATUS "Failed to find Ceres - " ${REASON_MSG} ${ARGN})
  else (Ceres_FIND_REQUIRED)
    message(FATAL_ERROR "Failed to find Ceres - " ${REASON_MSG} ${ARGN})
  else()
    # Neither QUIETLY nor REQUIRED, use SEND_ERROR which emits an error
    # that prevents generation, but continues configuration.
    message(SEND_ERROR "Failed to find Ceres - " ${REASON_MSG} ${ARGN})
  endif ()
  return()
endmacro(CERES_REPORT_NOT_FOUND)

# If Ceres was not installed, then by definition it was exported
# from a build directory.
set(CERES_WAS_INSTALLED TRUE)

# Record the state of the CMake module path when this script was
# called so that we can ensure that we leave it in the same state on
# exit as it was on entry, but modify it locally.
set(CALLERS_CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH})

# Get the (current, i.e. installed) directory containing this file.
get_filename_component(CERES_CURRENT_CONFIG_DIR
  "${CMAKE_CURRENT_LIST_FILE}" PATH)

if (CERES_WAS_INSTALLED)
  # Reset CMake module path to the installation directory of this
  # script, thus we will use the FindPackage() scripts shipped with
  # Ceres to find Ceres' dependencies, even if the user has equivalently
  # named FindPackage() scripts in their project.
  set(CMAKE_MODULE_PATH ${CERES_CURRENT_CONFIG_DIR})

  # Build the absolute root install directory as a relative path
  # (determined when Ceres was configured & built) from the current
  # install directory for this this file.  This allows for the install
  # tree to be relocated, after Ceres was built, outside of CMake.
  get_filename_component(CURRENT_ROOT_INSTALL_DIR
    ${CERES_CURRENT_CONFIG_DIR}/../
    ABSOLUTE)
  if (NOT EXISTS ${CURRENT_ROOT_INSTALL_DIR})
    ceres_report_not_found(
      "Ceres install root: ${CURRENT_ROOT_INSTALL_DIR}, "
      "determined from relative path from CeresConfig.cmake install location: "
      "${CERES_CURRENT_CONFIG_DIR}, does not exist. Either the install "
      "directory was deleted, or the install tree was only partially relocated "
      "outside of CMake after Ceres was built.")
  endif (NOT EXISTS ${CURRENT_ROOT_INSTALL_DIR})

  # Set the include directories for Ceres (itself).
  set(CERES_INCLUDE_DIR "${CURRENT_ROOT_INSTALL_DIR}/include")
  if (NOT EXISTS ${CERES_INCLUDE_DIR}/ceres/ceres.h)
    ceres_report_not_found(
      "Ceres install root: ${CURRENT_ROOT_INSTALL_DIR}, "
      "determined from relative path from CeresConfig.cmake install location: "
      "${CERES_CURRENT_CONFIG_DIR}, does not contain Ceres headers. "
      "Either the install directory was deleted, or the install tree was only "
      "partially relocated outside of CMake after Ceres was built.")
  endif (NOT EXISTS ${CERES_INCLUDE_DIR}/ceres/ceres.h)
  list(APPEND CERES_INCLUDE_DIRS ${CERES_INCLUDE_DIR})

else(CERES_WAS_INSTALLED)
  # Ceres was exported from the build tree.
  set(CERES_EXPORTED_BUILD_DIR ${CERES_CURRENT_CONFIG_DIR})
  get_filename_component(CERES_EXPORTED_SOURCE_DIR
    ${CERES_EXPORTED_BUILD_DIR}/../
    ABSOLUTE)
  if (NOT EXISTS ${CERES_EXPORTED_SOURCE_DIR})
    ceres_report_not_found(
      "Ceres exported source directory: ${CERES_EXPORTED_SOURCE_DIR}, "
      "determined from relative path from CeresConfig.cmake exported build "
      "directory: ${CERES_EXPORTED_BUILD_DIR} does not exist.")
  endif()

  # Reset CMake module path to the cmake directory in the Ceres source
  # tree which was exported, thus we will use the FindPackage() scripts shipped
  # with Ceres to find Ceres' dependencies, even if the user has equivalently
  # named FindPackage() scripts in their project.
  set(CMAKE_MODULE_PATH ${CERES_EXPORTED_SOURCE_DIR}/cmake)

  # Set the include directories for Ceres (itself).
  set(CERES_INCLUDE_DIR "${CERES_EXPORTED_SOURCE_DIR}/include")
  if (NOT EXISTS ${CERES_INCLUDE_DIR}/ceres/ceres.h)
    ceres_report_not_found(
      "Ceres exported source directory: ${CERES_EXPORTED_SOURCE_DIR}, "
      "determined from relative path from CeresConfig.cmake exported build "
      "directory: ${CERES_EXPORTED_BUILD_DIR}, does not contain Ceres headers.")
  endif (NOT EXISTS ${CERES_INCLUDE_DIR}/ceres/ceres.h)
  list(APPEND CERES_INCLUDE_DIRS ${CERES_INCLUDE_DIR})

  # Append the path to the configured config.h in the exported build directory
  # to the Ceres include directories.
  set(CERES_CONFIG_FILE
    ${CERES_EXPORTED_BUILD_DIR}/config/ceres/internal/config.h)
  if (NOT EXISTS ${CERES_CONFIG_FILE})
    ceres_report_not_found(
      "Ceres exported build directory: ${CERES_EXPORTED_BUILD_DIR}, "
      "does not contain required configured Ceres config.h, it is not here: "
      "${CERES_CONFIG_FILE}.")
  endif (NOT EXISTS ${CERES_CONFIG_FILE})
  list(APPEND CERES_INCLUDE_DIRS ${CERES_EXPORTED_BUILD_DIR}/config)
endif(CERES_WAS_INSTALLED)

# Set the version.
set(CERES_VERSION 1.11.0 )

# Eigen.
# Flag set during configuration and build of Ceres.
set(CERES_EIGEN_VERSION 3.2.4)
# Append the locations of Eigen when Ceres was built to the search path hints.
list(APPEND EIGEN_INCLUDE_DIR_HINTS C:/Program Files/Eigen/include/eigen3)
# Search quietly to control the timing of the error message if not found. The
# search should be for an exact match, but for usability reasons do a soft
# match and reject with an explanation below.
find_package(Eigen ${CERES_EIGEN_VERSION} QUIET)
if (EIGEN_FOUND)
  if (NOT EIGEN_VERSION VERSION_EQUAL CERES_EIGEN_VERSION)
    # CMake's VERSION check in FIND_PACKAGE() will accept any version >= the
    # specified version. However, only version = is supported. Improve
    # usability by explaining why we don't accept non-exact version matching.
    ceres_report_not_found("Found Eigen dependency, but the version of Eigen "
      "found (${EIGEN_VERSION}) does not exactly match the version of Eigen "
      "Ceres was compiled with (${CERES_EIGEN_VERSION}). This can cause subtle "
      "bugs by triggering violations of the One Definition Rule. See the "
      "Wikipedia article http://en.wikipedia.org/wiki/One_Definition_Rule "
      "for more details")
  endif ()
  message(STATUS "Found required Ceres dependency: "
    "Eigen version ${CERES_EIGEN_VERSION} in ${EIGEN_INCLUDE_DIRS}")
else (EIGEN_FOUND)
  ceres_report_not_found("Missing required Ceres "
    "dependency: Eigen version ${CERES_EIGEN_VERSION}, please set "
    "EIGEN_INCLUDE_DIR.")
endif (EIGEN_FOUND)
list(APPEND CERES_INCLUDE_DIRS ${EIGEN_INCLUDE_DIRS})

# Glog.
# Flag set during configuration and build of Ceres.
set(CERES_USES_MINIGLOG OFF)
if (CERES_USES_MINIGLOG)
  set(MINIGLOG_INCLUDE_DIR ${CERES_INCLUDE_DIR}/ceres/internal/miniglog)
  if (NOT EXISTS ${MINIGLOG_INCLUDE_DIR})
    ceres_report_not_found(
      "Ceres include directory: "
      "${CERES_INCLUDE_DIR} does not include miniglog, but Ceres was "
      "compiled with MINIGLOG enabled (in place of Glog).")
  endif (NOT EXISTS ${MINIGLOG_INCLUDE_DIR})
  list(APPEND CERES_INCLUDE_DIRS ${MINIGLOG_INCLUDE_DIR})
  # Output message at standard log level (not the lower STATUS) so that
  # the message is output in GUI during configuration to warn user.
  message("-- Found Ceres compiled with miniglog substitute "
    "for glog, beware this will likely cause problems if glog is later linked.")
else (CERES_USES_MINIGLOG)
  # Append the locations of glog when Ceres was built to the search path hints.
  list(APPEND GLOG_INCLUDE_DIR_HINTS D:\upredictor\uPredictor\ThirdParty\glog\src\windows)
  get_filename_component(CERES_BUILD_GLOG_LIBRARY_DIR D:\upredictor\uPredictor\ThirdParty\glog\x64\Release PATH)
  list(APPEND GLOG_LIBRARY_DIR_HINTS ${CERES_BUILD_GLOG_LIBRARY_DIR})

  # Search quietly s/t we control the timing of the error message if not found.
  find_package(Glog QUIET)
  if (GLOG_FOUND)
    message(STATUS "Found required Ceres dependency: "
      "Glog in ${GLOG_INCLUDE_DIRS}")
  else (GLOG_FOUND)
    ceres_report_not_found("Missing required Ceres "
      "dependency: Glog, please set GLOG_INCLUDE_DIR.")
  endif (GLOG_FOUND)
  list(APPEND CERES_INCLUDE_DIRS ${GLOG_INCLUDE_DIRS})
endif (CERES_USES_MINIGLOG)

# Import exported Ceres targets, if they have not already been imported.
if (NOT TARGET ceres AND NOT Ceres_BINARY_DIR)
  include(${CERES_CURRENT_CONFIG_DIR}/CeresTargets.cmake)
endif (NOT TARGET ceres AND NOT Ceres_BINARY_DIR)
# Set the expected XX_LIBRARIES variable for FindPackage().
set(CERES_LIBRARIES ceres)

# Make user aware of any compile flags that will be added to their targets
# which use Ceres (i.e. flags exported in the Ceres target).  Only CMake
# versions >= 2.8.12 support target_compile_options().
if (TARGET ${CERES_LIBRARIES} AND
    NOT CMAKE_VERSION VERSION_LESS "2.8.12")
  get_target_property(CERES_INTERFACE_COMPILE_OPTIONS
    ${CERES_LIBRARIES} INTERFACE_COMPILE_OPTIONS)

  if (CERES_WAS_INSTALLED)
    set(CERES_LOCATION "${CURRENT_ROOT_INSTALL_DIR}")
  else()
    set(CERES_LOCATION "${CERES_EXPORTED_BUILD_DIR}")
  endif()

  # Check for -std=c++11 flags.
  if (CERES_INTERFACE_COMPILE_OPTIONS MATCHES ".*std=c\\+\\+11.*")
    message(STATUS "Ceres version ${CERES_VERSION} detected here: "
      "${CERES_LOCATION} was built with C++11. Ceres target will add "
      "C++11 flags to compile options for targets using it.")
  endif()
endif()

# Set legacy include directories variable for backwards compatibility.
set(CERES_INCLUDES ${CERES_INCLUDE_DIRS})

# Reset CMake module path to its state when this script was called.
set(CMAKE_MODULE_PATH ${CALLERS_CMAKE_MODULE_PATH})

# As we use CERES_REPORT_NOT_FOUND() to abort, if we reach this point we have
# found Ceres and all required dependencies.
if (CERES_WAS_INSTALLED)
  message(STATUS "Found Ceres version: ${CERES_VERSION} "
    "installed in: ${CURRENT_ROOT_INSTALL_DIR}")
else (CERES_WAS_INSTALLED)
  message(STATUS "Found Ceres version: ${CERES_VERSION} "
    "exported from build directory: ${CERES_EXPORTED_BUILD_DIR}")
endif()

# Set CERES_FOUND to be equivalent to Ceres_FOUND, which is set to
# TRUE by FindPackage() if this file is found and run, and after which
# Ceres_FOUND is not (explicitly, i.e. undefined does not count) set
# to FALSE.
set(CERES_FOUND TRUE)
