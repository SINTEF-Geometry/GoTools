# vcpkg setup for Windows builds. Include this file before project() in a
# top-level CMakeLists.txt. The third party libraries are then installed by
# vcpkg from the vcpkg.json manifest next to that CMakeLists.txt.
#
# - The vcpkg toolchain is used if VCPKG_ROOT is set (or if the toolchain is
#   given with -DCMAKE_TOOLCHAIN_FILE).
# - The triplet follows the toolset: "-T v142" uses x64-windows-static-md-v142
#   (from the vcpkg-triplets folder next to this file), otherwise the native
#   compiler is used with x64-windows-static-md. Override with
#   -DVCPKG_TARGET_TRIPLET=<triplet>.
#
# On other platforms this file does nothing: there the libraries are installed
# with the system package manager (apt, Homebrew).

if(NOT CMAKE_HOST_WIN32)
  return()
endif()

if(NOT DEFINED CMAKE_TOOLCHAIN_FILE AND DEFINED ENV{VCPKG_ROOT})
  file(TO_CMAKE_PATH "$ENV{VCPKG_ROOT}" _gotools_vcpkg_root)
  set(CMAKE_TOOLCHAIN_FILE "${_gotools_vcpkg_root}/scripts/buildsystems/vcpkg.cmake"
    CACHE FILEPATH "vcpkg toolchain file")
endif()

if(NOT CMAKE_TOOLCHAIN_FILE MATCHES "vcpkg\\.cmake$")
  message(STATUS "vcpkg is not used (set VCPKG_ROOT or CMAKE_TOOLCHAIN_FILE to use it)")
  return()
endif()

list(APPEND VCPKG_OVERLAY_TRIPLETS "${CMAKE_CURRENT_LIST_DIR}/vcpkg-triplets")

if(NOT DEFINED VCPKG_TARGET_TRIPLET)
  if(CMAKE_GENERATOR_TOOLSET MATCHES "^v142")
    set(_gotools_triplet x64-windows-static-md-v142)
  else()
    set(_gotools_triplet x64-windows-static-md)
  endif()
  set(VCPKG_TARGET_TRIPLET ${_gotools_triplet} CACHE STRING "vcpkg target triplet")
endif()
message(STATUS "vcpkg triplet: ${VCPKG_TARGET_TRIPLET}")

# A BOOST_ROOT pointing to a separate Boost install would make CMake pick that
# Boost instead of the one from vcpkg (possibly built for another toolset).
if(DEFINED ENV{BOOST_ROOT})
  message(STATUS "Ignoring BOOST_ROOT=$ENV{BOOST_ROOT}, Boost is taken from vcpkg")
  unset(ENV{BOOST_ROOT})
endif()
