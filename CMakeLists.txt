cmake_minimum_required(VERSION 3.16)
project(fluxval VERSION 1.0 LANGUAGES CXX)

# ---------------------------------------------------------------------------
# C++ standard.  ROOT and this project must agree, or you get ABI errors that
# look like missing symbols.  ROOT reports the standard it was built with in
# ROOT_CXX_STANDARD (6.28+); fall back to 17 otherwise, and let the user
# override with -DCMAKE_CXX_STANDARD=20.
# ---------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "" FORCE)
endif()
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(ROOT 6.22 REQUIRED
  COMPONENTS Core RIO Net Hist Graf Graf3d Gpad Tree MathCore ROOTDataFrame)

if(NOT CMAKE_CXX_STANDARD)
  if(DEFINED ROOT_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD ${ROOT_CXX_STANDARD})
  else()
    set(CMAKE_CXX_STANDARD 17)
  endif()
endif()
message(STATUS "fluxval: building against ROOT ${ROOT_VERSION}, C++${CMAKE_CXX_STANDARD}")

find_package(Threads REQUIRED)

# ---------------------------------------------------------------------------
# The statistical core is header-only and deliberately ROOT-free, so that it
# can be unit-tested without ROOT (see the selftest target below).
# ---------------------------------------------------------------------------
add_library(fluxval_core INTERFACE)
target_include_directories(fluxval_core INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(fluxval_core INTERFACE Threads::Threads)

# ---------------------------------------------------------------------------
# The analysis executable
# ---------------------------------------------------------------------------
add_executable(flux_compare
  src/main.cxx
  src/Load.cxx
  src/Plot.cxx
  src/Log.cxx)
target_link_libraries(flux_compare PRIVATE
  fluxval_core
  ROOT::Core ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Gpad ROOT::Tree
  ROOT::MathCore ROOT::ROOTDataFrame)

add_executable(flux_weight_audit
  src/weight_audit.cxx
  src/Load.cxx
  src/Plot.cxx
  src/Log.cxx)
target_link_libraries(flux_weight_audit PRIVATE
  fluxval_core
  ROOT::Core ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Gpad ROOT::Tree
  ROOT::MathCore ROOT::ROOTDataFrame)
install(TARGETS flux_weight_audit RUNTIME DESTINATION bin)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(flux_compare PRIVATE -Wall -Wextra -Wno-unused-parameter)
endif()

# ---------------------------------------------------------------------------
# Tests -- no ROOT needed, so they run anywhere
# ---------------------------------------------------------------------------
enable_testing()
add_executable(flux_selftest test/core_selftest.cxx)
target_link_libraries(flux_selftest PRIVATE fluxval_core)
add_test(NAME core_selftest COMMAND flux_selftest)

install(TARGETS flux_compare RUNTIME DESTINATION bin)
