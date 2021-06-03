# ==================================================================
#  tubex-solve - cmake configuration file
# ==================================================================

cmake_minimum_required (VERSION 3.5.1)

project (TUBEX_SOLVE VERSION 0.42 LANGUAGES CXX)
set (TUBEX_SOLVE_DESCRIPTION "Short description of the plugin")
set (IBEX_MIN_VERSION_REQ 2.8.9.1)

# If IBEX_DIR was given on the command-line and not in the environment, set it
# as an environment variable. Enable the usage of 'cmake -DIBEX_DIR=<path> ...'
if (NOT DEFINED ENV{IBEX_DIR} AND DEFINED IBEX_DIR)
  set (ENV{IBEX_DIR} "${IBEX_DIR}")
endif ()
find_package (IBEX ${IBEX_MIN_VERSION_REQ} REQUIRED)
message (STATUS "Found Ibex version ${IBEX_VERSION}")
ibex_init_common () # Ibex should have installed this function

# FIXME: the TUBEX package should have IBEX as a dependency and it should be
# sufficient to look for TUBEX.
find_package (TUBEX REQUIRED tubex)
message (STATUS "Found Tubex version ${TUBEX_VERSION}")

################################################################################
# Setting some paths
################################################################################
get_filename_component(TUBELIB_DATA_DIR    "./make"     ABSOLUTE)

set (CMAKE_BINDIR "./make/bin")

################################################################################
# Compile sources
################################################################################
add_subdirectory (src)
add_subdirectory (problems)

################################################################################
# Tests
################################################################################
option(BUILD_TESTS "Build test" OFF)
if(BUILD_TESTS)
  enable_testing(true)
  include(CTest)

  add_test(NAME solver_01
           COMMAND ./problems/01_picard/01_picard 0)
  add_test(NAME solver_02
           COMMAND ./problems/02_xmsin_fwd/02_xmsin_fwd 0)
  #add_test(NAME solver_03
  #         COMMAND ./problems/03_xmsin_bwd/03_xmsin_bwd 0)
  add_test(NAME solver_04
           COMMAND ./problems/04_bvp/04_bvp 0)
  # add_test(NAME solver_05
  #          COMMAND ./problems/05_delay/05_delay 0)
  # add_test(NAME solver_06
  #          COMMAND ./problems/06_bvp_delay/06_bvp_delay 0)
  add_test(NAME solver_07
           COMMAND ./problems/07_integro_diff/07_integro_diff 0)
  # add_test(NAME solver_08
  #          COMMAND ./problems/08_bvp_delay_2d/08_bvp_delay_2d 0)
  add_test(NAME solver_09
           COMMAND ./problems/09_csdp/09_csdp 0)
  add_test(NAME solver_10
           COMMAND ./problems/10_large_initvalue/10_large_initvalue 0)
endif()
