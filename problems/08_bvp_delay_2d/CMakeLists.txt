# ==================================================================
#  tubex-solve - Problems
# ==================================================================

add_executable (08_bvp_delay_2d ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
target_link_libraries (08_bvp_delay_2d PUBLIC tubex-solve)