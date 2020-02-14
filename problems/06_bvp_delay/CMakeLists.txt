# ==================================================================
#  tubex-solve - Problems
# ==================================================================

add_executable (06_bvp_delay ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
target_link_libraries (06_bvp_delay PUBLIC tubex-solve)