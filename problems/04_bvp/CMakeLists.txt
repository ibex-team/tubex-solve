# ==================================================================
#  tubex-solve - Problems
# ==================================================================

add_executable (04_bvp ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
target_link_libraries (04_bvp PUBLIC tubex-solve)