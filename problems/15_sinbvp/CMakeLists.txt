# ==================================================================
#  tubex-solve - Problems
# ==================================================================

add_executable (15_sinbvp ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
target_link_libraries (15_sinbvp PUBLIC tubex-solve)
