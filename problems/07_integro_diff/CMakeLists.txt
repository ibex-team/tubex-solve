# ==================================================================
#  tubex-solve - Problems
# ==================================================================

add_executable (07_integro_diff ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
target_link_libraries (07_integro_diff PUBLIC tubex-solve)