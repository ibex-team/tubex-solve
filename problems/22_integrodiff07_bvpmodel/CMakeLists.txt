# ==================================================================
#  tubex-solve - Problems
# ==================================================================

add_executable (22_integrodiff07_bvpmodel ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
target_link_libraries (22_integrodiff07_bvpmodel PUBLIC tubex-solve)