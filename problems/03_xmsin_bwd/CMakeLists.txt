# ==================================================================
#  tubex-solve - Problems
# ==================================================================

add_executable (03_xmsin_bwd ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
target_link_libraries (03_xmsin_bwd PUBLIC tubex-solve)