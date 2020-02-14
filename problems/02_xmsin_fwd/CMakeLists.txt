# ==================================================================
#  tubex-solve - Problems
# ==================================================================

add_executable (02_xmsin_fwd ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
target_link_libraries (02_xmsin_fwd PUBLIC tubex-solve)