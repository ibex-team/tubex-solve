# ==================================================================
#  tubex-solve - Problems
# ==================================================================

list (APPEND SRC ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
add_executable (03_xmsin_bwd ${SRC})
target_link_libraries (03_xmsin_bwd PUBLIC tubex-solve)