# ==================================================================
#  tubex-solve - Problems
# ==================================================================

list (APPEND SRC ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
add_executable (02_xmsin_fwd ${SRC})
target_link_libraries (02_xmsin_fwd PUBLIC tubex-solve)