# ==================================================================
#  tubex-solve - Problems
# ==================================================================

list (APPEND SRC ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
add_executable (04_bvp ${SRC})
target_link_libraries (04_bvp PUBLIC tubex-solve)