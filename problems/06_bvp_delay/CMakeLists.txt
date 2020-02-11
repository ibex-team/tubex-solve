# ==================================================================
#  tubex-solve - Problems
# ==================================================================

list (APPEND SRC ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
add_executable (06_bvp_delay ${SRC})
target_link_libraries (06_bvp_delay PUBLIC tubex-solve)