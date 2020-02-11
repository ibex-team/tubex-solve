# ==================================================================
#  tubex-solve - Problems
# ==================================================================

list (APPEND SRC ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
add_executable (08_bvp_delay_2d ${SRC})
target_link_libraries (08_bvp_delay_2d PUBLIC tubex-solve)