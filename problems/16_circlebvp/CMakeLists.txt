# ==================================================================
#  tubex-solve - Problems
# ==================================================================

add_executable (16_circlebvp ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
target_link_libraries (16_circlebvp PUBLIC tubex-solve)
