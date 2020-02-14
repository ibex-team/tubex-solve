# ==================================================================
#  tubex-solve - Problems
# ==================================================================

add_executable (01_picard ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
target_link_libraries (01_picard PUBLIC tubex-solve)