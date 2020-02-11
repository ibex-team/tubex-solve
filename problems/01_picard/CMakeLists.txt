# ==================================================================
#  tubex-solve - Problems
# ==================================================================

list (APPEND SRC ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
add_executable (01_picard ${SRC})
target_link_libraries (01_picard PUBLIC tubex-solve)