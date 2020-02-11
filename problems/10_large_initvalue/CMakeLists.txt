# ==================================================================
#  tubex-solve - Problems
# ==================================================================

list (APPEND SRC ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
add_executable (10_large_initvalue ${SRC})
target_link_libraries (10_large_initvalue PUBLIC tubex-solve)