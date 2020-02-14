# ==================================================================
#  tubex-solve - Problems
# ==================================================================

add_executable (10_large_initvalue ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
target_link_libraries (10_large_initvalue PUBLIC tubex-solve)