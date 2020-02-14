# ==================================================================
#  tubex-solve - Problems
# ==================================================================

add_executable (05_delay ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
target_link_libraries (05_delay PUBLIC tubex-solve)