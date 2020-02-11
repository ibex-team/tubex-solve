# ==================================================================
#  tubex-solve - Problems
# ==================================================================

list (APPEND SRC ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
add_executable (05_delay ${SRC})
target_link_libraries (05_delay PUBLIC tubex-solve)