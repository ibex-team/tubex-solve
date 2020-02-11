# ==================================================================
#  tubex-solve - Problems
# ==================================================================

list (APPEND SRC ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
add_executable (07_integro_diff ${SRC})
target_link_libraries (07_integro_diff PUBLIC tubex-solve)