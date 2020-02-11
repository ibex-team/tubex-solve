# ==================================================================
#  tubex-solve - Problems
# ==================================================================

list (APPEND SRC ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
add_executable (09_csdp ${SRC})
target_link_libraries (09_csdp PUBLIC tubex-solve)