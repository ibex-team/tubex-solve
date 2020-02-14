# ==================================================================
#  tubex-solve - Problems
# ==================================================================

add_executable (09_csdp ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
target_link_libraries (09_csdp PUBLIC tubex-solve)