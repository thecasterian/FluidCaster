cmake_minimum_required(VERSION 3.10)

include(../option.cmake)

add_executable(${LIBRARY_NAME}-bin main.c)
target_link_libraries(${LIBRARY_NAME}-bin ${LIBRARY_NAME})
set_target_properties(${LIBRARY_NAME}-bin
    PROPERTIES
    OUTPUT_NAME ${LIBRARY_NAME}
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}
)
