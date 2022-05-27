cmake_minimum_required(VERSION 3.10)

include(../option.cmake)

set(TARGET ${LIBRARY_NAME})
include(srclist.cmake)
set(LIBS petsc)

add_library(${TARGET} SHARED ${SRCS})
target_link_libraries(${TARGET} PUBLIC ${LIBS})
