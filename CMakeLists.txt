cmake_minimum_required (VERSION 3.0)

project(FreeFormHolders)

set (CMAKE_CXX_STANDARD 14)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL "Use OpenGL" ON)
option(LIBIGL_WITH_OPENGL_GLFW "Use GLFW" ON)

find_package(LIBIGL REQUIRED QUIET)
include_directories(include ../libigl/include)

set(SOURCES
        src/main.cpp
        )
set(HEADERS
        include/main.h
        )

add_executable(FreeFormHolders ${SOURCES} ${HEADERS})

target_link_libraries(FreeFormHolders
        igl::core
        igl::opengl
        igl::opengl_glfw
        igl::opengl_glfw_imgui
        )
