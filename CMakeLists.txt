cmake_minimum_required (VERSION 3.0)

project(FreeFormHolders)

set (CMAKE_CXX_STANDARD 14)

#set(Boost_DEBUG ON)
set(BOOST_ROOT "/usr/include/boost/")
set(BOOST_INCLUDEDIR "/usr/include")
set(BOOST_LIBRARYDIR "/usr/lib/x86_64-linux-gnu/")

find_package(Boost REQUIRED COMPONENTS
        date_time
        filesystem
        )


if(Boost_FOUND)
    include_directories(include ${Boost_INCLUDE_DIRS})
endif()


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL "Use OpenGL" ON)
option(LIBIGL_WITH_OPENGL_GLFW "Use GLFW" ON)

find_package(LIBIGL REQUIRED QUIET)
include_directories(include ../libigl/include)

set(SOURCES
        ../src/main.cpp
#        ../src/geodesic_distance/geodesic_matlab_api.cpp
#        ../src/geodesic_distance/example0.cpp
#        ../src/geodesic_distance/example1.cpp
        )
set(HEADERS
        ../include/main.h
        ../include/geodesic_distance/geodesic_algorithm_base.h
        ../include/geodesic_distance/geodesic_algorithm_dijkstra.h
        ../include/geodesic_distance/geodesic_algorithm_dijkstra_alternative.h
        ../include/geodesic_distance/geodesic_algorithm_exact.h
        ../include/geodesic_distance/geodesic_algorithm_exact_elements.h
        ../include/geodesic_distance/geodesic_algorithm_graph_base.h
        ../include/geodesic_distance/geodesic_algorithm_subdivision.h
        ../include/geodesic_distance/geodesic_constants_and_simple_functions.h
        ../include/geodesic_distance/geodesic_matlab_api.h
        ../include/geodesic_distance/geodesic_memory.h
        ../include/geodesic_distance/geodesic_mesh.h
        ../include/geodesic_distance/geodesic_mesh_elements.h
        )

add_executable(FreeFormHolders ${SOURCES} ${HEADERS})

target_link_libraries(FreeFormHolders ${Boost_LIBRARIES} igl::core igl::opengl_glfw)
