cmake_minimum_required (VERSION 3.0)

project(FreeFormHolders)

find_package(Boost COMPONENTS date_time filesystem system)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()


set(SOURCES
        ../src/main.cpp
        ../src/geodesic_distance/geodesic_matlab_api.cpp
        ../src/geodesic_distance/example0.cpp
        ../src/geodesic_distance/example1.cpp
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

target_link_libraries(FreeFormHolders ${Boost_LIBRARIES})
