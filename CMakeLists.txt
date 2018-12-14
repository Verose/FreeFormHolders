cmake_minimum_required (VERSION 3.0)

project(FreeFormHolders)

set(SOURCES
        ../src/bla.cpp
        )
set(HEADERS
        ../include/bla.h
        )

add_executable(FreeFormHolders ${SOURCES} ${HEADERS})

