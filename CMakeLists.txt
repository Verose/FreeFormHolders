cmake_minimum_required (VERSION 3.0)

project(FreeFormHolders)

set(SOURCES
        ../src/main.cpp
        )
set(HEADERS
        ../include/main.h
        )

add_executable(FreeFormHolders ${SOURCES} ${HEADERS})

