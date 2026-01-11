cmake_minimum_required(VERSION 3.16)
project(MultiExeProject LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Collect all cpp files recursively from src/
file(GLOB_RECURSE ALL_SOURCES CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
)

# For each .cpp, create a standalone executable
foreach(SRC_FILE ${ALL_SOURCES})

    # Get filename without extension
    get_filename_component(EXE_NAME ${SRC_FILE} NAME_WE)

    # Make sure subfolder names don’t collide
    file(RELATIVE_PATH REL_PATH ${CMAKE_CURRENT_SOURCE_DIR}/src ${SRC_FILE})
    string(REPLACE "/" "_" REL_PATH ${REL_PATH})
    string(REPLACE ".cpp" "" REL_PATH ${REL_PATH})
    set(EXE_NAME ${REL_PATH})

    add_executable(${EXE_NAME} ${SRC_FILE})

    target_include_directories(${EXE_NAME}
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
    )

endforeach()
