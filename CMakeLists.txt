cmake_minimum_required(VERSION 3.16)

project(GraphReliabilityAnalysis 
        VERSION 1.0.0
        DESCRIPTION "Graph Reliability Analysis Library"
        LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler-specific options
if(MSVC)
    add_compile_options(/W4 /permissive-)
else()
    add_compile_options(-Wall -Wextra -Wpedantic -Werror)
endif()

# Include directories
include_directories(include)

# Source files
set(SOURCES
    src/main.cpp
    src/Graph.cpp
    src/ReliabilityGraph.cpp
    src/DataImporter.cpp
    src/TestSuite.cpp
    src/GraphOperations.cpp
)

# Header files
set(HEADERS
    include/graph_reliability.h
    include/graph_reliability/Graph.h
    include/graph_reliability/ReliabilityGraph.h
    include/graph_reliability/DataImporter.h
    include/graph_reliability/TestSuite.h
    include/graph_reliability/GraphOperations.h
)

# Create executable
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

# Set target properties
set_target_properties(${PROJECT_NAME} PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
    OUTPUT_NAME "graph_reliability"
)

# Find OpenMP if available
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)
    target_compile_definitions(${PROJECT_NAME} PRIVATE HAVE_OPENMP)
endif()

# Find filesystem library (needed for some compilers)
find_library(FILESYSTEM_LIB stdc++fs)
if(FILESYSTEM_LIB)
    target_link_libraries(${PROJECT_NAME} ${FILESYSTEM_LIB})
endif()

# Installation
install(TARGETS ${PROJECT_NAME}
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(FILES ${HEADERS}
    DESTINATION include/graph_reliability
)

# Create data directory
install(DIRECTORY graphs_data/
    DESTINATION share/${PROJECT_NAME}/data
    FILES_MATCHING PATTERN "*.txt"
)

# CPack configuration
set(CPACK_PACKAGE_NAME "GraphReliabilityAnalysis")
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY ${PROJECT_DESCRIPTION})
set(CPACK_PACKAGE_VENDOR "Graduate Work Project")
set(CPACK_PACKAGE_CONTACT "graduate.work@example.com")

include(CPack)

# Testing
enable_testing()

# Add test for basic functionality
add_test(NAME BasicFunctionality
         COMMAND ${PROJECT_NAME} --help)

# Add test for data directory
add_test(NAME DataDirectoryExists
         COMMAND test -d graphs_data)

# Custom target for running tests
add_custom_target(run_tests
    COMMAND ${PROJECT_NAME} --test 0 test_output.csv
    DEPENDS ${PROJECT_NAME}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Running graph reliability tests"
)

# Documentation target
find_package(Doxygen)
if(DOXYGEN_FOUND)
    configure_file(${CMAKE_SOURCE_DIR}/docs/Doxyfile.in
                   ${CMAKE_BINARY_DIR}/Doxyfile @ONLY)
    
    add_custom_target(docs
        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM
    )
endif()

# Print configuration summary
message(STATUS "Graph Reliability Analysis Configuration:")
message(STATUS "  Version: ${PROJECT_VERSION}")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  OpenMP: ${OpenMP_CXX_FOUND}")
message(STATUS "  Filesystem library: ${FILESYSTEM_LIB}")
message(STATUS "  Doxygen: ${DOXYGEN_FOUND}")
