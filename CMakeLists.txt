cmake_minimum_required(VERSION 3.16)

# Detect and prefer MinGW/GCC over MSVC if available (more portable)
if(NOT DEFINED CMAKE_C_COMPILER AND NOT DEFINED CMAKE_CXX_COMPILER)
    # Try to find MinGW/GCC first (more portable, doesn't require Visual Studio)
    find_program(CMAKE_C_COMPILER NAMES gcc mingw32-gcc clang)
    find_program(CMAKE_CXX_COMPILER NAMES g++ mingw32-g++ clang++)
endif()

project(GraphReliabilityAnalysis 
        VERSION 1.0.0
        DESCRIPTION "Graph Reliability Analysis Library"
        LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Print compiler information
message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "C++ Compiler Path: ${CMAKE_CXX_COMPILER}")

# Set build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Option for sanitizers (only for Debug builds)
option(ENABLE_SANITIZERS "Enable AddressSanitizer and UndefinedBehaviorSanitizer" OFF)

# Variable to track if logging should be enabled
set(LOGGING_ENABLED OFF)

# Build type specific options
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_definitions(GRAPH_RELIABILITY_DEBUG)
    add_compile_definitions(GRAPH_RELIABILITY_ENABLE_LOGGING)
    set(LOGGING_ENABLED ON)
    if(MSVC)
        add_compile_options(/W4 /permissive- /Od)
    else()
        add_compile_options(-Wall -Wno-invalid-utf8 -Wno-sign-compare -Wno-unused-variable -g -O0)
        # Add sanitizers if enabled
        if(ENABLE_SANITIZERS)
            add_compile_options(-fsanitize=address,undefined)
            add_link_options(-fsanitize=address,undefined)
            message(STATUS "Sanitizers enabled: AddressSanitizer, UndefinedBehaviorSanitizer")
        endif()
    endif()
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    # Release: NO logging enabled - zero overhead
    if(MSVC)
        add_compile_options(/W4 /permissive- /O2)
    else()
        add_compile_options(-Wall -Wno-invalid-utf8 -Wno-sign-compare -Wno-unused-variable -O3 -DNDEBUG)
    endif()
elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    # RelWithDebInfo: optional logging via option
    option(ENABLE_LOGGING_IN_RELWITHDEBINFO "Enable logging in RelWithDebInfo build" OFF)
    if(ENABLE_LOGGING_IN_RELWITHDEBINFO)
        add_compile_definitions(GRAPH_RELIABILITY_ENABLE_LOGGING)
        set(LOGGING_ENABLED ON)
    endif()
    if(MSVC)
        add_compile_options(/W4 /permissive- /O2)
    else()
        add_compile_options(-Wall -Wno-invalid-utf8 -Wno-sign-compare -Wno-unused-variable -O2 -g -DNDEBUG)
    endif()
elseif(CMAKE_BUILD_TYPE STREQUAL "MinSizeRel")
    # MinSizeRel: minimize size
    if(MSVC)
        add_compile_options(/W4 /permissive- /O1)
    else()
        add_compile_options(-Wall -Wno-invalid-utf8 -Wno-sign-compare -Wno-unused-variable -Os -DNDEBUG)
    endif()
endif()

# Compiler-specific options (fallback if build type not set)
if(MSVC)
    if(NOT CMAKE_BUILD_TYPE)
        add_compile_options(/W4 /permissive-)
    endif()
else()
    if(NOT CMAKE_BUILD_TYPE)
        add_compile_options(-Wall -Wno-invalid-utf8 -Wno-sign-compare -Wno-unused-variable)
    endif()
endif()

# Include directories
include_directories(include)
include_directories(.)

# Source files
set(SOURCES
    src/main.cpp
    src/CliHandlers.cpp
    src/Graph.cpp
    src/ReliabilityGraph.cpp
    src/DataImporter.cpp
    src/GraphOperations.cpp
    src/TestSuite.cpp
    src/Logger.cpp
    src/PathEnumerator.cpp
    src/CancelaPetingiState.cpp
    src/SeriesParallelTransform.cpp
    src/GraphVisualizer.cpp
)

# Header files
set(HEADERS
    include/graph_reliability.h
    include/graph_reliability/Graph.h
    include/graph_reliability/ReliabilityGraph.h
    include/graph_reliability/DataImporter.h
    include/graph_reliability/GraphOperations.h
    include/graph_reliability/TestSuite.h
    include/graph_reliability/CliHandlers.h
    include/graph_reliability/Logger.h
    include/graph_reliability/Exceptions.h
    include/graph_reliability/PathEnumerator.h
    include/graph_reliability/CancelaPetingiState.h
    include/graph_reliability/SeriesParallelTransform.h
    include/graph_reliability/GraphVisualizer.h
)

# Create executable
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

# Set target properties
set_target_properties(${PROJECT_NAME} PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
    OUTPUT_NAME "graph_reliability"
)

# Find filesystem library (needed for some compilers, e.g. older GCC)
find_library(FILESYSTEM_LIB stdc++fs)
if(FILESYSTEM_LIB)
    target_link_libraries(${PROJECT_NAME} PRIVATE ${FILESYSTEM_LIB})
endif()

# Add spdlog dependency (only if logging is enabled)
if(LOGGING_ENABLED)
    include(FetchContent)
    FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG v1.12.0
    )
    FetchContent_MakeAvailable(spdlog)
    target_link_libraries(${PROJECT_NAME} PRIVATE spdlog::spdlog)
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

# Add subdirectory for unit tests
add_subdirectory(tests)

# Add test for basic functionality
add_test(NAME BasicFunctionality
         COMMAND ${PROJECT_NAME} --help)

# Cross-platform test for data directory (single COMMAND per add_test)
if(WIN32)
    add_test(NAME DataDirectoryExists
             COMMAND powershell -Command "if (Test-Path '${CMAKE_SOURCE_DIR}/graphs_data') { exit 0 } else { exit 1 }")
else()
    add_test(NAME DataDirectoryExists
             COMMAND test -d ${CMAKE_SOURCE_DIR}/graphs_data)
endif()

# Custom target for running tests
# Ensure graphs_data is accessible from build directory
add_custom_target(run_tests
    COMMAND ${PROJECT_NAME} --test 0 test_output.csv
    DEPENDS ${PROJECT_NAME}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Running graph reliability tests"
)

# Add target for running unit tests
add_custom_target(run_unit_tests
    COMMAND graph_reliability_tests
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Running unit tests"
    DEPENDS graph_reliability_tests
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
message(STATUS "  Logging enabled: ${LOGGING_ENABLED}")
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "  Sanitizers enabled: ${ENABLE_SANITIZERS}")
endif()
message(STATUS "  Filesystem library: ${FILESYSTEM_LIB}")
message(STATUS "  Doxygen: ${DOXYGEN_FOUND}")
