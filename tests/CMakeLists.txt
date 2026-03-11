# Google Test setup
include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.14.0
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# Enable testing
enable_testing()

# Test sources
set(TEST_SOURCES
    test_main.cpp
    test_Graph.cpp
    test_ReliabilityGraph.cpp
    test_DataImporter.cpp
    test_GraphOperations.cpp
)

# Create test executable first, then set include dirs
add_executable(graph_reliability_tests ${TEST_SOURCES})
target_include_directories(graph_reliability_tests PRIVATE ${CMAKE_SOURCE_DIR}/include)

# Link libraries
target_link_libraries(graph_reliability_tests
    PRIVATE
    gtest_main
    gtest
)

# When parent has logging enabled, tests need spdlog too (Logger.h includes it)
if(TARGET spdlog::spdlog)
    target_compile_definitions(graph_reliability_tests PRIVATE GRAPH_RELIABILITY_ENABLE_LOGGING)
    target_link_libraries(graph_reliability_tests PRIVATE spdlog::spdlog)
endif()

# Link against the main library sources (since it's not a separate library)
target_sources(graph_reliability_tests PRIVATE
    ${CMAKE_SOURCE_DIR}/src/Graph.cpp
    ${CMAKE_SOURCE_DIR}/src/ReliabilityGraph.cpp
    ${CMAKE_SOURCE_DIR}/src/DataImporter.cpp
    ${CMAKE_SOURCE_DIR}/src/GraphOperations.cpp
    ${CMAKE_SOURCE_DIR}/src/Logger.cpp
    ${CMAKE_SOURCE_DIR}/src/PathEnumerator.cpp
    ${CMAKE_SOURCE_DIR}/src/CancelaPetingiState.cpp
    ${CMAKE_SOURCE_DIR}/src/SeriesParallelTransform.cpp
)

# Set C++ standard
set_target_properties(graph_reliability_tests PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
)

# Add tests
add_test(NAME GraphTests COMMAND graph_reliability_tests --gtest_filter=GraphTest.*)
add_test(NAME ReliabilityGraphTests COMMAND graph_reliability_tests --gtest_filter=ReliabilityGraphTest.*)
add_test(NAME DataImporterTests COMMAND graph_reliability_tests --gtest_filter=DataImporterTest.*)
add_test(NAME GraphOperationsTests COMMAND graph_reliability_tests --gtest_filter=GraphOperationsTest.*)

# Set working directory for tests (so they can find graphs_data/)
set_tests_properties(
    GraphTests
    ReliabilityGraphTests
    DataImporterTests
    GraphOperationsTests
    PROPERTIES
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
