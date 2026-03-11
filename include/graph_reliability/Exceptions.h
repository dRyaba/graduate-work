/**
 * @file Exceptions.h
 * @brief Custom exception classes for Graph Reliability Analysis Library
 * @author Graduate Work Project
 * @date 2026
 * 
 * This header defines a hierarchy of custom exceptions for better error handling
 * and more specific error reporting throughout the library.
 */

#pragma once

#include <stdexcept>
#include <string>

namespace graph_reliability {

/**
 * @brief Base exception class for all graph reliability errors
 * 
 * All custom exceptions inherit from this class, which itself inherits from
 * std::runtime_error. This allows catching all library-specific exceptions
 * with a single catch block if needed.
 */
class GraphReliabilityException : public std::runtime_error {
public:
    /**
     * @brief Constructor
     * @param message Error message
     */
    explicit GraphReliabilityException(const std::string& message)
        : std::runtime_error("GraphReliability: " + message) {}
    
    /**
     * @brief Constructor with C-string
     * @param message Error message
     */
    explicit GraphReliabilityException(const char* message)
        : std::runtime_error(std::string("GraphReliability: ") + message) {}
    
    virtual ~GraphReliabilityException() = default;
};

/**
 * @brief Exception thrown when graph structure is invalid
 * 
 * This exception is thrown when:
 * - Graph has invalid CSR format
 * - Graph has edges without vertices
 * - Graph validation fails
 */
class InvalidGraphException : public GraphReliabilityException {
public:
    explicit InvalidGraphException(const std::string& message)
        : GraphReliabilityException("Invalid graph: " + message) {}
    
    explicit InvalidGraphException(const char* message)
        : GraphReliabilityException(std::string("Invalid graph: ") + message) {}
};

/**
 * @brief Exception thrown when a file cannot be found
 * 
 * This exception is thrown when:
 * - Graph file does not exist
 * - Data directory is missing
 * - Configuration file is not found
 */
class FileNotFoundException : public GraphReliabilityException {
public:
    explicit FileNotFoundException(const std::string& filepath)
        : GraphReliabilityException("File not found: " + filepath) {}
    
    explicit FileNotFoundException(const char* filepath)
        : GraphReliabilityException(std::string("File not found: ") + filepath) {}
};

/**
 * @brief Exception thrown when file format is invalid or corrupted
 * 
 * This exception is thrown when:
 * - KAO format is malformed
 * - Edge List format is invalid
 * - File parsing fails
 * - Required data fields are missing
 */
class InvalidFormatException : public GraphReliabilityException {
public:
    explicit InvalidFormatException(const std::string& message)
        : GraphReliabilityException("Invalid format: " + message) {}
    
    explicit InvalidFormatException(const char* message)
        : GraphReliabilityException(std::string("Invalid format: ") + message) {}
};

/**
 * @brief Exception thrown when function parameters are invalid
 * 
 * This exception is thrown when:
 * - Vertex indices are out of range
 * - Diameter is negative or zero
 * - Method ID is invalid
 * - Reliability values are outside [0, 1]
 */
class InvalidParameterException : public GraphReliabilityException {
public:
    explicit InvalidParameterException(const std::string& message)
        : GraphReliabilityException("Invalid parameter: " + message) {}
    
    explicit InvalidParameterException(const char* message)
        : GraphReliabilityException(std::string("Invalid parameter: ") + message) {}
};

/**
 * @brief Exception thrown when calculation fails
 * 
 * This exception is thrown when:
 * - Reliability calculation encounters an error
 * - Decomposition fails
 * - Algorithm cannot proceed
 */
class CalculationException : public GraphReliabilityException {
public:
    explicit CalculationException(const std::string& message)
        : GraphReliabilityException("Calculation error: " + message) {}
    
    explicit CalculationException(const char* message)
        : GraphReliabilityException(std::string("Calculation error: ") + message) {}
};

} // namespace graph_reliability
