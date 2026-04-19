#pragma once

#include "Graph.h"
#include <string>
#include <vector>

namespace graph_reliability {

/** Options for GraphVisualizer — defined outside the class to allow default arguments. */
struct GraphVisualizerOptions {
    int  width         = 900;   ///< SVG canvas width in pixels
    int  height        = 700;   ///< SVG canvas height in pixels
    int  iterations    = 500;   ///< FR algorithm iterations
    int  vertex_radius = 14;    ///< Vertex circle radius in pixels
    int  source        = -1;    ///< Source vertex (green); -1 = none
    int  target        = -1;    ///< Target vertex (red);   -1 = none
    bool show_probs    = false; ///< Label edges with their probabilities
};

/**
 * @brief Graph visualization: Fruchterman-Reingold layout → SVG / DOT export.
 *
 * Usage:
 *   GraphVisualizerOptions opts;
 *   opts.source = 11; opts.target = 95;
 *   GraphVisualizer::exportSVG(*graph, "output.svg", opts);
 */
class GraphVisualizer {
public:
    using Options = GraphVisualizerOptions;

    /**
     * @brief Compute Fruchterman-Reingold layout and write SVG.
     * @param g     Graph in CSR format
     * @param path  Output file path (e.g. "graph.svg")
     * @param opts  Rendering options
     */
    static void exportSVG(const Graph& g,
                          const std::string& path,
                          const Options& opts = Options{});

    /**
     * @brief Write Graphviz DOT file.
     *        Run: neato -Tsvg graph.dot -o graph.svg
     * @param g     Graph in CSR format
     * @param path  Output file path (e.g. "graph.dot")
     * @param opts  Options (source/target highlighting in node attributes)
     */
    static void exportDot(const Graph& g,
                          const std::string& path,
                          const Options& opts = Options{});

private:
    struct Vec2 { double x = 0.0, y = 0.0; };

    /** Run Fruchterman-Reingold spring-embedder. Returns positions in [0,W]×[0,H]. */
    static std::vector<Vec2> computeLayout(const Graph& g,
                                           double W, double H,
                                           int iterations);
};

} // namespace graph_reliability
