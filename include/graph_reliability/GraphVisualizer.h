#pragma once

#include "ReliabilityGraph.h"
#include <string>
#include <vector>

namespace graph_reliability {

/** Options for GraphVisualizer — defined outside the class to allow default arguments. */
struct GraphVisualizerOptions {
    int  width         = 900;   ///< SVG canvas width in pixels (auto-scaled for large graphs)
    int  height        = 700;   ///< SVG canvas height in pixels
    int  iterations    = 500;   ///< FR algorithm iterations
    int  vertex_radius = 14;    ///< Vertex circle radius in pixels
    int  source        = -1;    ///< Source vertex (green); -1 = none
    int  target        = -1;    ///< Target vertex (red);   -1 = none
    bool show_probs    = false; ///< Label edges with their probabilities
};

/**
 * @brief Graph visualization: two-level block-aware FR layout → SVG / DOT.
 *
 * Layout algorithm:
 *   1. Decompose graph into biconnected blocks.
 *   2. Run FR on the block graph to place block centres.
 *   3. Arrange vertices within each block on a local circle.
 *   4. Short FR refinement pass on the full graph.
 *   5. Iterative overlap resolution (min gap between vertex boundaries).
 *
 * Vertex colours: blue = regular, green = source, red = target, orange = AP.
 * Block backgrounds drawn as coloured circles.
 */
class GraphVisualizer {
public:
    using Options = GraphVisualizerOptions;

    /** Two-level block FR layout → SVG file. */
    static void exportSVG(const ReliabilityGraph& g,
                          const std::string& path,
                          const Options& opts = Options{});

    /** Graphviz DOT file with cluster subgraphs per block.
     *  Run: neato -Tsvg graph.dot -o graph.svg */
    static void exportDot(const ReliabilityGraph& g,
                          const std::string& path,
                          const Options& opts = Options{});

    struct Vec2 { double x = 0.0, y = 0.0; };

private:
    /** Two-level block-aware FR layout. Returns positions in canvas space. */
    static std::vector<Vec2> computeLayout(const ReliabilityGraph& g,
                                           double W, double H,
                                           int iterations);
};

} // namespace graph_reliability
