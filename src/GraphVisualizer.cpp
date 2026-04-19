#include "graph_reliability/GraphVisualizer.h"

#define _USE_MATH_DEFINES
#include <cmath>
#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace graph_reliability {

// ---------------------------------------------------------------------------
// Fruchterman-Reingold layout
// ---------------------------------------------------------------------------

std::vector<GraphVisualizer::Vec2>
GraphVisualizer::computeLayout(const Graph& g, double W, double H, int iterations)
{
    const int n = static_cast<int>(g.numVertices());
    std::vector<Vec2> pos(n);

    if (n == 0) return pos;
    if (n == 1) { pos[0] = {W / 2.0, H / 2.0}; return pos; }

    // Initial positions on a circle so the algorithm converges quickly
    const double cx = W / 2.0, cy = H / 2.0;
    const double r0 = std::min(W, H) * 0.4;
    for (int i = 0; i < n; ++i) {
        double angle = 2.0 * M_PI * i / n;
        pos[i] = {cx + r0 * std::cos(angle), cy + r0 * std::sin(angle)};
    }

    // Optimal distance between vertices
    const double k = std::sqrt(W * H / n);

    // Starting temperature (max displacement per step)
    double temp = W / 10.0;
    const double cooling = temp / iterations;  // linear cooling

    std::vector<Vec2> disp(n);

    for (int iter = 0; iter < iterations; ++iter) {
        // Reset displacements
        for (auto& d : disp) d = {0.0, 0.0};

        // --- Repulsive forces between all pairs O(n²) ---
        for (int v = 0; v < n; ++v) {
            for (int u = 0; u < n; ++u) {
                if (u == v) continue;
                double dx = pos[v].x - pos[u].x;
                double dy = pos[v].y - pos[u].y;
                double dist = std::hypot(dx, dy);
                if (dist < 0.01) dist = 0.01;
                double f = k * k / dist;        // fr(d) = k² / d
                disp[v].x += dx / dist * f;
                disp[v].y += dy / dist * f;
            }
        }

        // --- Attractive forces along edges ---
        for (int u = 0; u < n; ++u) {
            for (size_t i = g.kao_[u]; i < g.kao_[u + 1]; ++i) {
                int v = g.fo_[i];
                if (v <= u) continue;  // process each undirected edge once
                double dx = pos[v].x - pos[u].x;
                double dy = pos[v].y - pos[u].y;
                double dist = std::hypot(dx, dy);
                if (dist < 0.01) dist = 0.01;
                double f = dist * dist / k;     // fa(d) = d² / k
                double fx = dx / dist * f;
                double fy = dy / dist * f;
                disp[v].x -= fx;
                disp[v].y -= fy;
                disp[u].x += fx;
                disp[u].y += fy;
            }
        }

        // --- Apply displacement, clamp to canvas ---
        const double pad = 60.0;  // keep vertices away from the border
        for (int v = 0; v < n; ++v) {
            double d = std::hypot(disp[v].x, disp[v].y);
            if (d > 0.001) {
                double move = std::min(d, temp);
                pos[v].x += disp[v].x / d * move;
                pos[v].y += disp[v].y / d * move;
            }
            pos[v].x = std::clamp(pos[v].x, pad, W - pad);
            pos[v].y = std::clamp(pos[v].y, pad, H - pad);
        }

        temp -= cooling;
        if (temp < 0.01) temp = 0.01;
    }

    return pos;
}

// ---------------------------------------------------------------------------
// Overlap resolution (post-processing after FR)
// ---------------------------------------------------------------------------
// Iteratively pushes pairs of vertices apart until all boundaries are at
// least min_gap pixels away from each other.
// min_dist = 2*radius + min_gap (center-to-center threshold).

static void resolveOverlaps(std::vector<GraphVisualizer::Vec2>& pos,
                             double min_dist,
                             double W, double H, double pad,
                             int max_iter = 300)
{
    const int n = static_cast<int>(pos.size());
    for (int iter = 0; iter < max_iter; ++iter) {
        bool any = false;
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                double dx = pos[j].x - pos[i].x;
                double dy = pos[j].y - pos[i].y;
                double d  = std::hypot(dx, dy);
                if (d >= min_dist) continue;
                any = true;
                // Push each vertex half the required separation
                double need = (min_dist - d) * 0.5 + 0.5;  // +0.5 避免 float drift
                if (d < 0.01) { dx = 1.0; dy = 0.0; d = 1.0; }  // coincident
                double nx = dx / d, ny = dy / d;
                pos[i].x = std::clamp(pos[i].x - nx * need, pad, W - pad);
                pos[i].y = std::clamp(pos[i].y - ny * need, pad, H - pad);
                pos[j].x = std::clamp(pos[j].x + nx * need, pad, W - pad);
                pos[j].y = std::clamp(pos[j].y + ny * need, pad, H - pad);
            }
        }
        if (!any) break;
    }
}

// ---------------------------------------------------------------------------
// SVG export
// ---------------------------------------------------------------------------

static std::string fmtProb(double p) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(2) << p;
    return ss.str();
}

void GraphVisualizer::exportSVG(const Graph& g,
                                const std::string& path,
                                const Options& opts)
{
    const int n = static_cast<int>(g.numVertices());
    if (n == 0) throw std::invalid_argument("Graph has no vertices");

    const double r        = opts.vertex_radius;
    const double min_gap  = 5.0;                // minimum pixels between boundaries
    const double min_dist = 2.0 * r + min_gap;  // minimum center-to-center distance
    const double pad      = r + min_gap;

    // Auto-scale canvas so there is theoretically enough room for all vertices
    // without overlap (each needs a square of side min_dist).
    double W = opts.width, H = opts.height;
    if (n > 1) {
        double needed = static_cast<double>(n) * min_dist * min_dist * 1.8;
        if (needed > W * H) {
            double scale = std::sqrt(needed / (W * H));
            W *= scale;
            H *= scale;
        }
    }

    auto pos = computeLayout(g, W, H, opts.iterations);

    // Post-process: push overlapping vertices apart until all gaps ≥ min_gap.
    resolveOverlaps(pos, min_dist, W, H, pad + r);

    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot open output file: " + path);

    // SVG header — use actual (possibly scaled) dimensions
    f << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
      << "<svg xmlns=\"http://www.w3.org/2000/svg\""
      << " width=\"" << static_cast<int>(W) << "\" height=\"" << static_cast<int>(H) << "\""
      << " viewBox=\"0 0 " << static_cast<int>(W) << " " << static_cast<int>(H) << "\">\n"
      << "<rect width=\"100%\" height=\"100%\" fill=\"#f8f9fa\"/>\n";

    // --- Edges (trimmed to vertex boundary so lines never pass through circles) ---
    f << "<g stroke=\"#adb5bd\" stroke-width=\"1.5\">\n";
    for (int u = 0; u < n; ++u) {
        for (size_t i = g.kao_[u]; i < g.kao_[u + 1]; ++i) {
            int v = g.fo_[i];
            if (v <= u) continue;
            double dx   = pos[v].x - pos[u].x;
            double dy   = pos[v].y - pos[u].y;
            double dist = std::hypot(dx, dy);
            if (dist < 0.01) continue;  // skip coincident (shouldn't happen after resolveOverlaps)
            double nx = dx / dist, ny = dy / dist;
            // Start/end at circle boundary, not at center
            double x1 = pos[u].x + nx * r, y1 = pos[u].y + ny * r;
            double x2 = pos[v].x - nx * r, y2 = pos[v].y - ny * r;
            double mx  = (pos[u].x + pos[v].x) / 2.0;
            double my  = (pos[u].y + pos[v].y) / 2.0;
            f << "  <line x1=\"" << x1 << "\" y1=\"" << y1 << "\""
              << " x2=\"" << x2 << "\" y2=\"" << y2 << "\"/>\n";
            if (opts.show_probs) {
                double prob = g.p_array_[i];
                f << "  <text x=\"" << mx << "\" y=\"" << my
                  << "\" font-size=\"9\" fill=\"#6c757d\""
                  << " text-anchor=\"middle\" dominant-baseline=\"middle\">"
                  << fmtProb(prob) << "</text>\n";
            }
        }
    }
    f << "</g>\n";

    // --- Vertices ---
    f << "<g font-family=\"sans-serif\" font-size=\"" << (r - 2) << "\">\n";
    for (int v = 0; v < n; ++v) {
        std::string fill;
        if (v == opts.source)      fill = "#2ecc71";  // green
        else if (v == opts.target) fill = "#e74c3c";  // red
        else                       fill = "#4a90d9";  // blue

        f << "  <circle cx=\"" << pos[v].x << "\" cy=\"" << pos[v].y
          << "\" r=\"" << r << "\""
          << " fill=\"" << fill << "\""
          << " stroke=\"#2c3e50\" stroke-width=\"1.5\"/>\n"
          << "  <text x=\"" << pos[v].x << "\" y=\"" << pos[v].y
          << "\" text-anchor=\"middle\" dominant-baseline=\"middle\""
          << " fill=\"white\" font-weight=\"bold\">" << v << "</text>\n";
    }
    f << "</g>\n";

    // --- Legend ---
    f << "<g font-family=\"sans-serif\" font-size=\"11\" fill=\"#495057\">\n"
      << "  <text x=\"10\" y=\"20\">V=" << n
      << "  E=" << (g.fo_.size() / 2) << "</text>\n";
    if (opts.source >= 0)
        f << "  <circle cx=\"10\" cy=\"35\" r=\"5\" fill=\"#2ecc71\"/>"
          << "<text x=\"20\" y=\"39\">s=" << opts.source << "</text>\n";
    if (opts.target >= 0)
        f << "  <circle cx=\"10\" cy=\"52\" r=\"5\" fill=\"#e74c3c\"/>"
          << "<text x=\"20\" y=\"56\">t=" << opts.target << "</text>\n";
    f << "</g>\n";

    f << "</svg>\n";
}

// ---------------------------------------------------------------------------
// DOT export
// ---------------------------------------------------------------------------

void GraphVisualizer::exportDot(const Graph& g,
                                const std::string& path,
                                const Options& opts)
{
    const int n = static_cast<int>(g.numVertices());

    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot open output file: " + path);

    f << "graph G {\n"
      << "  overlap=false;\n"
      << "  splines=true;\n"
      << "  node [shape=circle, fixedsize=true, width=0.4, fontsize=10];\n\n";

    // Node attributes (color overrides for s/t)
    for (int v = 0; v < n; ++v) {
        if (v == opts.source)
            f << "  " << v << " [style=filled, fillcolor=\"#2ecc71\", fontcolor=white];\n";
        else if (v == opts.target)
            f << "  " << v << " [style=filled, fillcolor=\"#e74c3c\", fontcolor=white];\n";
    }
    f << "\n";

    // Edges
    for (int u = 0; u < n; ++u) {
        for (size_t i = g.kao_[u]; i < g.kao_[u + 1]; ++i) {
            int v = g.fo_[i];
            if (v <= u) continue;
            f << "  " << u << " -- " << v;
            if (opts.show_probs)
                f << " [label=\"" << fmtProb(g.p_array_[i]) << "\"]";
            f << ";\n";
        }
    }

    f << "}\n";
}

} // namespace graph_reliability
