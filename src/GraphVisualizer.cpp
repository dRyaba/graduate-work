#define _USE_MATH_DEFINES
#include <cmath>
#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

#include "graph_reliability/GraphVisualizer.h"
#include "graph_reliability/ReliabilityGraph.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <set>
#include <sstream>
#include <stdexcept>

namespace graph_reliability {

using V2 = GraphVisualizer::Vec2;

// ---------------------------------------------------------------------------
// Generic Fruchterman-Reingold
// No border clamping — caller normalises afterwards.
// start_temp: initial max displacement (default W/10).
// ---------------------------------------------------------------------------
static std::vector<V2> runFR(
    int n,
    const std::vector<std::pair<int,int>>& edges,
    std::vector<V2> pos,
    double W, double H,
    int iterations,
    double start_temp = -1.0)
{
    if (n <= 1) return pos;

    double k    = std::sqrt(W * H / n);
    double temp = (start_temp > 0) ? start_temp : W / 10.0;
    double cool = temp / iterations;
    std::vector<V2> disp(n);

    for (int iter = 0; iter < iterations; ++iter) {
        for (auto& d : disp) d = {0.0, 0.0};

        // Repulsion between all pairs
        for (int v = 0; v < n; ++v) {
            for (int u = 0; u < n; ++u) {
                if (u == v) continue;
                double dx = pos[v].x - pos[u].x;
                double dy = pos[v].y - pos[u].y;
                double d  = std::max(std::hypot(dx, dy), 0.01);
                double f  = k * k / d;
                disp[v].x += dx / d * f;
                disp[v].y += dy / d * f;
            }
        }

        // Attraction along edges
        for (auto& [u, v] : edges) {
            double dx = pos[v].x - pos[u].x;
            double dy = pos[v].y - pos[u].y;
            double d  = std::max(std::hypot(dx, dy), 0.01);
            double f  = d * d / k;
            double fx = dx / d * f, fy = dy / d * f;
            disp[v].x -= fx; disp[v].y -= fy;
            disp[u].x += fx; disp[u].y += fy;
        }

        // Apply
        for (int v = 0; v < n; ++v) {
            double d    = std::max(std::hypot(disp[v].x, disp[v].y), 0.01);
            double move = std::min(d, temp);
            pos[v].x += disp[v].x / d * move;
            pos[v].y += disp[v].y / d * move;
        }

        temp = std::max(temp - cool, 0.01);
    }
    return pos;
}

// ---------------------------------------------------------------------------
// Normalize positions to fit canvas [pad, W-pad] x [pad, H-pad]
// ---------------------------------------------------------------------------
static void normalise(std::vector<V2>& pos, double W, double H, double pad)
{
    if (pos.empty()) return;
    double minx = pos[0].x, maxx = pos[0].x;
    double miny = pos[0].y, maxy = pos[0].y;
    for (auto& p : pos) {
        minx = std::min(minx, p.x); maxx = std::max(maxx, p.x);
        miny = std::min(miny, p.y); maxy = std::max(maxy, p.y);
    }
    double rx = maxx - minx, ry = maxy - miny;
    double sx = (rx > 1) ? (W - 2 * pad) / rx : 1.0;
    double sy = (ry > 1) ? (H - 2 * pad) / ry : 1.0;
    double s  = std::min(sx, sy);
    double ox = (W - s * rx) / 2.0;
    double oy = (H - s * ry) / 2.0;
    for (auto& p : pos) {
        p.x = ox + (p.x - minx) * s;
        p.y = oy + (p.y - miny) * s;
    }
}

// ---------------------------------------------------------------------------
// Overlap resolution — push pairs apart until boundaries are min_gap apart.
// min_dist = 2*radius + min_gap (center-to-center).
// ---------------------------------------------------------------------------
static void resolveOverlaps(std::vector<V2>& pos, double min_dist,
                             double W, double H, double pad,
                             int max_iter = 500)
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
                double need = (min_dist - d) * 0.5 + 0.5;
                if (d < 0.01) { dx = 1.0; dy = 0.0; d = 1.0; }
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
// Two-level block-aware layout
// Level 1: FR on block graph  →  block centres
// Level 2: circular arrangement within each block
// Level 3: short FR refinement on full graph
// ---------------------------------------------------------------------------
std::vector<V2> GraphVisualizer::computeLayout(
    const ReliabilityGraph& g, double W, double H, int iterations)
{
    const int n = static_cast<int>(g.numVertices());
    if (n == 0) return {};
    if (n == 1) return {{W / 2.0, H / 2.0}};

    // Build full-graph edge list
    std::vector<std::pair<int,int>> all_edges;
    for (int u = 0; u < n; ++u)
        for (size_t i = g.kao_[u]; i < g.kao_[u + 1]; ++i) {
            int v = g.fo_[i];
            if (v > u) all_edges.push_back({u, v});
        }

    auto decomp     = g.decomposeIntoBlocks();
    int  num_blocks = decomp.back();

    // ---- Single-block fallback: plain FR ----
    if (num_blocks <= 1) {
        std::vector<V2> init(n);
        double cx = W / 2, cy = H / 2, r0 = std::min(W, H) * 0.4;
        for (int i = 0; i < n; ++i) {
            double a = 2.0 * M_PI * i / n;
            init[i] = {cx + r0 * std::cos(a), cy + r0 * std::sin(a)};
        }
        auto pos = runFR(n, all_edges, init, W, H, iterations);
        normalise(pos, W, H, 50.0);
        return pos;
    }

    // ---- Build block membership ----
    std::vector<std::vector<int>> block_verts(num_blocks + 1); // 1-indexed
    std::vector<std::vector<int>> vertex_blocks(n);
    std::vector<bool> is_ap(n, false);

    for (int v = 0; v < n; ++v) {
        auto bs = g.getBlocksContainingVertex(v, decomp);
        vertex_blocks[v] = bs;
        for (int b : bs) block_verts[b].push_back(v);
        if (bs.size() > 1) is_ap[v] = true;
    }

    // ---- Level 1: FR on block graph ----
    // Build block-level edges via shared APs
    std::set<std::pair<int,int>> beset;
    for (int v = 0; v < n; ++v) {
        if (!is_ap[v]) continue;
        auto& bs = vertex_blocks[v];
        for (size_t i = 0; i < bs.size(); ++i)
            for (size_t j = i + 1; j < bs.size(); ++j) {
                int a = bs[i] - 1, b = bs[j] - 1; // 0-indexed
                if (a > b) std::swap(a, b);
                beset.insert({a, b});
            }
    }
    std::vector<std::pair<int,int>> block_edges(beset.begin(), beset.end());

    // Initial block positions: circle
    std::vector<V2> binit(num_blocks);
    {
        double cx = W / 2, cy = H / 2;
        double r0 = std::min(W, H) * 0.38;
        for (int b = 0; b < num_blocks; ++b) {
            double a = 2.0 * M_PI * b / num_blocks;
            binit[b] = {cx + r0 * std::cos(a), cy + r0 * std::sin(a)};
        }
    }

    auto block_centers = runFR(num_blocks, block_edges, binit, W, H, iterations / 2);
    normalise(block_centers, W, H, 80.0);

    // ---- Level 2: arrange vertices within each block ----
    const double vr      = 14.0; // vertex radius
    const double min_gap = 6.0;
    const double step    = 2.0 * vr + min_gap; // centre-to-centre step

    std::vector<V2> pos(n, {0.0, 0.0});

    for (int b = 1; b <= num_blocks; ++b) {
        std::vector<int> reg; // non-AP vertices in this block
        for (int v : block_verts[b])
            if (!is_ap[v]) reg.push_back(v);

        int m = static_cast<int>(reg.size());
        // Radius of the local circle: enough so adjacent vertices don't overlap
        double local_r = (m <= 1) ? 0.0
                                  : (step / (2.0 * std::sin(M_PI / m)));
        local_r = std::max(local_r, step * 0.5); // minimum for m=1 case

        V2 bc = block_centers[b - 1];
        for (int i = 0; i < m; ++i) {
            double a = 2.0 * M_PI * i / std::max(m, 1);
            pos[reg[i]] = {bc.x + local_r * std::cos(a),
                           bc.y + local_r * std::sin(a)};
        }
    }

    // APs: weighted average of adjacent block centres
    for (int v = 0; v < n; ++v) {
        if (!is_ap[v]) continue;
        double ax = 0, ay = 0;
        for (int b : vertex_blocks[v]) {
            ax += block_centers[b - 1].x;
            ay += block_centers[b - 1].y;
        }
        int cnt = static_cast<int>(vertex_blocks[v].size());
        pos[v] = {ax / cnt, ay / cnt};
    }

    // ---- Level 3: short FR refinement on full graph ----
    // Use a small starting temperature so the block structure is preserved.
    double fine_temp = std::min(W, H) / 20.0;
    pos = runFR(n, all_edges, pos, W, H, iterations / 4, fine_temp);
    normalise(pos, W, H, 50.0);

    // ---- Resolve any remaining overlaps ----
    double min_dist = 2.0 * vr + min_gap;
    resolveOverlaps(pos, min_dist, W, H, vr + min_gap, 600);

    return pos;
}

// ---------------------------------------------------------------------------
// SVG helpers
// ---------------------------------------------------------------------------
static std::string fmtProb(double p) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(2) << p;
    return ss.str();
}

// Pastel block background colours (cycle through them)
static const char* BLOCK_FILL[] = {
    "#dce8f5", "#d5f0e0", "#f5e8dc", "#ece0f5",
    "#f5f0d5", "#d5f0f0", "#f5d5e8", "#e8f5d5",
    "#f5d5d5", "#d5d5f5"
};
static const int NUM_BLOCK_COLORS = 10;

// ---------------------------------------------------------------------------
// SVG export
// ---------------------------------------------------------------------------
void GraphVisualizer::exportSVG(const ReliabilityGraph& g,
                                const std::string& path,
                                const Options& opts)
{
    const int n = static_cast<int>(g.numVertices());
    if (n == 0) throw std::invalid_argument("Graph has no vertices");

    const double vr      = opts.vertex_radius;
    const double min_gap = 6.0;

    // Auto-scale canvas when there are many vertices
    double W = opts.width, H = opts.height;
    if (n > 1) {
        double min_dist = 2.0 * vr + min_gap;
        double needed   = static_cast<double>(n) * min_dist * min_dist * 2.0;
        if (needed > W * H) {
            double s = std::sqrt(needed / (W * H));
            W *= s; H *= s;
        }
    }

    auto pos = computeLayout(g, W, H, opts.iterations);

    // Block membership for background drawing
    auto decomp     = g.decomposeIntoBlocks();
    int  num_blocks = decomp.back();

    std::vector<std::vector<int>> block_verts(num_blocks + 1);
    std::vector<bool> is_ap(n, false);
    for (int v = 0; v < n; ++v) {
        auto bs = g.getBlocksContainingVertex(v, decomp);
        for (int b : bs) block_verts[b].push_back(v);
        if (bs.size() > 1) is_ap[v] = true;
    }

    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot open output file: " + path);

    int iW = static_cast<int>(W), iH = static_cast<int>(H);
    f << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
      << "<svg xmlns=\"http://www.w3.org/2000/svg\""
      << " width=\"" << iW << "\" height=\"" << iH << "\""
      << " viewBox=\"0 0 " << iW << " " << iH << "\">\n"
      << "<rect width=\"100%\" height=\"100%\" fill=\"#f8f9fa\"/>\n";

    // --- Block backgrounds ---
    if (num_blocks > 1) {
        f << "<g opacity=\"0.55\">\n";
        for (int b = 1; b <= num_blocks; ++b) {
            auto& verts = block_verts[b];
            if (verts.empty()) continue;
            // Centroid
            double cx = 0, cy = 0;
            for (int v : verts) { cx += pos[v].x; cy += pos[v].y; }
            cx /= verts.size(); cy /= verts.size();
            // Max radius from centroid
            double maxr = 0;
            for (int v : verts)
                maxr = std::max(maxr, std::hypot(pos[v].x - cx, pos[v].y - cy));
            double bgr = maxr + vr + 6.0;
            const char* fill = BLOCK_FILL[(b - 1) % NUM_BLOCK_COLORS];
            f << "  <circle cx=\"" << cx << "\" cy=\"" << cy
              << "\" r=\"" << bgr << "\""
              << " fill=\"" << fill << "\" stroke=\"#ccc\" stroke-width=\"1\"/>\n";
        }
        f << "</g>\n";
    }

    // --- Edges (trimmed to vertex boundary) ---
    f << "<g stroke=\"#8e9aaf\" stroke-width=\"1.5\">\n";
    for (int u = 0; u < n; ++u) {
        for (size_t i = g.kao_[u]; i < g.kao_[u + 1]; ++i) {
            int v = g.fo_[i];
            if (v <= u) continue;
            double dx   = pos[v].x - pos[u].x;
            double dy   = pos[v].y - pos[u].y;
            double dist = std::hypot(dx, dy);
            if (dist < 0.01) continue;
            double nx = dx / dist, ny = dy / dist;
            double x1 = pos[u].x + nx * vr, y1 = pos[u].y + ny * vr;
            double x2 = pos[v].x - nx * vr, y2 = pos[v].y - ny * vr;
            f << "  <line x1=\"" << x1 << "\" y1=\"" << y1 << "\""
              << " x2=\"" << x2 << "\" y2=\"" << y2 << "\"/>\n";
            if (opts.show_probs) {
                double mx = (pos[u].x + pos[v].x) / 2.0;
                double my = (pos[u].y + pos[v].y) / 2.0;
                f << "  <text x=\"" << mx << "\" y=\"" << my
                  << "\" font-size=\"9\" fill=\"#6c757d\""
                  << " text-anchor=\"middle\" dominant-baseline=\"middle\">"
                  << fmtProb(g.p_array_[i]) << "</text>\n";
            }
        }
    }
    f << "</g>\n";

    // --- Vertices ---
    int fs = std::max(8, static_cast<int>(vr) - 3);
    f << "<g font-family=\"sans-serif\" font-size=\"" << fs << "\" font-weight=\"bold\">\n";
    for (int v = 0; v < n; ++v) {
        std::string fill, stroke = "#2c3e50";
        if (v == opts.source)       { fill = "#27ae60"; }   // source: green
        else if (v == opts.target)  { fill = "#e74c3c"; }   // target: red
        else if (is_ap[v])          { fill = "#f39c12"; stroke = "#8a5500"; } // AP: orange
        else                        { fill = "#4a90d9"; }   // regular: blue

        f << "  <circle cx=\"" << pos[v].x << "\" cy=\"" << pos[v].y
          << "\" r=\"" << vr << "\""
          << " fill=\"" << fill << "\""
          << " stroke=\"" << stroke << "\" stroke-width=\"1.5\"/>\n"
          << "  <text x=\"" << pos[v].x << "\" y=\"" << pos[v].y
          << "\" text-anchor=\"middle\" dominant-baseline=\"middle\""
          << " fill=\"white\">" << v << "</text>\n";
    }
    f << "</g>\n";

    // --- Legend ---
    f << "<g font-family=\"sans-serif\" font-size=\"11\" fill=\"#495057\">\n"
      << "  <text x=\"10\" y=\"18\">V=" << n
      << "  E=" << (g.fo_.size() / 2) << "</text>\n";
    int ly = 34;
    if (opts.source >= 0) {
        f << "  <circle cx=\"10\" cy=\"" << ly << "\" r=\"5\" fill=\"#27ae60\"/>"
          << "<text x=\"20\" y=\"" << (ly + 4) << "\">s=" << opts.source << "</text>\n";
        ly += 16;
    }
    if (opts.target >= 0) {
        f << "  <circle cx=\"10\" cy=\"" << ly << "\" r=\"5\" fill=\"#e74c3c\"/>"
          << "<text x=\"20\" y=\"" << (ly + 4) << "\">t=" << opts.target << "</text>\n";
        ly += 16;
    }
    if (num_blocks > 1) {
        f << "  <circle cx=\"10\" cy=\"" << ly << "\" r=\"5\" fill=\"#f39c12\"/>"
          << "<text x=\"20\" y=\"" << (ly + 4) << "\">точка сочленения</text>\n";
    }
    f << "</g>\n";

    f << "</svg>\n";
}

// ---------------------------------------------------------------------------
// DOT export
// ---------------------------------------------------------------------------
void GraphVisualizer::exportDot(const ReliabilityGraph& g,
                                const std::string& path,
                                const Options& opts)
{
    const int n = static_cast<int>(g.numVertices());

    auto decomp     = g.decomposeIntoBlocks();
    int  num_blocks = decomp.back();
    std::vector<bool> is_ap(n, false);
    for (int v = 0; v < n; ++v)
        if (g.getBlocksContainingVertex(v, decomp).size() > 1) is_ap[v] = true;

    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot open output file: " + path);

    f << "graph G {\n"
      << "  overlap=false;\n"
      << "  splines=true;\n"
      << "  node [shape=circle, fixedsize=true, width=0.4, fontsize=10];\n\n";

    for (int v = 0; v < n; ++v) {
        if (v == opts.source)
            f << "  " << v << " [style=filled, fillcolor=\"#27ae60\", fontcolor=white];\n";
        else if (v == opts.target)
            f << "  " << v << " [style=filled, fillcolor=\"#e74c3c\", fontcolor=white];\n";
        else if (is_ap[v])
            f << "  " << v << " [style=filled, fillcolor=\"#f39c12\", fontcolor=white];\n";
    }

    // Group vertices by block using subgraph clusters
    if (num_blocks > 1) {
        for (int b = 1; b <= num_blocks; ++b) {
            auto verts = g.getBlocksContainingVertex(0, decomp); // dummy
            // Collect non-AP vertices of this block for cluster grouping
            std::vector<int> reg;
            for (int v = 0; v < n; ++v) {
                auto bs = g.getBlocksContainingVertex(v, decomp);
                if (bs.size() == 1 && bs[0] == b) reg.push_back(v);
            }
            if (reg.empty()) continue;
            f << "  subgraph cluster_" << b << " {\n"
              << "    style=filled; color=\"" << BLOCK_FILL[(b-1) % NUM_BLOCK_COLORS] << "\";\n";
            for (int v : reg) f << "    " << v << ";\n";
            f << "  }\n";
        }
    }
    f << "\n";

    for (int u = 0; u < n; ++u)
        for (size_t i = g.kao_[u]; i < g.kao_[u + 1]; ++i) {
            int v = g.fo_[i];
            if (v <= u) continue;
            f << "  " << u << " -- " << v;
            if (opts.show_probs)
                f << " [label=\"" << fmtProb(g.p_array_[i]) << "\"]";
            f << ";\n";
        }

    f << "}\n";
}

} // namespace graph_reliability
