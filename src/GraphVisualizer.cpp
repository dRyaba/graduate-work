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
// Level 1: FR on block graph (3W×3H virtual) → push apart → normalise to canvas
// Level 2: non-APs on local circles, APs at centroid of their block centres
// Level 3: overlap resolution only (no FR — preserves block groupings)
// ---------------------------------------------------------------------------
std::vector<V2> GraphVisualizer::computeLayout(
    const ReliabilityGraph& g, double W, double H, int iterations)
{
    const int n = static_cast<int>(g.numVertices());
    if (n == 0) return {};
    if (n == 1) return {{W / 2.0, H / 2.0}};

    std::vector<std::pair<int,int>> all_edges;
    for (int u = 0; u < n; ++u)
        for (size_t i = g.kao_[u]; i < g.kao_[u + 1]; ++i) {
            int v = g.fo_[i];
            if (v > u) all_edges.push_back({u, v});
        }

    auto decomp     = g.decomposeIntoBlocks();
    int  num_blocks = decomp.back();

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

    const double vr      = 14.0;
    const double min_gap = 6.0;
    const double step    = 2.0 * vr + min_gap; // 34 px — centre-to-centre budget

    // ---- Block membership ----
    std::vector<std::vector<int>> block_verts(num_blocks + 1); // 1-indexed
    std::vector<std::vector<int>> vertex_blocks(n);
    std::vector<bool>             is_ap(n, false);

    for (int v = 0; v < n; ++v) {
        auto bs = g.getBlocksContainingVertex(v, decomp);
        vertex_blocks[v] = bs;
        for (int b : bs) block_verts[b].push_back(v);
        if (bs.size() > 1) is_ap[v] = true;
    }

    // ---- Block-graph edges (via shared APs) ----
    std::set<std::pair<int,int>> beset;
    for (int v = 0; v < n; ++v) {
        if (!is_ap[v]) continue;
        auto& bs = vertex_blocks[v];
        for (size_t i = 0; i < bs.size(); ++i)
            for (size_t j = i + 1; j < bs.size(); ++j) {
                int a = bs[i] - 1, b = bs[j] - 1;
                if (a > b) std::swap(a, b);
                beset.insert({a, b});
            }
    }
    std::vector<std::pair<int,int>> block_edges(beset.begin(), beset.end());

    // ---- Level 1: FR on block graph in virtual 3W×3H canvas ----
    double BW = W * 3.0, BH = H * 3.0;
    std::vector<V2> binit(num_blocks);
    {
        double cx = BW / 2, cy = BH / 2;
        double r0 = std::min(BW, BH) * 0.42;
        for (int b = 0; b < num_blocks; ++b) {
            double a = 2.0 * M_PI * b / num_blocks;
            binit[b] = {cx + r0 * std::cos(a), cy + r0 * std::sin(a)};
        }
    }
    auto block_centers = runFR(num_blocks, block_edges, binit, BW, BH,
                               iterations * 2 / 3);

    // Compute per-block circle radius in virtual space (scaled so after
    // normalise to W×H the vertex circles have correct step = 34 px).
    double scale_v = std::sqrt(BW * BH / (W * H)); // ≈ 3
    std::vector<double> bR(num_blocks); // circle radius per block in virtual space
    for (int b = 1; b <= num_blocks; ++b) {
        std::vector<int> reg;
        for (int v : block_verts[b])
            if (!is_ap[v]) reg.push_back(v);
        int m = static_cast<int>(reg.size());
        double lr = (m <= 1) ? step * 0.5 : step / (2.0 * std::sin(M_PI / m));
        lr = std::max(lr, step * 0.5);
        bR[b - 1] = lr * scale_v; // in virtual space
    }

    // Push block circles apart in virtual space so they don't overlap
    {
        double max_r = *std::max_element(bR.begin(), bR.end());
        double bc_min = 2.0 * max_r + min_gap * scale_v * 4.0;
        for (int iter = 0; iter < 600; ++iter) {
            bool any = false;
            for (int i = 0; i < num_blocks; ++i) {
                for (int j = i + 1; j < num_blocks; ++j) {
                    double dx = block_centers[j].x - block_centers[i].x;
                    double dy = block_centers[j].y - block_centers[i].y;
                    double d  = std::hypot(dx, dy);
                    double need = bR[i] + bR[j] + min_gap * scale_v * 4.0;
                    if (d >= need) continue;
                    any = true;
                    double push = (need - d) * 0.51;
                    if (d < 0.01) { dx = 1.0; dy = 0.0; d = 1.0; }
                    double nx = dx / d, ny = dy / d;
                    block_centers[i].x -= nx * push;
                    block_centers[i].y -= ny * push;
                    block_centers[j].x += nx * push;
                    block_centers[j].y += ny * push;
                }
            }
            if (!any) break;
        }
    }

    // Normalise block centres to canvas
    normalise(block_centers, W, H, 70.0);

    // ---- Level 2: place vertices ----
    std::vector<V2> pos(n, {0.0, 0.0});

    // Non-APs on local circles (in canvas space, step = 34 px)
    for (int b = 1; b <= num_blocks; ++b) {
        std::vector<int> reg;
        for (int v : block_verts[b])
            if (!is_ap[v]) reg.push_back(v);
        int m = static_cast<int>(reg.size());
        double lr = (m <= 1) ? step * 0.5 : step / (2.0 * std::sin(M_PI / m));
        lr = std::max(lr, step * 0.5);

        V2 bc = block_centers[b - 1];
        for (int i = 0; i < m; ++i) {
            double a = 2.0 * M_PI * i / std::max(m, 1);
            pos[reg[i]] = {bc.x + lr * std::cos(a), bc.y + lr * std::sin(a)};
        }
    }

    // APs at centroid of their block centres (block centres now in canvas space)
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

    // Final normalise + overlap resolution (no FR — preserves block groupings)
    normalise(pos, W, H, 50.0);
    double min_dist = 2.0 * vr + min_gap;
    resolveOverlaps(pos, min_dist, W, H, vr + min_gap, 800);

    // ---- Edge-vertex repulsion: push non-endpoint vertices off edges ----
    // Runs alternately with resolveOverlaps to clear the "cobweb" look where
    // edges visually cross foreign vertices.
    const double ev_clear = vr * 1.8; // perpendicular clearance target
    const double pad = vr + min_gap;
    for (int outer = 0; outer < 40; ++outer) {
        bool any = false;
        for (auto& [u, v] : all_edges) {
            double dx = pos[v].x - pos[u].x;
            double dy = pos[v].y - pos[u].y;
            double L  = std::hypot(dx, dy);
            if (L < 0.01) continue;
            double nx = dx / L, ny = dy / L;
            for (int w = 0; w < n; ++w) {
                if (w == u || w == v) continue;
                double wx = pos[w].x - pos[u].x;
                double wy = pos[w].y - pos[u].y;
                double t  = wx * nx + wy * ny;
                if (t < vr * 0.5 || t > L - vr * 0.5) continue;
                double perp = -wx * ny + wy * nx;
                double ap = std::fabs(perp);
                if (ap >= ev_clear) continue;
                any = true;
                double need = (ev_clear - ap) * 0.5 + 0.3;
                double sign = (perp >= 0) ? 1.0 : -1.0;
                // Push w away; push endpoints slightly toward opposite side
                pos[w].x = std::clamp(pos[w].x + (-ny) * sign * need * 0.7,
                                      pad, W - pad);
                pos[w].y = std::clamp(pos[w].y + ( nx) * sign * need * 0.7,
                                      pad, H - pad);
                pos[u].x = std::clamp(pos[u].x - (-ny) * sign * need * 0.15,
                                      pad, W - pad);
                pos[u].y = std::clamp(pos[u].y - ( nx) * sign * need * 0.15,
                                      pad, H - pad);
                pos[v].x = std::clamp(pos[v].x - (-ny) * sign * need * 0.15,
                                      pad, W - pad);
                pos[v].y = std::clamp(pos[v].y - ( nx) * sign * need * 0.15,
                                      pad, H - pad);
            }
        }
        // Re-resolve vertex overlaps after each pass
        resolveOverlaps(pos, min_dist, W, H, pad, 100);
        if (!any) break;
    }

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

// 20-colour palette for block graph colouring (avoids source green / target red).
// Adjacent blocks always get different colours via greedy graph colouring.
static const char* BLOCK_PALETTE[] = {
    "#3a86ff", // blue
    "#8338ec", // violet
    "#fb8500", // orange
    "#023047", // navy
    "#219ebc", // sky blue
    "#8ecae6", // light blue
    "#6a4c93", // purple
    "#f4a261", // peach
    "#264653", // dark teal
    "#2a9d8f", // teal
    "#e9c46a", // gold
    "#457b9d", // steel blue
    "#a8dadc", // pale cyan
    "#c77dff", // lavender
    "#780000", // dark red (distinct from #e74c3c target)
    "#606c38", // olive
    "#bc6c25", // brown
    "#48cae4", // cerulean
    "#7209b7", // deep purple
    "#4cc9f0", // cyan
};
static const int NUM_BLOCK_PALETTE = 20;

// Used for DOT cluster backgrounds (still pastel fill there)
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

    // Auto-scale: give every vertex a (2.5 × step) square of private space.
    // For n=103 this expands the canvas to ~1600×1250, eliminating the
    // "cobweb" look caused by edges passing through foreign vertices.
    double W = opts.width, H = opts.height;
    if (n > 1) {
        double step_sz = 2.0 * vr + min_gap;
        double cell    = step_sz * 4.0; // per-vertex area side
        double needed  = static_cast<double>(n) * cell * cell;
        if (needed > W * H) {
            double s = std::sqrt(needed / (W * H));
            W *= s; H *= s;
        }
    }

    auto pos = computeLayout(g, W, H, opts.iterations);

    // ---- Block membership + greedy graph colouring ----
    auto decomp     = g.decomposeIntoBlocks();
    int  num_blocks = decomp.back();

    std::vector<int>  vertex_block(n, 0); // primary block index (0-based)
    std::vector<bool> is_ap(n, false);
    for (int v = 0; v < n; ++v) {
        auto bs = g.getBlocksContainingVertex(v, decomp);
        if (!bs.empty()) vertex_block[v] = bs[0] - 1;
        if (bs.size() > 1) is_ap[v] = true;
    }

    // Build block adjacency list for graph colouring
    std::vector<std::vector<int>> block_adj(num_blocks);
    for (int v = 0; v < n; ++v) {
        if (!is_ap[v]) continue;
        auto bs = g.getBlocksContainingVertex(v, decomp);
        for (size_t i = 0; i < bs.size(); ++i)
            for (size_t j = i + 1; j < bs.size(); ++j) {
                int a = bs[i] - 1, b = bs[j] - 1;
                block_adj[a].push_back(b);
                block_adj[b].push_back(a);
            }
    }

    // Greedy colouring: adjacent blocks get different palette indices
    std::vector<int> block_color(num_blocks, 0);
    for (int b = 0; b < num_blocks; ++b) {
        std::set<int> used;
        for (int nb : block_adj[b])
            if (nb < b) used.insert(block_color[nb]);
        int c = 0;
        while (used.count(c)) ++c;
        block_color[b] = c % NUM_BLOCK_PALETTE;
    }

    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot open output file: " + path);

    int iW = static_cast<int>(W), iH = static_cast<int>(H);
    f << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
      << "<svg xmlns=\"http://www.w3.org/2000/svg\""
      << " width=\"" << iW << "\" height=\"" << iH << "\""
      << " viewBox=\"0 0 " << iW << " " << iH << "\">\n"
      << "<rect width=\"100%\" height=\"100%\" fill=\"#f8f9fa\"/>\n";

    // --- Edges (straight, trimmed to vertex boundary) ---
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
    // Colour: source=green, target=red, others=block palette colour.
    // APs are indicated by a white outer ring (not a different fill colour).
    int fs = std::max(8, static_cast<int>(vr) - 3);
    f << "<g font-family=\"sans-serif\" font-size=\"" << fs << "\" font-weight=\"bold\">\n";
    for (int v = 0; v < n; ++v) {
        std::string fill, stroke = "#1a1a2e";
        if (v == opts.source) {
            fill = "#27ae60"; stroke = "#145a32";
        } else if (v == opts.target) {
            fill = "#e74c3c"; stroke = "#7b241c";
        } else {
            fill = BLOCK_PALETTE[block_color[vertex_block[v]]];
        }

        // AP indicator: white outer ring so the block colour stays visible
        if (is_ap[v] && v != opts.source && v != opts.target) {
            f << "  <circle cx=\"" << pos[v].x << "\" cy=\"" << pos[v].y
              << "\" r=\"" << (vr + 3.0) << "\""
              << " fill=\"white\" stroke=\"" << fill
              << "\" stroke-width=\"2\"/>\n";
        }

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
        f << "  <circle cx=\"10\" cy=\"" << ly << "\" r=\"8\" fill=\"#3a86ff\""
          << " stroke=\"#1a1a2e\" stroke-width=\"1.5\"/>"
          << "<circle cx=\"10\" cy=\"" << ly << "\" r=\"11\" fill=\"white\""
          << " stroke=\"#3a86ff\" stroke-width=\"2\"/>"
          << "<text x=\"26\" y=\"" << (ly + 4) << "\">точка сочленения</text>\n";
        ly += 16;
        f << "  <text x=\"10\" y=\"" << (ly + 4) << "\">цвет = блок (жадная раскраска)</text>\n";
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
