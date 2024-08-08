#pragma once
#include <bits/stdc++.h>
#include <boost/graph/maximum_weighted_matching.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/property_map/function_property_map.hpp>

#include <fmt/core.h>
#include <fmt/ranges.h>

using weight_t = double;

using std::pair;
using std::vector;
using std::max;
using std::min;
using boost::num_vertices;
using boost::add_edge;
using boost::maximum_weighted_matching;
using boost::matching_weight_sum;
using boost::prim_minimum_spanning_tree;
using boost::print_graph;
using boost::filtered_graph;

using edge_property_t = boost::property<boost::edge_weight_t, weight_t>;
using graph_t = boost::adjacency_matrix<boost::undirectedS, boost::no_property, edge_property_t>;
using vertex_t = graph_t::vertex_descriptor;

using pvv = pair<vertex_t, vertex_t>;

vector<pvv> collect_tree_edges(const vector<vertex_t>& pre) {
    vector<pvv> res;
    for (vertex_t i = 0; i < pre.size(); ++i) {
        if (pre[i] != i) {
            res.emplace_back(pre[i], i);
        }
    }
    return res;
}

vector<size_t> calc_deg(size_t n, const vector<pvv>& edges) {
    vector<vertex_t> deg(n, 0);
    for (auto p : edges) {
        deg[p.first]++;
        deg[p.second]++;
    }
    return deg;
}

vector<vertex_t> collect_odd_deg_vertices(const vector<vertex_t>& deg) {
    const size_t n = deg.size();
    vector<vertex_t> res;
    for (vertex_t i = 0; i < n; ++i) {
        if (deg[i] % 2 == 1) 
            res.push_back(i);
    }
    return res;
}

graph_t induce(const graph_t& g, const vector<vertex_t>& id) {
    graph_t h(id.size());
    for (size_t i = 0; i < id.size(); ++i) {
        for (size_t j = i + 1; j < id.size(); ++j) {
            auto e = g.get_edge(id[i], id[j]);
            if (!e.first) continue;
            add_edge(i, j, e.second, h);
        }
    }
    return h;
}

vector<pvv> collect_matching_edges(const vector<vertex_t>& mate) {
    vector<pvv> res;
    const int n = mate.size();
    for (size_t i = 0; i < n; ++i) {
        if (mate[i] < i) {
            res.emplace_back(mate[i], i);
        }
    }
    return res;
}

vector<pvv> restore(const vector<pvv>& edges, const vector<vertex_t>& id) {
    vector<pvv> res;
    res.reserve(edges.size());
    for (pvv p : edges) {
        res.emplace_back(id[p.first], id[p.second]);
    }
    return res;
}

vector<vertex_t> eulerian_circuit(size_t n, const vector<pvv>& edges, vertex_t s) {
    vector<vector<pvv>> g(n);
    for (pvv p : edges) {
        vertex_t u = p.first, v = p.second;
        g[u].push_back({ v, g[v].size() });
        g[v].push_back({ u, g[u].size() - 1 });
    }

    vector<vertex_t> stk, res;
    stk.push_back(s);
    while (!stk.empty()) {
        int u = stk.back();
        if (g[u].empty()) {
            res.push_back(u);
            stk.pop_back();
        }
        else {
            pvv p = g[u].back();
            int v = p.first, i = p.second;
            if (v != u) {
                if (i + 1 != g[v].size()) {
                    pvv& q = g[v].back();
                    g[q.first][q.second].second = i;
                    swap(g[v][i], q);
                }
                g[u].pop_back();
                g[v].pop_back();
                stk.push_back(v);
            }
            else {
                g[u].pop_back();
                stk.push_back(u);
            }
        }
    }

    return res;
}

void rotate_eul_tour(vector<vertex_t>& tour) {
    rotate(tour.begin(), tour.begin() + 1, tour.end());
    tour.back() = tour.front();
}

vector<vertex_t> filter_eul_tour(size_t n, const vector<vertex_t>& tour) {
    vector<bool> appeared(n, false);
    vector<vertex_t> res;
    for (vertex_t v : tour) {
        if (appeared[v]) continue;
        appeared[v] = true;
        res.push_back(v);
    }
    vertex_t v = res.front();
    res.push_back(v);
    return res;
}

weight_t tsp_tour_weight(const graph_t& g, const vector<vertex_t>& tour) {
    const size_t n = num_vertices(g);
    weight_t sum = 0;
    for (size_t i = 0; i < n; ++i) {
        size_t u = tour[i], v = tour[i + 1];
        sum += g.get_edge(u,v).second.m_value;
    }
    return sum;
}

vector<vertex_t> minimum_weighted_matching(const graph_t& g) {
    const size_t n = num_vertices(g);
    vector<vertex_t> mate(n);

    weight_t max_w = 0;

    for (size_t i = 0; i < n; ++i)
        for (size_t j = i + 1; j < n; ++j) 
            max_w = max(max_w, g.get_edge(i, j).second.m_value);

    graph_t h(n);
    for (size_t i = 0; i < n; ++i)
        for (size_t j = i + 1; j < n; ++j) 
            add_edge(i, j, max_w - g.get_edge(i, j).second.m_value, h);

    maximum_weighted_matching(g, mate.data());
    return mate;
}