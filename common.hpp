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
constexpr weight_t inf = 1e20;

using boost::add_edge;
using boost::edge_weight;
using boost::edges;
using boost::filtered_graph;
using boost::get;
using boost::matching_weight_sum;
using boost::maximum_weighted_matching;
using boost::num_vertices;
using boost::prim_minimum_spanning_tree;
using boost::print_graph;
using boost::source;
using boost::target;
using std::make_pair;
using std::lower_bound;
using std::unordered_map;
using std::map;
using std::max;
using std::min;
using std::move;
using std::pair;
using std::rotate;
using std::swap;
using std::vector;
using std::unordered_set;
using std::remove_if;

using edge_property_t = boost::property<boost::edge_weight_t, weight_t>;
using graph_t = boost::adjacency_matrix<boost::undirectedS, boost::no_property, edge_property_t>;
using vertex_t = graph_t::vertex_descriptor;

using pvv = pair<vertex_t, vertex_t>;

void normalize(pvv& p) {
    if (p.first > p.second)
        swap(p.first, p.second);
}

struct vmap : vector<vertex_t> {
    void build() {
        sort(begin(), end());
        erase(unique(begin(), end()), end());
    }
    vertex_t id(vertex_t v) {
        return lower_bound(begin(), end(), v) - begin();
    }
};

struct hash_pvv {
    size_t operator()(pvv p) const {
        return ((size_t)(p.first) << 32) | p.second;
    }
};

pair<weight_t, pvv> maximum_edge_weight(const graph_t &g) {
    pair<weight_t, pvv> result(
        std::numeric_limits<weight_t>::min(),
        pvv(graph_t::null_vertex(), graph_t::null_vertex()));
    for (auto [it, end] = edges(g); it != end; ++it) {
        if (!it->exists())
            continue;
        auto src = boost::source(*it, g);
        auto dst = boost::target(*it, g);
        auto w = boost::get(edge_weight, g, *it);
        result = max(result, make_pair(w, pvv(src, dst)));
    }
    return result;
}

vector<pvv> collect_tree_edges(const vector<vertex_t> &pre) {
    vector<pvv> res;
    for (vertex_t i = 0; i < pre.size(); ++i) {
        if (pre[i] != i) {
            res.emplace_back(pre[i], i);
        }
    }
    return res;
}

vector<size_t> calc_deg(size_t n, const vector<pvv> &edges) {
    vector<vertex_t> deg(n, 0);
    for (auto p : edges) {
        deg[p.first]++;
        deg[p.second]++;
    }
    return deg;
}

vector<vertex_t> collect_odd_deg_vertices(const vector<vertex_t> &deg) {
    const size_t n = deg.size();
    vector<vertex_t> res;
    for (vertex_t i = 0; i < n; ++i) {
        if (deg[i] % 2 == 1)
            res.push_back(i);
    }
    return res;
}

graph_t induce(const graph_t &g, const vector<vertex_t> &ids) {
    graph_t h(ids.size());
    for (size_t i = 0; i < ids.size(); ++i) {
        for (size_t j = i + 1; j < ids.size(); ++j) {
            auto e = g.get_edge(ids[i], ids[j]);
            if (!e.first)
                continue;
            add_edge(i, j, e.second, h);
        }
    }
    return h;
}

vector<pvv> collect_matching_edges(const vector<vertex_t> &mate) {
    vector<pvv> res;
    const int n = mate.size();
    for (size_t i = 0; i < n; ++i) {
        if (mate[i] < i) {
            res.emplace_back(mate[i], i);
        }
    }
    return res;
}

vector<pvv> restore(const vector<pvv> &edges, const vector<vertex_t> &ids) {
    vector<pvv> res;
    res.reserve(edges.size());
    for (pvv p : edges) {
        res.emplace_back(ids[p.first], ids[p.second]);
    }
    return res;
}

vector<vertex_t> restore(const vector<vertex_t>& tour, const vector<vertex_t> &ids) {
    vector<vertex_t> res;
    for (vertex_t v : tour)
        res.push_back(ids[v]);
    return res;
}

vector<vertex_t> to_mate(const graph_t &g, const vector<pvv> &edges) {
    vector<vertex_t> mate(num_vertices(g), g.null_vertex());
    for (pvv p : edges) {
        mate[p.first] = p.second;
        mate[p.second] = p.first;
    }
    return mate;
}

//  returns a sequence of vertices starting from s.
vector<vertex_t> eulerian_path(size_t n, const vector<pvv> &edges, vertex_t s) {
    vector<vector<pvv>> g(n);
    for (pvv p : edges) {
        vertex_t u = p.first, v = p.second;
        g[u].push_back({v, g[v].size()});
        g[v].push_back({u, g[u].size() - 1});
    }

    vector<vertex_t> stk, res;
    stk.push_back(s);
    while (!stk.empty()) {
        size_t u = stk.back();
        if (g[u].empty()) {
            res.push_back(u);
            stk.pop_back();
        }
        else {
            pvv p = g[u].back();
            size_t v = p.first, i = p.second;
            if (v != u) {
                if (i + 1 != g[v].size()) {
                    pvv &q = g[v].back();
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
    reverse(res.begin(), res.end());
    return res;
}

void rotate_eul_tour(vector<vertex_t> &tour) {
    rotate(tour.begin(), tour.begin() + 1, tour.end());
    tour.back() = tour.front();
}

vector<vertex_t> filter_eul_tour(size_t n, const vector<vertex_t> &tour) {
    vector<bool> appeared(n, false);
    vector<vertex_t> res;
    for (vertex_t v : tour) {
        if (appeared[v])
            continue;
        appeared[v] = true;
        res.push_back(v);
    }
    vertex_t v = res.front();
    res.push_back(v);
    return res;
}


weight_t tsp_weight(const graph_t& g, const vector<vertex_t>& tour, bool path = false) {
    const size_t n = num_vertices(g);
    weight_t sum = 0;
    for (size_t i = 0; i < n - (path ? 1 : 0); ++i) {
        vertex_t u = tour[i], v = tour[(i + 1) % n];
        auto e = boost::edge(u, v, g);
        if (!e.second) continue;
        sum += boost::get(edge_weight, g, e.first);
    }
    return sum;
}

bool is_perm(size_t n, const vector<vertex_t> &perm) {
    if (perm.size() != n)
        return false;
    vector<bool> appeared(n, false);
    for (vertex_t v : perm) {
        if (v >= n)
            return false;
        appeared[v] = true;
    }
    for (bool b : appeared)
        if (!b)
            return false;
    return true;
}

bool is_perfect(size_t n, const vector<vertex_t> &mate) {
    assert(n % 2 == 0);
    if (!is_perm(n, mate))
        return false;
    for (size_t u = 0; u < n; ++u) {
        if (mate[u] == u)
            return false;
        if (mate[mate[u]] != u)
            return false;
    }
    return true;
}

bool is_spp(size_t n, const vector<vertex_t>& tour, const vector<vertex_t>& mate) {
    assert(n % 2 == 0);
    if (!is_perm(n, tour))
        return false;
    for (size_t i = 0; i < n; i += 2) {
        vertex_t u = tour[i], v = tour[i + 1];
        if (mate[u] != v || mate[v] != u)
            return false;
    }
    return true;
}

//  tour must be restored
vector<vertex_t> splice_spp(const map<pvv, vector<vertex_t>>& hs, const vector<vertex_t>& tour) {
    const size_t n = tour.size();
    assert(n % 2 == 0);
    vector<vertex_t> res;
    for (size_t i = 0; i < n; i += 2) {
        pvv p(tour[i], tour[i + 1]);
        normalize(p);
        auto it = hs.find(p);
        assert(it != hs.end());
        auto subtour = it->second;
        if (subtour.front() != tour[i])
            reverse(subtour.begin(), subtour.end());
        assert(subtour.front() == tour[i] && subtour.back() == tour[i + 1]);
        for (vertex_t v : subtour)
            res.push_back(v);
    }
    return res;
}

using std::chrono::system_clock;
using std::chrono::seconds;
using std::chrono::milliseconds;
using std::chrono::duration_cast;

struct timer {
    system_clock::time_point t_start;
    auto now() { return system_clock::now(); }
    auto dur() { return duration_cast<milliseconds>(now() - t_start); }
    auto dur_s() { return duration_cast<milliseconds>(now() - t_start).count() / 1000.0; }
    void reset() { t_start = now(); }
    timer() : t_start(now()) {}
    template<class Duration>
    bool check(Duration m) { return dur() < m; }
};