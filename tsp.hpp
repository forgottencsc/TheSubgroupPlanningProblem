#pragma once
#include "matching.hpp"

vector<vertex_t> tsp(const graph_t& g) {
    const auto n = boost::num_vertices(g);
    vector<vertex_t> pre(n);
    prim_minimum_spanning_tree(g, pre.data());

    auto tree_edges = collect_tree_edges(pre);

    auto deg = calc_deg(n, tree_edges);

    auto ids = collect_odd_deg_vertices(deg);

    graph_t h = induce(g, ids);

    vector<vertex_t> mate = minimum_weighted_matching(h);

    auto matching_edges = restore(collect_matching_edges(mate), ids);

    vector<pvv> edges(move(tree_edges));
    edges.insert(edges.end(), matching_edges.begin(), matching_edges.end());

    vector<vertex_t> eul_tour = eulerian_path(n, edges, 0);

    vector<bool> appeared(n, false);
    vector<vertex_t> res;
    for (vertex_t v : eul_tour) {
        if (appeared[v]) continue;
        appeared[v] = true;
        res.push_back(v);
    }
    return res;
}

vector<vertex_t> spp(const graph_t& g, const vector<vertex_t>& mate) {
    const auto n = boost::num_vertices(g);
    vector<vertex_t> pre(n);
    auto new_weight = boost::make_function_property_map<graph_t::edge_descriptor, weight_t>([&](const graph_t::edge_descriptor& e) {
        auto src = boost::source(e, g);
        auto dst = boost::target(e, g);
        auto w = boost::get(boost::edge_weight, g, e);
        if (mate[src] == dst && mate[dst] == src)
            return weight_t(0);
        return w;
    });
    prim_minimum_spanning_tree(g, pre.data(), boost::weight_map(new_weight));


    auto tree_edges = collect_tree_edges(pre);

    auto deg = calc_deg(n, tree_edges);

    auto id = collect_odd_deg_vertices(deg);

    graph_t h = induce(g, id);

    vector<vertex_t> mate2 = minimum_weighted_matching(h);
    auto matching_edges = restore(collect_matching_edges(mate2), id);

    vector<pvv> edges(move(tree_edges));
    edges.insert(edges.end(), matching_edges.begin(), matching_edges.end());

    vector<vertex_t> eul_tour = eulerian_path(n, edges, 0);

    vector<bool> appeared(n, false);
    vector<vertex_t> res;
    for (size_t i = 0; i + 1 < eul_tour.size(); ++i) {
        size_t u = eul_tour[i], v = eul_tour[i + 1];
        if (appeared[u]) continue;
        if (mate[u] != v) continue;
        res.push_back(u);
        res.push_back(v);
        appeared[u] = appeared[v] = true;
    }

    return res;
}

vector<vertex_t> tspp(const graph_t& g, vertex_t s, vertex_t t) {
    const auto n = boost::num_vertices(g);
    vector<vertex_t> pre(n);
    prim_minimum_spanning_tree(g, pre.data());

    auto tree_edges = collect_tree_edges(pre);

    auto deg = calc_deg(n, tree_edges);

    deg[s] ^= 1;
    deg[t] ^= 1;

    const auto ids = collect_odd_deg_vertices(deg);
    
    graph_t h = induce(g, ids);

    vector<vertex_t> mate2 = minimum_weighted_matching(h);
    auto matching_edges = restore(collect_matching_edges(mate2), ids);

    vector<pvv> edges(move(tree_edges));
    edges.insert(edges.end(), matching_edges.begin(), matching_edges.end());
    deg = calc_deg(n, edges);
    
    vector<vertex_t> eul_tour = eulerian_path(n, edges, t);
    vector<bool> appeared(n, false);
    appeared[t] = true;
    vector<vertex_t> res;
    for (vertex_t v : eul_tour) {
        if (appeared[v]) continue;
        appeared[v] = true;
        res.push_back(v);
    }
    res.push_back(t);
    return res;
}

vector<vector<weight_t>> tspp_dp(const graph_t& g, vertex_t s) {
    const vertex_t n = num_vertices(g);
    vector<vector<weight_t>> dp(1 << n, vector<weight_t>(n, inf));
    dp[1 << s][s] = 0;
    for (vertex_t S = (1 << s) + 1; S < (1 << n) - 1; ++S) {
        for (vertex_t u = 0; u < n; ++u) {
            if (!(S & (1 << u)))
                continue;
            for (vertex_t v = 0; v < n; ++v) {
                if (S & (1 << v))
                    continue;
                auto e = boost::edge(u, v, g);
                if (!e.second) continue;
                auto w = boost::get(edge_weight, g, e.first);
                dp[S | (1 << v)][v] = min(dp[S | (1 << v)][v], dp[S][u] + w);
            }
        }
    }
    return dp;
}

vector<vertex_t> tspp_dp_route(const graph_t& g, const vector<vector<weight_t>>& dp, vertex_t s, vertex_t t) {
    const vertex_t n = num_vertices(g);
    vector<vertex_t> res;
    vertex_t S = (1 << n) - 1;
    while (t != s) {
        res.push_back(t);
        vertex_t T = S ^ (1 << t);
        for (vertex_t u = 0; u < n; ++u) {
            if (!(T & (1 << u)))
                continue;
            auto e = boost::edge(u, t, g);
            if (!e.second) continue;
            auto w = boost::get(edge_weight, g, e.first);
            if (dp[T][u] + w != dp[S][t])
                continue;
            S = T;
            t = u;
            break;
        }
    }
    res.push_back(t);
    return res;
}