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
    // print_graph(h);
    // assert(false);
    vector<vertex_t> mate = minimum_weighted_matching(h);

    auto matching_edges = restore(collect_matching_edges(mate), ids);

    vector<pvv> edges(move(tree_edges));
    edges.insert(edges.end(), matching_edges.begin(), matching_edges.end());

    vector<vertex_t> eul_tour = eulerian_path(n, edges, 0);

    pair<weight_t, vector<vertex_t>> res(inf, vector<vertex_t>());
    for (size_t i = 0; i < eul_tour.size(); ++i) {
        vector<vertex_t> tour;
        vector<bool> appeared(n, false);
        for (size_t j = 0; j < eul_tour.size(); ++j) {
            vertex_t v = eul_tour[(i + j) % eul_tour.size()];
            if (appeared[v]) continue;
            appeared[v] = true;
            tour.push_back(v);
        }
        res = min(res, make_pair(tsp_weight(g, tour, false), tour));
    }
    return res.second;
}

vector<vertex_t> spp(const graph_t& g, const vector<vertex_t>& mate) {
    const auto n = boost::num_vertices(g);
    vector<vertex_t> pre(n);
    auto new_weight = boost::make_function_property_map<graph_t::edge_descriptor, weight_t>([&](const graph_t::edge_descriptor& e) {
        auto src = boost::source(e, g);
        auto dst = boost::target(e, g);
        auto w = boost::get(edge_weight, g, e);
        if (mate[src] == dst && mate[dst] == src)
            return weight_t(0);
        return w;
    });
    prim_minimum_spanning_tree(g, pre.data(), boost::weight_map(new_weight));

    //  check whether every matching edge is inside the MST
    for (const vertex_t& v : mate)
        assert(pre[v] == mate[v] || pre[mate[v]] == v);

    auto tree_edges = collect_tree_edges(pre);


    auto deg = calc_deg(n, tree_edges);

    auto id = collect_odd_deg_vertices(deg);

    graph_t h = induce(g, id);
    // print_graph(h);
    // assert(false);

    vector<vertex_t> mate2 = minimum_weighted_matching(h);
    auto matching_edges = restore(collect_matching_edges(mate2), id);

    vector<pvv> edges(move(tree_edges));
    edges.insert(edges.end(), matching_edges.begin(), matching_edges.end());

    vector<vertex_t> eul_tour = eulerian_path(n, edges, 0);
    assert(eul_tour.size() == edges.size() + 1);
    
    pair<weight_t, vector<vertex_t>> res(inf, vector<vertex_t>());
    for (size_t i = 0; i < eul_tour.size(); ++i) {
        vector<vertex_t> tour;
        vector<bool> appeared(n, false);
        for (size_t j = 0; j < eul_tour.size(); ++j) {
            size_t u = eul_tour[(i + j) % eul_tour.size()], v = eul_tour[(i + j + 1) % eul_tour.size()];
            if (appeared[u]) continue;
            if (mate[u] != v) continue;
            tour.push_back(u);
            tour.push_back(v);
            appeared[u] = appeared[v] = true;
        }
        res = min(res, make_pair(tsp_weight(g, tour, false), tour));
    }
    return res.second;
}

vector<vertex_t> spp_induce(const graph_t& g, const map<pvv, vector<vertex_t>>& subtours) {
    vector<pvv> edges;
    for (const auto& p : subtours)
        edges.push_back(p.first);

    vmap gids;
    for (pvv p : edges) {
        gids.push_back(p.first);
        gids.push_back(p.second);
    }
    gids.build();
    for (pvv& p : edges) {
        p.first = gids.id(p.first);
        p.second = gids.id(p.second);
    }

    graph_t h = induce(g, gids);
    assert(num_vertices(h) == gids.size());
    vector<vertex_t> mate = to_mate(h, edges);

    assert(is_perfect(num_vertices(h), mate));
    vector<vertex_t> tour = spp(h, mate);
    assert(is_perm(num_vertices(h), tour));
    tour = restore(tour, gids);
    return splice_spp(subtours, tour);
}

vector<vertex_t> eul_to_tsp(size_t n, const vector<pvv>& edges, vertex_t s, vertex_t t) {
    vector<vertex_t> eul_tour = eulerian_path(n, edges, s);
    vector<bool> appeared(n, false);
    appeared[t] = true;
    vector<vertex_t> res;
    for (vertex_t v : eul_tour) {
        if (appeared[v]) continue;
        appeared[v] = true;
        res.push_back(v);
    }
    res.push_back(t);
    assert(is_perm(n, res) && res.front() == s && res.back() == t);
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
    // deg = calc_deg(n, edges);
    return eul_to_tsp(n, edges, s, t);
}

vector<vertex_t> tspp(const graph_t& g) {
    const size_t n = boost::num_vertices(g);
    vector<vertex_t> pre(n);
    prim_minimum_spanning_tree(g, pre.data());

    auto tree_edges = collect_tree_edges(pre);

    auto deg = calc_deg(n, tree_edges);

    auto ids = collect_odd_deg_vertices(deg);
    const size_t m = ids.size();
    graph_t h(m + 2);
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = i + 1; j < m; ++j) {
            auto e = g.get_edge(ids[i], ids[j]);
            if (!e.first)
                continue;
            add_edge(i, j, e.second, h);
        }
        add_edge(i, m, 0, h);
        add_edge(i, m + 1, 0, h);
    }
    ids.push_back(n);
    ids.push_back(n + 1);
    
    vector<vertex_t> mate2 = minimum_weighted_matching(h);
    auto matching_edges = restore(collect_matching_edges(mate2), ids);
    vector<pvv> edges(move(tree_edges));
    edges.insert(edges.end(), matching_edges.begin(), matching_edges.end());
    auto res = eul_to_tsp(n + 2, edges, n, n + 1);
    res.pop_back();
    res.erase(res.begin());
    return res;
}

//  dp[S][t] is the weight of the optimal tsp-(s-t)-path of G[S]
vector<vector<weight_t>> tspp_dp(const graph_t& g, vertex_t s) {
    const vertex_t n = num_vertices(g);
    vector<vector<weight_t>> dp(1 << n, vector<weight_t>(n, inf));
    dp[1 << s][s] = 0;
    for (vertex_t S = (1 << s); S < (1 << n) - 1; ++S) {
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
    for (vertex_t t = 0; t < n; ++t)
        if (t != s)
            assert(dp[(1 << n) - 1][t] != inf);
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
    reverse(res.begin(), res.end());
    return res;
}