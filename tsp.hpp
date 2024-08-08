#pragma once
#include "common.hpp"
#include <boost/graph/metric_tsp_approx.hpp>
vector<vertex_t> tsp(const graph_t& g) {
    const auto n = boost::num_vertices(g);
    vector<vertex_t> pre(n);
    prim_minimum_spanning_tree(g, pre.data());

    auto tree_edges = collect_tree_edges(pre);

    auto deg = calc_deg(n, tree_edges);

    auto id = collect_odd_deg_vertices(deg);

    graph_t h = induce(g, id);

    vector<vertex_t> mate = minimum_weighted_matching(h);

    auto matching_edges = restore(collect_matching_edges(mate), id);

    vector<pvv> edges(move(tree_edges));
    edges.insert(edges.end(), matching_edges.begin(), matching_edges.end());

    vector<vertex_t> eul_tour = eulerian_circuit(n, edges, 0);

    pair<weight_t, vector<vertex_t>> res(DBL_MAX, {});
    for (size_t i = 0; i < edges.size(); ++i) {
        vector<vertex_t> tsp_tour = filter_eul_tour(n, eul_tour);
        rotate_eul_tour(eul_tour);
        res = min(res, make_pair(tsp_tour_weight(g, tsp_tour), tsp_tour));
    }
    return res.second;
}