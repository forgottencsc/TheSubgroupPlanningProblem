#pragma once 
#include "tsp.hpp"

vector<vertex_t> alg1(const graph_t& g, const vector<vector<vertex_t>>& c) {
    map<pvv, vector<vertex_t>> hs;
    vector<pvv> ps;
    vmap gids;
    for (const vector<vertex_t>& hids : c) {
        graph_t h = induce(g, hids);
        auto[w, p] = maximum_edge_weight(h);
        auto tour = tspp(h, p.first, p.second);

        p.first = hids[p.first];
        p.second = hids[p.second];
        hs[p] = tour;
        swap(p.first, p.second);
        hs[p] = tour;
        gids.push_back(p.first);
        gids.push_back(p.second);
        ps.push_back(p);
    }
    gids.build();
    ps = restore(ps, gids);
    auto mate = to_mate(g, ps);
    graph_t h = induce(g, gids);
    vector<vertex_t> tour = spp(h, mate);
    
}

vector<vertex_t> alg2(const graph_t& g, const vector<vector<vertex_t>>& c) {
    const int m = c.size();
    vector<graph_t> hs;
    for (const vector<vertex_t>& hids: c) 
        hs.emplace_back(induce(g, hids));

    vector<pair<weight_t, pvv>> mw;
    for (const graph_t& h : hs)
        mw.emplace_back(maximum_edge_weight(h));

    vector<vector<pvv>> he(m, vector<pvv>(m));

    graph_t h(m);
    for (size_t i = 0; i < m; ++i)
        for (size_t j = i + 1; j < m; ++j) {
            pair<weight_t, pvv> min_w(inf, pvv(-1, -1));
            for (size_t ii = 0; ii < c[i].size(); ++ii)
                for (size_t jj = 0; jj < c[j].size(); ++jj) {
                    auto p = boost::edge(c[i][ii], c[j][jj], g);
                    min_w = min(min_w, make_pair(get(edge_weight, g, p.first), pvv(ii, jj)));
                }
            add_edge(i, j, min_w.first + (mw[i].first + mw[j].first) / 2, h);
            he[i][j] = he[j][i] = min_w.second;
            swap(he[j][i].first, he[j][i].second);
        }
    

    vector<vertex_t> ans;

    vector<vertex_t> res = tsp(h);
    for (size_t i = 0; i < m; ++i) {
        const vertex_t u = res[i], v = res[(i + 1) % m], w = res[(i + 2) % m];
        const auto puv = he[u][v], pvw = he[v][w];
        vector<vertex_t> tour;
        if (puv.second != pvw.first)
            tour = tspp(hs[v], puv.second, pvw.first);
        else {
            tour = tsp(hs[v]);
            auto it = find(tour.begin(), tour.end(), puv.second);
            rotate(tour.begin(), it, tour.end());
        }
        tour = restore(tour, c[v]);
        for (vertex_t v : tour)
            ans.push_back(v);
    }

    return ans;
}

// vector<vertex_t> alg3(const graph_t& g, const vector<vector<vertex_t>>& c) {
//     vector<pair<weight_t, vector<vertex_t>>> res;
//     vector<pvv> ps;
//     vmap gids;
//     for (const vector<vertex_t>& hids : c) {
//         const graph_t& h = induce(g, hids);
//         const vertex_t m = num_vertices(h);
//         res.emplace_back(inf, vector<vertex_t>{});
//         for (vertex_t s = 0; s < m; ++s) {
//             const auto& dp = tspp_dp(h, s);
//             for (vertex_t t = s + 1; t < m; ++t) {
//                 auto e = boost::edge(s, t, h);
//                 if (!e.second) continue;
//                 auto w = dp[(1z << m) - 1][t] - boost::get(edge_weight, h, e.first);
//                 vector<vertex_t> tour = tspp_dp_route(h, dp, s, t);
//                 tour = restore(tour, hids);
//                 res.back() = min(res.back(), make_pair(w, tour));
//             }
//         }
//         const auto& tour = res.back().second;
//         pvv p(tour.front(), tour.back());
//         gids.push_back(p.first, p.second);
//     }



// }