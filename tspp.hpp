#pragma once 
#include "tsp.hpp"

vector<vertex_t> alg1(const graph_t& g, const vector<vector<vertex_t>>& c) {
    vector<vector<vertex_t>> hs;
    vector<pvv> ps;
    vmap gids;
    for (const vector<vertex_t>& hids : c) {
        graph_t h = induce(g, hids);
        auto[w, p] = maximum_edge_weight(h);
        gids.push_back(hids[p.first]);
        gids.push_back(hids[p.second]);
        ps.push_back(p);
        hs.push_back(tspp(h, p.first, p.second));
    }
    gids.build();
    graph_t h = induce(g, gids);
    vector<vertex_t> mate(gids.size());
    for (pvv& p : ps) {
        const vertex_t u = gids.id(p.first);
        const vertex_t v = gids.id(p.second);
        mate[u] = v;
        mate[v] = u;
    }
    return spp(h, mate);
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

vector<vertex_t> alg3(const graph_t& g, const vector<vector<vertex_t>>& c) {
    
}