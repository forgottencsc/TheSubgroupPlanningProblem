#pragma once 
#include "tsp.hpp"

vector<vertex_t> alg1(const graph_t& g, const vector<vector<vertex_t>>& c) {
    map<pvv, vector<vertex_t>> subtours;
    for (const vector<vertex_t>& hids : c) {
        graph_t h = induce(g, hids);
        if (num_vertices(h) == 1) {
            subtours[make_pair(hids[0], hids[0])] = vector<vertex_t>(hids[0]);
            continue;
        }
        auto[w, p] = maximum_edge_weight(h);
        auto tour = tspp(h, p.first, p.second);
        assert(tour.front() == p.first && tour.back() == p.second);
        tour = restore(tour, hids);
        p.first = hids[p.first];
        p.second = hids[p.second];
        normalize(p);
        subtours[p] = tour;
    }
    return spp_induce(g, subtours);
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
        if (puv.second != pvw.first) {
            tour = tspp(hs[v], puv.second, pvw.first);
            assert(tour.front() == puv.second);
            assert(tour.back() == pvw.first);
        }
        else {
            tour = tsp(hs[v]);
            auto it = find(tour.begin(), tour.end(), puv.second);
            rotate(tour.begin(), it, tour.end());
            assert(tour.front() == puv.second);
        }
        tour = restore(tour, c[v]);
        for (vertex_t v : tour)
            ans.push_back(v);
    }
    assert(is_perm(num_vertices(g), ans));
    return ans;
}

vector<vertex_t> alg3(const graph_t& g, const vector<vector<vertex_t>>& c) {
    map<pvv, vector<vertex_t>> subtours;
    for (const vector<vertex_t>& hids : c) {
        const graph_t& h = induce(g, hids);
        const vertex_t m = num_vertices(h);
        pair<weight_t, vector<vertex_t>> tmp(inf, vector<vertex_t>());
        for (vertex_t s = 0; s < m; ++s) {
            const auto& dp = tspp_dp(h, s);
            for (vertex_t t = s + 1; t < m; ++t) {
                auto e = boost::edge(s, t, h);
                if (!e.second) continue;
                auto w = dp[(1 << m) - 1][t] - boost::get(edge_weight, h, e.first);
                vector<vertex_t> tour = tspp_dp_route(h, dp, s, t);
                tmp = min(tmp, make_pair(w, tour));
            }
        }
        tmp.second = restore(tmp.second, hids);
        pvv p(tmp.second.front(), tmp.second.back());
        normalize(p);
        subtours[p] = tmp.second;
    }
    return spp_induce(g, subtours);
}

vector<vertex_t> alg4(const graph_t& g, const vector<vector<vertex_t>>& c) {
    const size_t n = num_vertices(g), k = c.size();
    vector<unordered_map<pvv, pair<weight_t, vector<vertex_t>>, hash_pvv>> dp(1 << k);

    for (size_t id = 0; id < k; ++id) {
        const graph_t& h = induce(g, c[id]);
        const size_t m = num_vertices(h);
        for (size_t i = 0; i < m; ++i)
            for (size_t j = i + 1; j < m; ++j) {
                auto tour = tspp(h, i, j);
                auto w = tsp_weight(h, tour, true);
                tour = restore(tour, c[id]);
                dp[1 << id][make_pair(c[id][i], c[id][j])] = make_pair(w, tour);
                reverse(tour.begin(), tour.end());
                dp[1 << id][make_pair(c[id][j], c[id][i])] = make_pair(w, tour);
            }
    }
    for (size_t S = 1; S < (1 << k); ++S) {
        for (const auto& p : dp[S]) {
            for (size_t x = 0; x < k; ++x) {
                if (S & (1 << x))
                    continue;
                for (const auto& q : dp[1 << x]) {
                    const auto[sp, tp] = p.first;
                    const auto[sq, tq] = q.first;
                    const auto lp = p.second;
                    const auto lq = q.second;
                    auto r = make_pair(sp, tq);
                    auto e = boost::edge(sq, tp, g);
                    auto w = boost::get(edge_weight, g, e.first);
                    vector<vertex_t> t(lp.second);
                    t.insert(t.end(), lq.second.begin(), lq.second.end());
                    pair<weight_t, vector<vertex_t>> l(lp.first+ w + lq.first, t);

                    auto ret = dp[S | (1 << x)].insert_or_assign(r, l);
                    if (ret.second)
                        continue;
                    auto it = ret.first;
                    it->second = min(it->second, l);
                }
            }
        }
    }
    pair<weight_t, vector<vertex_t>> ans(inf, vector<vertex_t>());
    for (size_t u = 0; u < n; ++u)
        for (size_t v = u + 1; v < n; ++v) {
            auto it = dp[(1 << k) - 1].find(make_pair(u, v));
            if (it != dp[(1 << k) - 1].end())
                ans = min(ans, it->second);
        }
    return ans.second;
}

vector<vertex_t> algc(const graph_t& g, const vector<vector<vertex_t>>& c) {
    map<pvv, vector<vertex_t>> subtours;
    for (const vector<vertex_t>& hids : c) {
        const size_t m = hids.size();
        graph_t h(m + 1);
        for (size_t i = 0; i < m; ++i) {
            add_edge(i, m, 0, h);
            for (size_t j = i + 1; j < m; ++j) {
                auto e = g.get_edge(hids[i], hids[j]);
                if (!e.first)
                    continue;
                add_edge(i, j, e.second, h);
            }
        }
        const auto& dp = tspp_dp(h, m);
        pair<weight_t, vertex_t> tmp(inf, m);
        for (vertex_t t = 0; t < m; ++t)
            tmp = min(tmp, make_pair(dp[(1 << (m + 1)) - 1][t], t));

        vector<vertex_t> tour = tspp_dp_route(h, dp, m, tmp.second);
        tour.erase(tour.begin());
        tour = restore(tour, hids);
        pvv p(tour.front(), tour.back());
        normalize(p);
        subtours[p] = tour;
    }
    return spp_induce(g, subtours);
}