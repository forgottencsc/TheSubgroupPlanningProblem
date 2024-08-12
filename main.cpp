#include <bits/stdc++.h>
#include "common.hpp"

#include "tsp.hpp"

using namespace std;

graph_t read_simiple_graph(istream& is) {
    size_t n, m;
    is >> n >> m;

    graph_t g(n);

    for (int i = 0; i < m; ++i) {
        vertex_t u, v;
        int64_t w;
        is >> u >> v >> w;
        u--;
        v--;
        boost::add_edge(u, v, edge_property_t(w), g);
    }
    return g;
}

graph_t read_vlsi(istream& is) {
    typedef pair<double, double> pdd;
    vector<pdd> ps;
    string line;
    while (getline(is, line)) {
        if (line.empty() || !isdigit(line[0])) continue;
        int id;
        pdd p;
        sscanf(line.c_str(), "%d %lf %lf", &id, &p.first, &p.second);
        ps.push_back(p);
    }
    const size_t n = ps.size();
    graph_t g(n);
    for (size_t i = 0; i < n; ++i)
        for (size_t j = i + 1; j < n; ++j)
            boost::add_edge(i, j, hypot(ps[i].first - ps[j].first, ps[i].second - ps[j].second), g);
    return g;
}

int main(void) {
    // const graph_t& g = read_vlsi(ifs);
    // const auto v = tspp(g, 0, 1);
    // fmt::print("{}", v);
    ifstream ifs("/root/tspp/vlsi/xql662.tsp");
    const graph_t& g = read_vlsi(ifs);
    auto mate = minimum_weighted_matching(g);
    auto tour = spp(g, mate);
    fmt::print("{}", is_spp(num_vertices(g), tour, mate));
    return 0;
}