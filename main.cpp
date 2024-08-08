#include <bits/stdc++.h>
#include "common.hpp"

#include "tsp.hpp"

using namespace std;

int main(void) {
    size_t n, m;
    cin >> n >> m;

    graph_t g(n);

    for (int i = 0; i < m; ++i) {
        vertex_t u, v;
        int64_t w;
        cin >> u >> v >> w;
        u--;
        v--;
        boost::add_edge(u, v, edge_property_t(w), g);
    }

    auto res = tsp(g);


    fmt::print("{}\n", res);
    

}