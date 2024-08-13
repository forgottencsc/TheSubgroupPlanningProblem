#include "tspp.hpp"

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

graph_t apply_scale(graph_t g, const vector<vector<vertex_t>>& c, double x) {
    const size_t n = num_vertices(g);
    graph_t h(n);
    vector<vertex_t> color(n, g.null_vertex());
    for (size_t i = 0; i < c.size(); ++i)
        for (vertex_t v : c[i])
            color[v] = i;

    for (auto [it, end] = edges(g); it != end; ++it) {
        if (!it->exists())
            continue;
        auto src = boost::source(*it, g);
        auto dst = boost::target(*it, g);
        auto w = boost::get(edge_weight, g, *it);
        if (color[src] != color[dst])
            w *= x;
        add_edge(src, dst, edge_property_t(w), h);
    }

    return h;
}

vector<vector<vertex_t>> gen_subgroup(const string& group_by, const graph_t& g, size_t k) {
    const size_t n = num_vertices(g);
    if (group_by == "modulo") {
        const size_t m = (n - 1) / k + 1;
        vector<vector<vertex_t>> c(m);
        for (size_t i = 0; i < n; ++i)
            c[i % m].push_back(i);
        return c;
    }
    else if (group_by == "random") {
        const size_t seed = 0xdeadbeef;
        mt19937_64 mt(seed);
        vector<vertex_t> perm(n);
        iota(perm.begin(), perm.end(), vertex_t(0));
        vector<vector<vertex_t>> c;
        for (vertex_t v : perm) {
            if (c.empty() || c.back().size() == k)
                c.emplace_back();
            c.back().push_back(v);
        }
        return c;
    }
    else if (group_by == "greedy") {
        unordered_set<vertex_t> s;
        for (vertex_t i = 0; i < n; ++i)
            s.insert(i);
        
        vector<vector<vertex_t>> c;
        vector<weight_t> dis;
        while (!s.empty()) {
            vertex_t v = g.null_vertex();
            if (c.empty() || c.back().size() == k) {
                c.emplace_back();
                dis.assign(n, inf);
                v = *s.begin();
            }
            else {
                for (vertex_t u : s)
                    if (v == g.null_vertex() || dis[u] < dis[v])
                        v = u;
            }
            for (vertex_t u : s)
                dis[u] = min(dis[u], boost::get(edge_weight, g, boost::edge(u, v, g).first));
            c.back().push_back(v);
            s.erase(v);
        }
        return c;
    }
    else {
        assert(false);  // not implemented
    }
}

//  tspp <instance xxx.tsp> <alg: 1|2|3|c> <scale: X> <group by: modulo|random|greedy(todo)> <group size: k>
int main(int argc, char** argv) {
    string path = argv[1];
    string alg = argv[2];
    double x = strtod(argv[3], nullptr);
    string group_by = argv[4];
    size_t k = strtoull(argv[5], nullptr, 10);

    ifstream ifs(argv[1]);
    const graph_t g = read_vlsi(ifs);
    const size_t n = num_vertices(g);
    auto c = gen_subgroup(group_by, n, k);
    const graph_t h = apply_scale(g, c, x);

    vector<vertex_t> tour;
    if (alg == "1") {
        tour = alg1(h, c);
    }
    else if (alg == "2") {
        tour = alg2(h, c);
    }
    else if (alg == "3") {
        tour = alg3(h, c);
    }
    else if (alg == "4") {
        tour = alg4(h, c);
    }
    else if (alg == "c") {
        tour = algc(h, c);
    }

    weight_t ans = tsp_weight(h, tour, false);
    printf("%f\n", ans);

    return 0;
}