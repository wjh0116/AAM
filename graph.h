#ifndef __GRAPH_H__
#define __GRAPH_H__

#include "head.h"

class Graph {
public:
    u64 n = 0, m = 0;

    vector<u64> in_off, out_off;
    vector<u32> in_src, out_dst;

    vector<float> in_p;
    vector<double> node_p;
    vector<double> logq;

    string folder;

    Graph(string folder) : folder(folder) {
        readAttribute();
        readGraph();
    }

    inline u32 indeg(u32 v) const { return (u32)(in_off[v + 1] - in_off[v]); }

    inline u32 outdeg(u32 u) const { return (u32)(out_off[u + 1] - out_off[u]); }

    void readAttribute() {
        string attr_file = folder + "/attribute.txt";
        ifstream fin(attr_file);
        if (!fin.is_open()) {
            cerr << "Error: cannot open " << attr_file << endl;
            exit(1);
        }
        string line;
        while (getline(fin, line)) {
            if (line.substr(0, 2) == "n=") {
                n = stoull(line.substr(2));
            } else if (line.substr(0, 2) == "m=") {
                m = stoull(line.substr(2));
            }
        }
        fin.close();
        cout << "Graph: n=" << n << ", m=" << m << endl;
    }

    void readGraph() {
        string gfile = folder + "/graph.txt";
        ifstream fin(gfile);
        if (!fin.is_open()) {
            cerr << "Error: cannot open " << gfile << endl;
            exit(1);
        }

        struct Arc { u32 u, v; float p; };
        vector<Arc> arcs;
        arcs.reserve(m);

        string first_line;
        getline(fin, first_line);
        istringstream iss(first_line);
        u64 u, v;
        float p;
        if (iss >> u >> v >> p) arcs.push_back({(u32)u, (u32)v, p});
        while (fin >> u >> v >> p) arcs.push_back({(u32)u, (u32)v, p});
        fin.close();

        cout << "Loaded " << arcs.size() << " edges from " << gfile << endl;
        assert(arcs.size() == m);

        sort(arcs.begin(), arcs.end(), [](const Arc& x, const Arc& y) {
            return x.u < y.u || (x.u == y.u && x.v < y.v);
        });

        in_off.assign(n + 1, 0);
        out_off.assign(n + 1, 0);
        for (auto& e : arcs) {
            out_off[e.u + 1]++;
            in_off[e.v + 1]++;
        }
        for (u64 i = 0; i < n; i++) {
            in_off[i + 1] += in_off[i];
            out_off[i + 1] += out_off[i];
        }

        in_src.resize(m);
        out_dst.resize(m);
        vector<float> prob_in(m);
        vector<u64> pi(in_off.begin(), in_off.end() - 1), po(out_off.begin(), out_off.end() - 1);
        for (auto& e : arcs) {
            out_dst[po[e.u]++] = e.v;
            prob_in[pi[e.v]] = e.p;
            in_src[pi[e.v]++] = e.u;
        }
        vector<Arc>().swap(arcs);

        bool uniform = true;
        for (u64 x = 0; x < n && uniform; x++) {
            for (u64 i = in_off[x] + 1; i < in_off[x + 1]; i++) {
                if (prob_in[i] != prob_in[in_off[x]]) { uniform = false; break; }
            }
        }

        node_p.assign(n, 0.0);
        logq.assign(n, 0.0);
        if (uniform) {
            for (u64 x = 0; x < n; x++) {
                u32 d = indeg((u32)x);
                if (d == 0) continue;
                double q = prob_in[in_off[x]];
                if (fabs(q - 1.0 / d) <= 1e-6) q = 1.0 / d;
                node_p[x] = q;
                logq[x] = (q >= 1.0) ? -INFINITY : log1p(-q);
            }
            cout << "Probabilities: shared by the in-edges of each node" << endl;
        } else {
            in_p.swap(prob_in);
            cout << "Probabilities: per edge" << endl;
        }
    }

    template <class F>
    inline void sampleLiveIn(u32 v, Rng& rng, F&& f) const {
        u64 b = in_off[v], e = in_off[v + 1];
        i64 d = (i64)(e - b);
        if (d == 0) return;
        if (!in_p.empty()) {
            for (u64 i = b; i < e; i++) {
                if (rng.randomDouble() < in_p[i]) f(in_src[i]);
            }
            return;
        }
        if (node_p[v] >= 1.0) {
            for (u64 i = b; i < e; i++) f(in_src[i]);
            return;
        }
        double lq = logq[v];
        i64 i = -1;
        while (true) {
            double r = rng.randomDoubleOpen0();
            double sk = floor(log(r) / lq);
            if (sk >= (double)d) break;
            i += (i64)sk + 1;
            if (i >= d) break;
            f(in_src[b + i]);
        }
    }
};

#endif
