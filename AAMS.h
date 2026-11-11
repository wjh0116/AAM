#ifndef __AAMS_H__
#define __AAMS_H__

#include "MODE.h"

struct Realization {
    vector<u64> off;
    vector<u32> dst;
};

class AAMS {
public:
    const Graph& G;
    string algo;
    AdoptionFunction fe, fo;
    SampleKind kind;
    u32 k, b;
    u64 seed;
    State st;
    MODE* mode = nullptr;

    vector<BatchStat> stats;
    vector<vector<u32>> seeds;
    double adoption = 0, running_time = 0;
    u64 activated = 0, total_sets = 0, total_seqs = 0;

    AAMS(const Graph& G_, const string& algo_, const AdoptionFunction& fe_, const AdoptionFunction& fo_,
         u32 k_, u32 b_, double eps, u64 seed_)
        : G(G_), algo(algo_), fe(fe_), fo(fo_), k(k_), b(b_), seed(seed_) {
        if (algo == "AdaMODE-") kind = CAP;
        else if (algo == "HLA") kind = RAW;
        else kind = RES;
        if (algo == "MODE" || algo == "HLA") b = k;
        b = min(b, k);
        st.init(G.n);
        double t0 = nowSec();
        if (algo != "AdaDegree") mode = new MODE(G, fo, kind, st, eps);
        running_time = nowSec() - t0;
    }

    ~AAMS() { delete mode; }

    Realization makeRealization(u64 id) const {
        u64 n = G.n;
        Rng rr(mix64(0xC0FFEEULL, id, 7));
        Realization R;
        R.off.assign(n + 1, 0);
        vector<pair<u32, u32>> live;
        for (u64 v = 0; v < n; v++) {
            G.sampleLiveIn((u32)v, rr, [&](u32 z) { live.push_back({z, (u32)v}); });
        }
        for (auto& e : live) R.off[e.first + 1]++;
        for (u64 i = 0; i < n; i++) R.off[i + 1] += R.off[i];
        R.dst.resize(live.size());
        vector<u64> pos(R.off.begin(), R.off.end() - 1);
        for (auto& e : live) R.dst[pos[e.first]++] = e.second;
        return R;
    }

    size_t deploy(const Realization& R, const vector<u32>& S) {
        vector<u32> q;
        for (u32 x : S) {
            if (x < G.n && !st.active[x]) {
                st.active[x] = 1;
                q.push_back(x);
            }
        }
        for (size_t h = 0; h < q.size(); h++) {
            u32 u = q[h];
            for (u64 j = R.off[u]; j < R.off[u + 1]; j++) {
                u32 v = R.dst[j];
                if (!st.active[v]) {
                    st.active[v] = 1;
                    q.push_back(v);
                }
            }
        }
        for (u32 u : q) {
            for (u64 j = G.out_off[u]; j < G.out_off[u + 1]; j++) st.c0[G.out_dst[j]]++;
        }
        return q.size();
    }

    vector<u32> degreeSeeds(u32 bb, BatchStat& bs) {
        double t0 = nowSec();
        vector<u32> inactive;
        for (u64 v = 0; v < G.n; v++) if (!st.active[v]) inactive.push_back((u32)v);
        bs.n_psi = inactive.size();
        if (inactive.size() <= bb) {
            bs.time = nowSec() - t0;
            return inactive;
        }
        vector<pair<u64, u32>> sc;
        sc.reserve(inactive.size());
        for (u32 q : inactive) {
            u64 d = 0;
            for (u64 i = G.out_off[q]; i < G.out_off[q + 1]; i++) d += !st.active[G.out_dst[i]];
            sc.push_back({d, q});
        }
        partial_sort(sc.begin(), sc.begin() + bb, sc.end(), [](const pair<u64, u32>& x, const pair<u64, u32>& y) {
            return x.first > y.first || (x.first == y.first && x.second < y.second);
        });
        vector<u32> S;
        for (u32 i = 0; i < bb; i++) S.push_back(sc[i].second);
        bs.time = nowSec() - t0;
        return S;
    }

    void run(u64 id) {
        Realization R = makeRealization(id);
        double t0 = nowSec();
        u32 r = (k + b - 1) / b;
        for (u32 i = 0; i < r; i++) {
            BatchStat bs;
            u32 bb = min<u32>(b, k - i * b);
            vector<u32> S;
            if (algo == "AdaDegree") S = degreeSeeds(bb, bs);
            else if (algo == "HLA") S = mode->select(bb, i, seed, bs, true, args.theta);
            else S = mode->select(bb, i, seed, bs);
            bs.newly = deploy(R, S);
            total_sets += bs.nsets;
            total_seqs += bs.nseqs;
            cout << "Batch " << i << ": |S|=" << S.size() << " newly activated=" << bs.newly
                 << " |V_psi|=" << bs.n_psi << " beta=" << (bs.F > 0 ? bs.H / bs.F : 0)
                 << " RR sets=" << bs.nsets << " RR sequences=" << bs.nseqs
                 << " time=" << bs.time << "s" << endl;
            stats.push_back(bs);
            seeds.push_back(S);
        }
        running_time += nowSec() - t0;

        adoption = 0;
        activated = 0;
        for (u64 v = 0; v < G.n; v++) {
            u32 d = G.indeg((u32)v);
            u32 c = st.c0[v];
            int s = st.active[v];
            if (!c && !s) continue;
            adoption += bv(fe, (u32)v) * fe.phi(d, c, s);
            activated += s;
        }
    }

    void printResults() const {
        cout << "======================================" << endl;
        cout << "Algorithm: " << algo << endl;
        cout << "k: " << k << ", b: " << b << endl;
        cout << "Adoption spread: " << adoption << endl;
        cout << "Activated nodes: " << activated << endl;
        cout << "Samples: RR sets=" << total_sets << " RR sequences=" << total_seqs << endl;
        cout << "Running time: " << running_time << " seconds" << endl;
        cout << "======================================" << endl;
    }

    void saveResults(const string& path) const {
        ofstream fout(path);
        fout << "algorithm=" << algo << endl;
        fout << "k=" << k << endl;
        fout << "b=" << b << endl;
        fout << "adoption_spread=" << fixed << setprecision(6) << adoption << endl;
        fout << "activated_nodes=" << activated << endl;
        fout << "rr_sets=" << total_sets << endl;
        fout << "rr_sequences=" << total_seqs << endl;
        fout << "running_time=" << running_time << endl;
        fout << "seeds:" << endl;
        for (size_t i = 0; i < seeds.size(); i++) {
            for (size_t j = 0; j < seeds[i].size(); j++) fout << (j ? " " : "") << seeds[i][j];
            fout << endl;
        }
        fout.close();
    }
};

#endif
