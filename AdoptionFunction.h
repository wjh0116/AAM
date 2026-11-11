#ifndef __ADOPTION_FUNCTION_H__
#define __ADOPTION_FUNCTION_H__

#include "graph.h"

enum FuncType { EQ1, THRESH, THRESH_TRUNC, THRESH_MOD, INFL };

class AdoptionFunction {
public:
    FuncType t = EQ1;
    double gamma = 0.5, lng = log(0.5);
    int T = 5;

    inline double phi(u32 d, u32 c, int s) const {
        switch (t) {
            case EQ1: return d ? (1.0 - exp(lng * (double)c / d) + gamma * s) : gamma * s;
            case THRESH: return (c + (u32)s >= (u32)T) ? 1.0 : 0.0;
            case THRESH_TRUNC: return (d + 1 >= (u32)T) ? min<double>(c + s, T) / T : 0.0;
            case THRESH_MOD: return (d + 1 >= (u32)T) ? (double)(c + s) / T : 0.0;
            case INFL: return s;
        }
        return 0;
    }

    inline double lin(u32 d) const { return d ? phi(d, d, 1) - phi(d, d - 1, 1) : 0.0; }

    inline double lself(u32 d) const { return phi(d, d, 1) - phi(d, d, 0); }

    static AdoptionFunction parse(const string& name, double gamma, int T) {
        AdoptionFunction f;
        f.gamma = gamma;
        f.lng = log(gamma);
        f.T = T;
        if (name == "eq1") f.t = EQ1;
        else if (name == "thresh") f.t = THRESH;
        else if (name == "thresh_trunc") f.t = THRESH_TRUNC;
        else if (name == "thresh_mod") f.t = THRESH_MOD;
        else if (name == "infl") f.t = INFL;
        else {
            cerr << "Unknown adoption function: " << name << endl;
            exit(1);
        }
        return f;
    }
};

vector<float> interest;

inline double bv(const AdoptionFunction& f, u32 v) {
    return (interest.empty() || f.t == INFL) ? 1.0 : (double)interest[v];
}

void generateInterest(const Graph& G, double frac, double low, u64 seed) {
    u64 n = G.n;
    interest.assign(n, (float)low);
    u64 target = (u64)ceil(frac * n), got = 0;
    vector<u8> seen(n, 0);
    vector<u32> q;
    Rng r(mix64(seed, 0xBA11));
    while (got < target) {
        u32 st0 = (u32)r.randomBelow(n);
        if (seen[st0]) continue;
        seen[st0] = 1;
        q.assign(1, st0);
        for (size_t h = 0; h < q.size() && got < target; h++) {
            u32 u = q[h];
            interest[u] = 1.0f;
            got++;
            auto visit = [&](u32 z) {
                if (!seen[z]) {
                    seen[z] = 1;
                    q.push_back(z);
                }
            };
            for (u64 i = G.out_off[u]; i < G.out_off[u + 1]; i++) visit(G.out_dst[i]);
            for (u64 i = G.in_off[u]; i < G.in_off[u + 1]; i++) visit(G.in_src[i]);
        }
    }
}

#endif
