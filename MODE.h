#ifndef __MODE_H__
#define __MODE_H__

#include "alias.h"
#include "HybridGen.h"
#include "DecGreedy.h"

struct BatchStat {
    double time = 0, W = 0, H = 0, F = 0, LB = 0, sl = 0, su = 0;
    u64 theta = 0, theta_max = 0, nsets = 0, nseqs = 0, n_psi = 0;
    int iters = 0;
    size_t newly = 0;
};

class MODE {
public:
    const Graph& G;
    AdoptionFunction fo;
    SampleKind kind;
    const State& st;
    double eps;

    vector<double> lin, lself, w;
    vector<double> eta, cap;
    vector<u32> topout;

    MODE(const Graph& G_, const AdoptionFunction& fo_, SampleKind kind_, const State& st_, double eps_)
        : G(G_), fo(fo_), kind(kind_), st(st_), eps(eps_) {
        u64 n = G.n;
        lin.resize(n);
        lself.resize(n);
        for (u64 v = 0; v < n; v++) {
            u32 d = G.indeg((u32)v);
            lin[v] = bv(fo, (u32)v) * fo.lin(d);
            lself[v] = bv(fo, (u32)v) * fo.lself(d);
        }
        w.assign(n, 0);
        for (u64 q = 0; q < n; q++) {
            double s = lself[q];
            for (u64 i = G.out_off[q]; i < G.out_off[q + 1]; i++) s += lin[G.out_dst[i]];
            w[q] = s;
        }
        topout.resize(n);
        iota(topout.begin(), topout.end(), 0);
        u64 keep = min<u64>(n, 2000);
        partial_sort(topout.begin(), topout.begin() + keep, topout.end(),
                     [&](u32 x, u32 y) { return G.outdeg(x) > G.outdeg(y); });
        topout.resize(keep);
        eta.assign(n, 0);
        cap.assign(n, 0);
    }

    double oneHopGain(u32 q) const {
        double g = bv(fo, q) * (fo.phi(G.indeg(q), st.c0[q], 1) - fo.phi(G.indeg(q), st.c0[q], 0));
        for (u64 i = G.out_off[q]; i < G.out_off[q + 1]; i++) {
            u32 v = G.out_dst[i];
            u32 d = G.indeg(v);
            g += bv(fo, v) * (fo.phi(d, st.c0[v] + 1, st.active[v]) - fo.phi(d, st.c0[v], st.active[v]));
        }
        return g;
    }

    void prepare(BatchStat& bs, vector<u32>& inactive, Alias& aW, Alias& aH) {
        u64 n = G.n;
        inactive.clear();
        vector<u32> idsW, idsH;
        vector<double> wW, wH;
        double W = 0, H = 0;
        for (u64 v = 0; v < n; v++) {
            u32 d = G.indeg((u32)v);
            u32 c0 = st.c0[v];
            int s0 = st.active[v];
            double b = bv(fo, (u32)v);
            double base = b * fo.phi(d, c0, s0);
            double cp = b * fo.phi(d, d, 1) - base;
            if (cp < 1e-15) cp = 0;
            double et = cp - (double)(d - c0) * lin[v] - (s0 ? 0.0 : lself[v]);
            if (et < 1e-12 * max(1.0, cp)) et = 0;
            cap[v] = cp;
            eta[v] = et;
            if (!s0) inactive.push_back((u32)v);
            if (kind == RES) {
                if (!s0 && w[v] > 0) { idsW.push_back((u32)v); wW.push_back(w[v]); W += w[v]; }
                if (et > 0) { idsH.push_back((u32)v); wH.push_back(et); H += et; }
            } else if (kind == CAP) {
                if (cp > 0) { idsH.push_back((u32)v); wH.push_back(cp); H += cp; }
            }
        }
        if (kind == RES) {
            aW.build(idsW, wW);
            aH.build(idsH, wH);
            bs.W = W;
            bs.H = H;
            bs.F = W + H;
        } else if (kind == CAP) {
            aW.item.clear();
            aH.build(idsH, wH);
            bs.W = 0;
            bs.H = H;
            bs.F = H;
        } else {
            vector<double> one(n, 1.0);
            vector<u32> all(n);
            iota(all.begin(), all.end(), 0);
            aH.build(all, one);
            aW.item.clear();
            bs.F = (double)n;
            bs.H = n;
            bs.W = 0;
        }
        double LB = 0;
        for (u32 q : topout) if (!st.active[q]) LB = max(LB, oneHopGain(q));
        if (LB <= 0) {
            for (size_t i = 0; i < inactive.size() && i < 100000; i++) LB = max(LB, oneHopGain(inactive[i]));
        }
        bs.LB = LB;
        bs.n_psi = inactive.size();
    }

    void draw(SamplePool& R, size_t target, Rng& rng, HybridGen& gen, const BatchStat& bs,
              const Alias& aW, const Alias& aH) {
        while (R.total() < target) {
            if (kind == RES) {
                if (bs.F <= 0) break;
                if (!aW.empty() && (aH.empty() || rng.randomDouble() * bs.F < bs.W)) gen.genRRSet(aW.draw(rng), rng, R.P);
                else gen.genRRSeq(aH.draw(rng), rng, R.Q);
            } else {
                if (aH.empty()) break;
                gen.genRRSeq(aH.draw(rng), rng, R.Q);
            }
        }
    }

    vector<u32> select(u32 b, u64 batch_idx, u64 seed, BatchStat& bs, bool literal = false, u64 theta_fixed = 0) {
        double t0 = nowSec();
        vector<u32> inactive;
        Alias aW, aH;
        prepare(bs, inactive, aW, aH);
        if (inactive.size() <= b) {
            bs.time = nowSec() - t0;
            return inactive;
        }
        HybridGen gen(G, fo, st, kind, lin, lself, kind == CAP ? cap : eta);
        gen.literal = literal;
        DecGreedy greedy(G.n, fo, kind);
        Rng r1(mix64(seed, batch_idx, 11)), r2(mix64(seed, batch_idx, 22));
        SamplePool R1, R2;

        if (theta_fixed > 0) {
            draw(R1, theta_fixed, r1, gen, bs, aW, aH);
            GreedyResult g = greedy.run(R1, b, inactive);
            bs.theta = R1.total();
            bs.nseqs = R1.Q.count();
            bs.iters = 1;
            bs.time = nowSec() - t0;
            return g.S;
        }

        double nps = (double)inactive.size();
        double delta = min(eps / 2, 1.0 / nps);
        double eps1 = (eps - delta) / (1 - delta);
        double rho = 1 - pow(1 - 1.0 / b, b);
        double epsc = rho * eps1 / (1 + rho);
        double lnC = lgamma(nps + 1) - lgamma((double)b + 1) - lgamma(nps - b + 1);
        double F = bs.F, LB = max(bs.LB, 1e-12);
        double tmax_d = ceil((2 + 2 * epsc / 3) * F * (lnC + log(8 / delta)) / (epsc * epsc * LB));
        u64 tmax = (u64)min(tmax_d, 1e18);
        u64 theta = (u64)max(1.0, ceil(tmax_d * LB / F));
        int T = (int)ceil(log2(max(2.0, F / LB))) + 1;
        double a = log(4.0 * T / delta);
        bs.theta_max = tmax;

        vector<u32> S;
        while (true) {
            bs.iters++;
            draw(R1, theta, r1, gen, bs, aW, aH);
            draw(R2, theta, r2, gen, bs, aW, aH);
            GreedyResult g = greedy.run(R1, b, inactive);
            S = g.S;
            double G2 = greedy.evaluate(R2, S);
            double th = (double)theta;
            double su = pow(sqrt(g.GammaU + a / 2) + sqrt(a / 2), 2) * F / th;
            double sl = max(0.0, pow(sqrt(G2 + 2 * a / 9) - sqrt(a / 2), 2) - a / 18) * F / th;
            bs.sl = sl;
            bs.su = su;
            bs.theta = theta;
            if (sl / su >= rho * (1 - eps1)) break;
            if (theta >= tmax) break;
            theta = min<u64>(2 * theta, tmax);
        }
        bs.nsets = R1.P.count() + R2.P.count();
        bs.nseqs = R1.Q.count() + R2.Q.count();
        bs.time = nowSec() - t0;
        return S;
    }
};

#endif
