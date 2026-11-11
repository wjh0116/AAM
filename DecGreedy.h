#ifndef __DEC_GREEDY_H__
#define __DEC_GREEDY_H__

#include "HybridSample.h"

struct GreedyResult {
    vector<u32> S;
    double GammaS = 0;
    double GammaU = INFINITY;
};

class DecGreedy {
public:
    u64 n;
    const AdoptionFunction& f;
    SampleKind kind;

    vector<double> mg;
    vector<u8> inS, touchedf, evmark;
    vector<u32> touched;
    vector<u32> cnt;
    vector<u64> offP, offQ, invP, invQ;
    vector<u8> coveredP;
    vector<u64> covw_off, covbits;
    vector<u32> curC;
    vector<u8> curS;
    vector<double> curVal, contrib, vals;

    DecGreedy(u64 n_, const AdoptionFunction& f_, SampleKind k) : n(n_), f(f_), kind(k) {
        mg.assign(n, 0);
        inS.assign(n, 0);
        touchedf.assign(n, 0);
        evmark.assign(n, 0);
    }

    inline void touch(u32 x) {
        if (!touchedf[x]) {
            touchedf[x] = 1;
            touched.push_back(x);
        }
    }

    double topB(u32 b) {
        vals.clear();
        for (u32 x : touched) if (!inS[x]) vals.push_back(mg[x]);
        double s = 0;
        if (vals.size() > b) {
            nth_element(vals.begin(), vals.begin() + b, vals.end(), greater<double>());
            for (u32 i = 0; i < b; i++) s += vals[i];
        } else {
            for (double x : vals) s += x;
        }
        return s;
    }

    inline bool isSelf(const RRSeqStore& Q, size_t j, u32 t) const { return Q.selfsrc[j] && t == 0; }

    inline void uncovered(const RRSeqStore& Q, size_t j, u64 e, u32& nc, int& ns) const {
        nc = 0;
        ns = 0;
        const u64* words = covbits.data() + covw_off[j];
        for (u64 i = Q.lab_off[e]; i < Q.lab_off[e + 1]; i++) {
            u32 t = Q.lab[i];
            if (!((words[t >> 6] >> (t & 63)) & 1)) {
                if (isSelf(Q, j, t)) ns = 1;
                else nc++;
            }
        }
    }

    GreedyResult run(const SamplePool& R, u32 b, const vector<u32>& fallback) {
        const RRSetStore& P = R.P;
        const RRSeqStore& Q = R.Q;
        for (u32 x : touched) {
            mg[x] = 0;
            inS[x] = 0;
            touchedf[x] = 0;
        }
        touched.clear();
        size_t nP = P.count(), nQ = Q.count();

        cnt.assign(n + 1, 0);
        for (u32 x : P.nodes) cnt[x + 1]++;
        offP.assign(n + 1, 0);
        for (u64 i = 0; i < n; i++) offP[i + 1] = offP[i] + cnt[i + 1];
        invP.resize(P.nodes.size());
        {
            vector<u64> pos(offP.begin(), offP.end() - 1);
            for (size_t j = 0; j < nP; j++) {
                for (u64 i = P.off[j]; i < P.off[j + 1]; i++) {
                    u32 x = P.nodes[i];
                    invP[pos[x]++] = j;
                    mg[x] += 1.0;
                    touch(x);
                }
            }
        }
        coveredP.assign(nP, 0);

        size_t nE = Q.ent_node.size();
        cnt.assign(n + 1, 0);
        for (u32 x : Q.ent_node) cnt[x + 1]++;
        offQ.assign(n + 1, 0);
        for (u64 i = 0; i < n; i++) offQ[i + 1] = offQ[i] + cnt[i + 1];
        invQ.resize(nE);
        {
            vector<u64> pos(offQ.begin(), offQ.end() - 1);
            for (u64 e = 0; e < nE; e++) invQ[pos[Q.ent_node[e]]++] = e;
        }
        covw_off.assign(nQ + 1, 0);
        for (size_t j = 0; j < nQ; j++) covw_off[j + 1] = covw_off[j] + (Q.nsrc[j] + 63) / 64;
        covbits.assign(covw_off[nQ], 0);
        curC.assign(nQ, 0);
        curS.assign(nQ, 0);
        curVal.assign(nQ, 0.0);
        contrib.assign(nE, 0.0);
        for (size_t j = 0; j < nQ; j++) {
            for (u64 e = Q.ent_off[j]; e < Q.ent_off[j + 1]; e++) {
                u32 nc;
                int ns;
                uncovered(Q, j, e, nc, ns);
                double c = seqValue(f, kind, Q, j, nc, ns);
                contrib[e] = c;
                u32 x = Q.ent_node[e];
                mg[x] += c;
                touch(x);
            }
        }

        GreedyResult out;
        for (u32 it = 0; it < b; it++) {
            out.GammaU = min(out.GammaU, out.GammaS + topB(b));
            u32 s = UINT32_MAX;
            double best = -1;
            for (u32 x : touched) {
                if (!inS[x] && mg[x] > best) {
                    best = mg[x];
                    s = x;
                }
            }
            if (s == UINT32_MAX || best <= 0) {
                for (u32 x : fallback) if (!inS[x] && !touchedf[x]) { s = x; break; }
                if (s == UINT32_MAX) for (u32 x : touched) if (!inS[x]) { s = x; break; }
                if (s == UINT32_MAX) break;
                if (!touchedf[s]) touch(s);
            }
            out.S.push_back(s);
            out.GammaS += max(0.0, mg[s]);
            inS[s] = 1;

            for (u64 i = offP[s]; i < offP[s + 1]; i++) {
                u64 j = invP[i];
                if (coveredP[j]) continue;
                coveredP[j] = 1;
                for (u64 k = P.off[j]; k < P.off[j + 1]; k++) mg[P.nodes[k]] -= 1.0;
            }
            for (u64 i = offQ[s]; i < offQ[s + 1]; i++) {
                u64 e = invQ[i];
                u32 j = Q.ent_seq[e];
                u64* words = covbits.data() + covw_off[j];
                bool changed = false;
                for (u64 q = Q.lab_off[e]; q < Q.lab_off[e + 1]; q++) {
                    u32 t = Q.lab[q];
                    if (!((words[t >> 6] >> (t & 63)) & 1)) {
                        words[t >> 6] |= (1ULL << (t & 63));
                        changed = true;
                        if (isSelf(Q, j, t)) curS[j] = 1;
                        else curC[j]++;
                    }
                }
                if (!changed) continue;
                curVal[j] = seqValue(f, kind, Q, j, curC[j], curS[j]);
                for (u64 e2 = Q.ent_off[j]; e2 < Q.ent_off[j + 1]; e2++) {
                    u32 nc;
                    int ns;
                    uncovered(Q, j, e2, nc, ns);
                    double c = (nc || ns) ? seqValue(f, kind, Q, j, curC[j] + nc, curS[j] | ns) - curVal[j] : 0.0;
                    if (c < 0) c = 0;
                    mg[Q.ent_node[e2]] += c - contrib[e2];
                    contrib[e2] = c;
                }
            }
        }
        out.GammaU = min(out.GammaU, out.GammaS + topB(b));
        return out;
    }

    double evaluate(const SamplePool& R, const vector<u32>& S) {
        const RRSetStore& P = R.P;
        const RRSeqStore& Q = R.Q;
        vector<u8>& mark = evmark;
        vector<u32> marked;
        for (u32 x : S) {
            if (!mark[x]) {
                mark[x] = 1;
                marked.push_back(x);
            }
        }
        double sum = 0;
        for (size_t j = 0; j < P.count(); j++) {
            for (u64 i = P.off[j]; i < P.off[j + 1]; i++) {
                if (mark[P.nodes[i]]) { sum += 1.0; break; }
            }
        }
        vector<u64> w;
        for (size_t j = 0; j < Q.count(); j++) {
            w.assign((Q.nsrc[j] + 63) / 64, 0);
            bool any = false;
            for (u64 e = Q.ent_off[j]; e < Q.ent_off[j + 1]; e++) {
                if (!mark[Q.ent_node[e]]) continue;
                any = true;
                for (u64 q = Q.lab_off[e]; q < Q.lab_off[e + 1]; q++) {
                    u32 t = Q.lab[q];
                    w[t >> 6] |= 1ULL << (t & 63);
                }
            }
            if (!any) continue;
            u32 c = 0;
            int s = 0;
            for (u32 t = 0; t < Q.nsrc[j]; t++) {
                if ((w[t >> 6] >> (t & 63)) & 1) {
                    if (Q.selfsrc[j] && t == 0) s = 1;
                    else c++;
                }
            }
            sum += seqValue(f, kind, Q, j, c, s);
        }
        for (u32 x : marked) mark[x] = 0;
        return sum;
    }
};

#endif
