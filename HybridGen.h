#ifndef __HYBRID_GEN_H__
#define __HYBRID_GEN_H__

#include "HybridSample.h"

class HybridGen {
public:
    const Graph& G;
    const AdoptionFunction& f;
    const State& st;
    SampleKind kind;
    const vector<double>& lin;
    const vector<double>& lself;
    const vector<double>& den;

    bool literal = false;

    vector<u32> stamp, slot;
    u32 cur = 0;
    vector<vector<u32>> label;
    vector<u32> node_of;
    vector<u8> inq, memo_done;
    vector<u64> memo_b, memo_e;
    vector<u32> memo, src, queue, tmp;
    vector<u64> lb;
    vector<u32> vis, bfs;
    u32 vcur = 0;

    HybridGen(const Graph& G_, const AdoptionFunction& f_, const State& st_, SampleKind k,
              const vector<double>& lin_, const vector<double>& lself_, const vector<double>& den_)
        : G(G_), f(f_), st(st_), kind(k), lin(lin_), lself(lself_), den(den_) {
        stamp.assign(G.n, 0);
        slot.assign(G.n, 0);
    }

    inline void newStamp() {
        if (++cur == 0) {
            fill(stamp.begin(), stamp.end(), 0);
            cur = 1;
        }
    }

    inline void newVisit() {
        if (vis.empty()) vis.assign(G.n, 0);
        if (++vcur == 0) {
            fill(vis.begin(), vis.end(), 0);
            vcur = 1;
        }
    }

    void genRRSet(u32 root, Rng& rng, RRSetStore& P) {
        newStamp();
        size_t b0 = P.nodes.size();
        P.nodes.push_back(root);
        stamp[root] = cur;
        for (size_t h = b0; h < P.nodes.size(); h++) {
            u32 w = P.nodes[h];
            G.sampleLiveIn(w, rng, [&](u32 z) {
                if (st.active[z] || stamp[z] == cur) return;
                stamp[z] = cur;
                P.nodes.push_back(z);
            });
        }
        P.off.push_back(P.nodes.size());
    }

    inline u32 getSlot(u32 x, u32& slots) {
        if (stamp[x] != cur) {
            stamp[x] = cur;
            slot[x] = slots;
            if (slots >= label.size()) {
                label.emplace_back();
                node_of.push_back(0);
                inq.push_back(0);
                memo_done.push_back(0);
                memo_b.push_back(0);
                memo_e.push_back(0);
            }
            label[slots].clear();
            node_of[slots] = x;
            inq[slots] = 0;
            memo_done[slots] = 0;
            return slots++;
        }
        return slot[x];
    }

    void collectSources(u32 v, u8& selfsrc) {
        src.clear();
        selfsrc = !st.active[v];
        if (selfsrc) src.push_back(v);
        for (u64 i = G.in_off[v]; i < G.in_off[v + 1]; i++) {
            u32 u = G.in_src[i];
            if (!st.active[u]) src.push_back(u);
        }
    }

    void genRRSeq(u32 v, Rng& rng, RRSeqStore& Q) {
        if (literal) genRRSeqLiteral(v, rng, Q);
        else genRRSeqLabel(v, rng, Q);
    }

    void genRRSeqLabel(u32 v, Rng& rng, RRSeqStore& Q) {
        newStamp();
        u32 slots = 0;
        memo.clear();
        queue.clear();
        u8 selfsrc;
        collectSources(v, selfsrc);
        if (src.size() <= 1024) {
            genRRSeqBits(v, rng, Q, selfsrc, slots);
            return;
        }
        for (u32 t = 0; t < src.size(); t++) {
            u32 sl = getSlot(src[t], slots);
            label[sl].push_back(t);
            if (!inq[sl]) {
                inq[sl] = 1;
                queue.push_back(sl);
            }
        }
        size_t head = 0;
        while (head < queue.size()) {
            u32 sl = queue[head++];
            inq[sl] = 0;
            expand(sl, rng);
            for (u64 i = memo_b[sl]; i < memo_e[sl]; i++) {
                u32 zs = getSlot(memo[i], slots);
                const vector<u32>& Lw = label[sl];
                vector<u32>& Lz = label[zs];
                if (!includes(Lz.begin(), Lz.end(), Lw.begin(), Lw.end())) {
                    tmp.clear();
                    set_union(Lz.begin(), Lz.end(), Lw.begin(), Lw.end(), back_inserter(tmp));
                    Lz.swap(tmp);
                    if (!inq[zs]) {
                        inq[zs] = 1;
                        queue.push_back(zs);
                    }
                }
            }
            compactQueue(head);
        }
        storeRRSeq(v, selfsrc, slots, Q);
    }

    void genRRSeqBits(u32 v, Rng& rng, RRSeqStore& Q, u8 selfsrc, u32& slots) {
        const u32 W = ((u32)src.size() + 63) / 64;
        auto bslot = [&](u32 x) -> u32 {
            u32 before = slots, sl = getSlot(x, slots);
            if (slots != before) {
                if (lb.size() < (size_t)slots * W) lb.resize((size_t)slots * W * 2);
                fill(lb.begin() + (size_t)sl * W, lb.begin() + (size_t)(sl + 1) * W, 0ULL);
            }
            return sl;
        };
        for (u32 t = 0; t < src.size(); t++) {
            u32 sl = bslot(src[t]);
            lb[(size_t)sl * W + (t >> 6)] |= 1ULL << (t & 63);
            if (!inq[sl]) {
                inq[sl] = 1;
                queue.push_back(sl);
            }
        }
        size_t head = 0;
        while (head < queue.size()) {
            u32 sl = queue[head++];
            inq[sl] = 0;
            expand(sl, rng);
            for (u64 i = memo_b[sl]; i < memo_e[sl]; i++) {
                u32 zs = bslot(memo[i]);
                const u64* Lw = lb.data() + (size_t)sl * W;
                u64* Lz = lb.data() + (size_t)zs * W;
                bool sub = true;
                for (u32 w = 0; w < W; w++) {
                    if (Lw[w] & ~Lz[w]) { sub = false; break; }
                }
                if (!sub) {
                    for (u32 w = 0; w < W; w++) Lz[w] |= Lw[w];
                    if (!inq[zs]) {
                        inq[zs] = 1;
                        queue.push_back(zs);
                    }
                }
            }
            compactQueue(head);
        }
        for (u32 sl = 0; sl < slots; sl++) {
            label[sl].clear();
            const u64* L = lb.data() + (size_t)sl * W;
            for (u32 w = 0; w < W; w++) {
                for (u64 x = L[w]; x; x &= x - 1) label[sl].push_back((w << 6) + (u32)__builtin_ctzll(x));
            }
        }
        storeRRSeq(v, selfsrc, slots, Q);
    }

    void genRRSeqLiteral(u32 v, Rng& rng, RRSeqStore& Q) {
        newStamp();
        u32 slots = 0;
        memo.clear();
        u8 selfsrc;
        collectSources(v, selfsrc);
        for (u32 t = 0; t < src.size(); t++) {
            newVisit();
            bfs.clear();
            bfs.push_back(src[t]);
            vis[src[t]] = vcur;
            for (size_t h = 0; h < bfs.size(); h++) {
                u32 sl = getSlot(bfs[h], slots);
                label[sl].push_back(t);
                expand(sl, rng);
                for (u64 i = memo_b[sl]; i < memo_e[sl]; i++) {
                    u32 z = memo[i];
                    if (vis[z] != vcur) {
                        vis[z] = vcur;
                        bfs.push_back(z);
                    }
                }
            }
        }
        storeRRSeq(v, selfsrc, slots, Q);
    }

    inline void expand(u32 sl, Rng& rng) {
        if (memo_done[sl]) return;
        memo_done[sl] = 1;
        memo_b[sl] = memo.size();
        G.sampleLiveIn(node_of[sl], rng, [&](u32 z) {
            if (!st.active[z]) memo.push_back(z);
        });
        memo_e[sl] = memo.size();
    }

    inline void compactQueue(size_t& head) {
        if (head > (1u << 20) && head * 2 > queue.size()) {
            queue.erase(queue.begin(), queue.begin() + head);
            head = 0;
        }
    }

    void storeRRSeq(u32 v, u8 selfsrc, u32 slots, RRSeqStore& Q) {
        size_t j = Q.root.size();
        u32 dv = G.indeg(v);
        Q.root.push_back(v);
        Q.d.push_back(dv);
        Q.cc.push_back(st.c0[v]);
        Q.ss.push_back(st.active[v]);
        Q.nsrc.push_back((u32)src.size());
        Q.selfsrc.push_back(selfsrc);
        Q.base.push_back(bv(f, v) * f.phi(dv, st.c0[v], st.active[v]));
        Q.lin.push_back(lin.empty() ? 0.0 : lin[v]);
        Q.lself.push_back(lself.empty() ? 0.0 : lself[v]);
        Q.den.push_back(kind == RAW ? 1.0 : den[v]);
        for (u32 sl = 0; sl < slots; sl++) {
            Q.ent_node.push_back(node_of[sl]);
            Q.ent_seq.push_back((u32)j);
            Q.lab.insert(Q.lab.end(), label[sl].begin(), label[sl].end());
            Q.lab_off.push_back(Q.lab.size());
        }
        Q.ent_off.push_back(Q.ent_node.size());
    }
};

#endif
