#ifndef __HYBRID_SAMPLE_H__
#define __HYBRID_SAMPLE_H__

#include "AdoptionFunction.h"

enum SampleKind { RES = 0, CAP = 1, RAW = 2 };

class RRSetStore {
public:
    vector<u64> off{0};
    vector<u32> nodes;

    size_t count() const { return off.size() - 1; }

    void clear() {
        off.assign(1, 0);
        nodes.clear();
        nodes.shrink_to_fit();
    }
};

class RRSeqStore {
public:
    vector<u32> root, d, cc, nsrc;
    vector<u8> ss, selfsrc;
    vector<double> base, lin, lself, den;
    vector<u64> ent_off{0};
    vector<u32> ent_node, ent_seq;
    vector<u64> lab_off{0};
    vector<u32> lab;

    size_t count() const { return root.size(); }

    void clear() {
        root.clear(); d.clear(); cc.clear(); nsrc.clear(); ss.clear(); selfsrc.clear();
        base.clear(); lin.clear(); lself.clear(); den.clear();
        ent_off.assign(1, 0);
        ent_node.clear(); ent_seq.clear();
        lab_off.assign(1, 0);
        lab.clear();
        ent_node.shrink_to_fit(); ent_seq.shrink_to_fit(); lab.shrink_to_fit(); lab_off.shrink_to_fit();
    }
};

class SamplePool {
public:
    RRSetStore P;
    RRSeqStore Q;

    size_t total() const { return P.count() + Q.count(); }

    void clear() {
        P.clear();
        Q.clear();
    }
};

inline double seqValue(const AdoptionFunction& f, SampleKind kind, const RRSeqStore& Q, size_t j, u32 c, int s) {
    double fp = bv(f, Q.root[j]) * f.phi(Q.d[j], Q.cc[j] + c, Q.ss[j] | s) - Q.base[j];
    if (kind == RES) fp -= c * Q.lin[j] + s * Q.lself[j];
    double x = fp / Q.den[j];
    return x < 0 ? 0 : x;
}

class State {
public:
    vector<u8> active;
    vector<u32> c0;

    void init(u64 n) {
        active.assign(n, 0);
        c0.assign(n, 0);
    }
};

#endif
