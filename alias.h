#ifndef __ALIAS_H__
#define __ALIAS_H__

#include "head.h"

class Alias {
public:
    vector<u32> item, alias;
    vector<double> prob;
    double total = 0;

    void build(const vector<u32>& ids, const vector<double>& w) {
        size_t k = ids.size();
        item = ids;
        alias.assign(k, 0);
        prob.assign(k, 0);
        total = 0;
        for (double x : w) total += x;
        if (k == 0 || total <= 0) {
            item.clear();
            return;
        }
        vector<double> q(k);
        vector<u32> small, large;
        for (size_t i = 0; i < k; i++) {
            q[i] = w[i] * k / total;
            (q[i] < 1.0 ? small : large).push_back((u32)i);
        }
        while (!small.empty() && !large.empty()) {
            u32 s = small.back();
            small.pop_back();
            u32 l = large.back();
            prob[s] = q[s];
            alias[s] = l;
            q[l] = (q[l] + q[s]) - 1.0;
            if (q[l] < 1.0) {
                large.pop_back();
                small.push_back(l);
            }
        }
        for (u32 l : large) prob[l] = 1.0;
        for (u32 s : small) prob[s] = 1.0;
    }

    inline u32 draw(Rng& r) const {
        size_t i = r.randomBelow(item.size());
        return r.randomDouble() < prob[i] ? item[i] : item[alias[i]];
    }

    bool empty() const { return item.empty(); }
};

#endif
