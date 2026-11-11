#ifndef __HEAD_H__
#define __HEAD_H__

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <ctime>
#include <cassert>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <queue>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <functional>
#include <numeric>
#include <chrono>

using namespace std;

typedef uint8_t u8;
typedef uint32_t u32;
typedef uint64_t u64;
typedef int64_t i64;

inline double nowSec() {
    return chrono::duration<double>(chrono::steady_clock::now().time_since_epoch()).count();
}

class Rng {
public:
    u64 s[4];

    static u64 splitmix(u64& x) {
        u64 z = (x += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }

    explicit Rng(u64 seed = 1) { setSeed(seed); }

    void setSeed(u64 seed) {
        u64 x = seed;
        for (int i = 0; i < 4; i++) s[i] = splitmix(x);
    }

    static inline u64 rotl(u64 x, int k) { return (x << k) | (x >> (64 - k)); }

    inline u64 next() {
        u64 r = rotl(s[1] * 5, 7) * 9, t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = rotl(s[3], 45);
        return r;
    }

    inline double randomDouble() { return (next() >> 11) * 0x1.0p-53; }

    inline double randomDoubleOpen0() { return ((next() >> 11) + 1) * 0x1.0p-53; }

    inline u64 randomBelow(u64 n) { return (u64)(((unsigned __int128)next() * n) >> 64); }
};

inline u64 mix64(u64 a, u64 b = 0, u64 c = 0) {
    u64 x = a * 0x9E3779B97F4A7C15ULL ^ (b + 0x632BE59BD9B4E019ULL) * 0xC2B2AE3D27D4EB4FULL ^
            (c + 0x165667B19E3779F9ULL) * 0x27D4EB2F165667C5ULL;
    return Rng::splitmix(x);
}

string __head_version;

class Argument {
public:
    string dataset;
    string algo;
    string func = "eq1";
    string res;
    u32 k = 100;
    u32 b = 5;
    double epsilon = 0.1;
    double gamma = 0.5;
    int T = 5;
    u64 realization = 0;
    u64 seed = 0;
    bool seed_given = false;
    u64 theta = 100000;
    double bv_frac = 0.1;
    double bv_low = 0.1;
};

Argument args;

void parseArgs(int argn, char** argv) {
    for (int i = 1; i < argn; i++) {
        string a = string(argv[i]);
        if (i + 1 >= argn) {
            cerr << "Missing value for " << a << endl;
            exit(1);
        }
        if (a == "-dataset") args.dataset = string(argv[++i]);
        else if (a == "-algo") args.algo = string(argv[++i]);
        else if (a == "-func") args.func = string(argv[++i]);
        else if (a == "-k") args.k = atoi(argv[++i]);
        else if (a == "-b") args.b = atoi(argv[++i]);
        else if (a == "-epsilon") args.epsilon = atof(argv[++i]);
        else if (a == "-gamma") args.gamma = atof(argv[++i]);
        else if (a == "-T") args.T = atoi(argv[++i]);
        else if (a == "-realization") args.realization = strtoull(argv[++i], nullptr, 10);
        else if (a == "-seed") { args.seed = strtoull(argv[++i], nullptr, 10); args.seed_given = true; }
        else if (a == "-theta") args.theta = strtoull(argv[++i], nullptr, 10);
        else if (a == "-bvFrac") args.bv_frac = atof(argv[++i]);
        else if (a == "-bvLow") args.bv_low = atof(argv[++i]);
        else if (a == "-res") args.res = string(argv[++i]);
        else {
            cerr << "Unknown argument: " << a << endl;
            exit(1);
        }
    }
    if (!args.seed_given) args.seed = 1000 + args.realization;
}

class OutputInfo {
public:
    OutputInfo(int argn, char** argv) {
        cout << "======================================" << endl;
        cout << "Adaptive Adoption Maximization via Objective Decomposition" << endl;
        cout << "Version: " << __head_version << endl;
        cout << "Arguments: ";
        for (int i = 0; i < argn; i++) cout << argv[i] << " ";
        cout << endl;
        cout << "======================================" << endl;
    }
};

#endif
