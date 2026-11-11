#include "AAMS.h"

void Run(int argn, char** argv) {
    parseArgs(argn, argv);

    const vector<string> algos = {"AdaMODE", "AdaMODE-", "MODE", "HLA", "AdaIM", "AdaDegree", "AdaMODE^s", "AdaMODE^m"};
    if (find(algos.begin(), algos.end(), args.algo) == algos.end()) {
        cerr << "Unknown algorithm: " << args.algo << endl;
        cerr << "Available: AdaMODE, AdaMODE-, MODE, HLA, AdaIM, AdaDegree, AdaMODE^s, AdaMODE^m" << endl;
        exit(1);
    }
    if (args.algo == "AdaMODE^s" || args.algo == "AdaMODE^m") args.func = "threshold";
    if (args.func != "eq1" && args.func != "threshold") {
        cerr << "Unknown adoption function: " << args.func << " (available: eq1, threshold)" << endl;
        exit(1);
    }

    cout << "Dataset: " << args.dataset << endl;
    cout << "Algorithm: " << args.algo << endl;
    cout << "Adoption function: " << args.func << endl;
    cout << "k: " << args.k << ", b: " << args.b << endl;
    cout << "Epsilon: " << args.epsilon << endl;
    cout << "Realization: " << args.realization << ", seed: " << args.seed << endl;

    Graph g(args.dataset);
    generateInterest(g, args.bv_frac, args.bv_low, 1);

    AdoptionFunction fe = AdoptionFunction::parse(args.func == "eq1" ? "eq1" : "thresh", args.gamma, args.T);
    string opt = args.func == "eq1" ? "eq1" : "thresh_trunc";
    if (args.func == "threshold" && args.algo == "AdaMODE^m") opt = "thresh_mod";
    if (args.func == "threshold" && args.algo == "HLA") opt = "thresh";
    if (args.algo == "AdaIM") opt = "infl";
    AdoptionFunction fo = AdoptionFunction::parse(opt, args.gamma, args.T);

    AAMS aams(g, args.algo, fe, fo, args.k, args.b, args.epsilon, args.seed);
    aams.run(args.realization);
    aams.printResults();

    if (!args.res.empty()) {
        aams.saveResults(args.res);
    }
}

int main(int argn, char** argv) {
    __head_version = "v0.1";
    OutputInfo info(argn, argv);
    Run(argn, argv);
    return 0;
}
