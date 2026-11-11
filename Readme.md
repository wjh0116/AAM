# AdaMODE: Adaptive Adoption Maximization via Objective Decomposition

## Build

Requires a C++17 compiler.

```bash
make        # produces ./AAM
```

## Usage

```bash
./AAM -dataset <path>/ -algo <name> -k <#seeds> -b <batch size> [options]
```

Examples (on the included toy dataset):

```bash
./AAM -dataset datasets/toy/ -algo AdaMODE -k 10 -b 2 -epsilon 0.1
./AAM -dataset datasets/toy/ -algo HLA -k 10 -theta 100000
./AAM -dataset datasets/toy/ -algo "AdaMODE^s" -k 10 -b 2
```

| Argument         | Description                                                                                     |
| ---------------- | ----------------------------------------------------------------------------------------------- |
| `-dataset`     | Dataset folder path                                                                             |
| `-algo`        | `AdaMODE`, `AdaMODE-`, `MODE`, `HLA`, `AdaIM`, `AdaDegree`, `AdaMODE^s`, `AdaMODE^m` |
| `-k`           | Number of seeds (default 100)                                                                   |
| `-b`           | Batch size (default 5)                                                                          |
| `-epsilon`     | Error tolerance ε (default 0.1)                                                                 |
| `-func`        | Adoption function: `eq1` (default) or `threshold`                                           |
| `-gamma`       | `eq1`: parameter γ (default 0.5)                                                              |
| `-T`           | `threshold`: threshold (default 5)                                                            |
| `-theta`       | HLA: sample size (default 100000)                                                               |
| `-realization` | Realization ID (default 0)                                                                      |
| `-seed`        | Random seed (default 1000 + realization ID)                                                     |
| `-bvFrac`      | Fraction of nodes with interest b_v = 1 (default 0.1)                                           |
| `-bvLow`       | Interest b_v of the other nodes (default 0.1)                                                   |
| `-res`         | Result output file path (optional)                                                              |

## Dataset Format

Each dataset folder contains (see `datasets/toy/` for an example):

- `attribute.txt` — two lines: `n=<#nodes>` and `m=<#edges>`
- `graph.txt` — one directed edge per line: `source target probability` (an optional `n m` header line is auto-detected and skipped)
