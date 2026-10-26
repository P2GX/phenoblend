#!/usr/bin/env python3
"""Counts how many phenopackets show enrichment in a phenoblend-batch table.

usage: python3 summarize_batch.py results.tsv [--mode observed|shuffled] [--one-per-study]
"""
import csv
import math
import random
import sys
from collections import defaultdict


def read_rows(path):
    with open(path, newline="") as f:
        lines = [line for line in f if not line.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def binomial_tail(n, p, k):
    """P(X >= k) for X ~ Binomial(n, p)."""
    return sum(math.comb(n, i) * p**i * (1 - p) ** (n - i) for i in range(k, n + 1))


def study_of(row):
    """The publication: PMID_12345678_... -> PMID_12345678."""
    return "_".join(row["phenopacket_id"].split("_")[:2])


def main():
    path = sys.argv[1]
    mode = "observed"
    if "--mode" in sys.argv:
        mode = sys.argv[sys.argv.index("--mode") + 1]
    one_per_study = "--one-per-study" in sys.argv

    rows = read_rows(path)
    status = defaultdict(int)
    for r in rows:
        status[(r["mode"], r["status"])] += 1
    print("rows by (mode, status):", dict(status))

    ok = [r for r in rows if r["mode"] == mode and r["status"] == "ok"]

    # one row per phenopacket for the main analysis: the pair test of two-gene individuals
    pairs = [r for r in ok if r["n_genes"] == "2"]
    by_ppkt = defaultdict(list)
    for r in pairs:
        by_ppkt[r["phenopacket_id"]].append(r)
    main_rows = [rs[0] for rs in by_ppkt.values() if len(rs) == 1]
    multi = {pid for pid, rs in by_ppkt.items() if len(rs) > 1}
    print(f"phenopackets with exactly two genes: {len(main_rows)}; with more genes: {len(multi)}")

    if one_per_study:
        best = {}
        for r in main_rows:
            key = (study_of(r), r["genes"])
            if key not in best or int(r["n_observed_terms"]) > int(best[key]["n_observed_terms"]):
                best[key] = r
        main_rows = list(best.values())
        print("after keeping one individual per (study, genes):", len(main_rows))

    n = len(main_rows)
    print(f"\nmain analysis ({mode}): {n} phenopackets with two genes")
    for column, label in (("p_value", "exact"), ("sim_p_value", "weighted simulation")):
        values = [float(r[column]) for r in main_rows if r[column] != ""]
        for alpha in (0.05, 0.01, 0.001):
            k = sum(v < alpha for v in values)
            print(f"  {label:20s} p < {alpha:<6}: {k:3d}  (expected at most {n * alpha:.1f}; "
                  f"P(X >= {k}) under the null = {binomial_tail(n, alpha, k):.1e})")

    informative = [r for r in main_rows if r["degenerate"] != "true"]
    print(f"  degenerate (no information) tests: {n - len(informative)} of {n}")

    # randomized p-values (exactly uniform under the null), one-sided KS for enrichment
    rng = random.Random(1)
    u = sorted(float(r["p_greater"]) + rng.random() * float(r["p_equal"]) for r in main_rows)
    d_plus = max((i + 1) / n - x for i, x in enumerate(u))
    print(f"\nrandomized exact p-values: KS D+ = {d_plus:.3f}, one-sided p ~ {math.exp(-2 * n * d_plus**2):.1e}")

    # three or more genes: the pairs of each individual
    if multi:
        print("\nindividuals with more than two genes (smallest pair p-value):")
        for pid in sorted(multi):
            best = min(by_ppkt[pid], key=lambda r: float(r["p_value"]))
            print(f"  {pid}: {best['genes']} p = {float(best['p_value']):.3g}")


if __name__ == "__main__":
    main()