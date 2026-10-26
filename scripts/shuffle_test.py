#!/usr/bin/env python3
"""Compares the real analysis with many shuffled controls.

usage: python3 shuffle_test.py observed.tsv shuffle_1.tsv shuffle_2.tsv ...
"""
import csv
import math
import random
import sys


def read_rows(path):
    with open(path, newline="") as f:
        lines = [line for line in f if not line.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def two_gene_rows(rows, mode):
    """One row per phenopacket that has exactly one two-gene result."""
    by_id = {}
    for r in rows:
        if r["mode"] == mode and r["status"] == "ok" and r["n_genes"] == "2":
            by_id.setdefault(r["phenopacket_id"], []).append(r)
    return [rs[0] for rs in by_id.values() if len(rs) == 1]


def statistics(rows, seed=1):
    n = len(rows)
    rng = random.Random(seed)
    u = sorted(float(r["p_greater"]) + rng.random() * float(r["p_equal"]) for r in rows)
    d_plus = max((i + 1) / n - x for i, x in enumerate(u))
    p = [float(r["p_value"]) for r in rows]
    return {
        "n": n,
        "below_0.05": sum(v < 0.05 for v in p),
        "below_0.01": sum(v < 0.01 for v in p),
        "d_plus": d_plus,
        "z": math.sqrt(n) * d_plus,  # comparable between tables of different size
        "degenerate": sum(r["degenerate"] == "true" for r in rows),
    }


def main():
    observed_path, shuffle_paths = sys.argv[1], sys.argv[2:]
    if not shuffle_paths:
        sys.exit("give at least one shuffled table")

    real = statistics(two_gene_rows(read_rows(observed_path), "observed"))
    shuffles = [statistics(two_gene_rows(read_rows(p), "shuffled"), seed=i + 2)
                for i, p in enumerate(shuffle_paths)]
    k = len(shuffles)

    print(f"real data: n = {real['n']}, p<0.05: {real['below_0.05']}, p<0.01: {real['below_0.01']}, "
          f"D+ = {real['d_plus']:.3f}, degenerate: {real['degenerate']}")
    print(f"{k} shuffles:")
    for key, label in (("below_0.05", "p<0.05"), ("below_0.01", "p<0.01"), ("d_plus", "D+"), ("degenerate", "degenerate")):
        values = [s[key] for s in shuffles]
        mean = sum(values) / k
        print(f"  {label:10s} mean {mean:7.3f}  min {min(values):7.3f}  max {max(values):7.3f}")

    print("\nhow often does a shuffle reach the real value? (empirical p = (1 + count) / (shuffles + 1))")
    for key, label in (("below_0.05", "number with p<0.05"), ("below_0.01", "number with p<0.01"), ("z", "sqrt(n) * D+")):
        count = sum(s[key] >= real[key] for s in shuffles)
        print(f"  {label:22s}: {count} of {k} shuffles; empirical p = {(1 + count) / (k + 1):.3f}")


if __name__ == "__main__":
    main()