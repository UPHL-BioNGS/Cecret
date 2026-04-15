#!/usr/bin/env python3

import argparse
import csv
from collections import defaultdict
from statistics import median

# metadata csv file should look something like this
# sample_id,cecret_fasta,viralrecon_fasta,artic_fasta,lineage_cecret,lineage_viralrecon,lineage_artic,sample_type
# S1,cecret/S1.fasta,viralrecon/S1.fasta,artic/S1.fasta,BA.2,BA.2,BA.2,high
# S2,cecret/S2.fasta,viralrecon/S2.fasta,artic/S2.fasta,BA.5,BA.5,BA.4,low

# usage
# python benchmarking.py --metadata samples.csv

VALID_BASES = set("ACGT")


def read_fasta(path):
    """Reads a single-sequence FASTA file."""
    seq = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq.append(line.strip().upper())
    return "".join(seq)


def snp_difference(seq1, seq2):
    """Count SNP differences ignoring Ns and gaps."""
    snps = 0
    for b1, b2 in zip(seq1, seq2):
        if b1 in VALID_BASES and b2 in VALID_BASES:
            if b1 != b2:
                snps += 1
    return snps


def percent_ns(seq):
    """Calculate % Ns over full sequence length."""
    n_count = seq.count("N")
    return (n_count / len(seq)) * 100 if len(seq) > 0 else 0


def load_metadata(path):
    """
    Expected CSV format:
    sample_id,cecret_fasta,viralrecon_fasta,artic_fasta,lineage_cecret,lineage_viralrecon,lineage_artic,sample_type
    """
    samples = []
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            samples.append(row)
    return samples


def compute_metrics(samples):
    results = []

    for s in samples:
        seq_c = read_fasta(s["cecret_fasta"])
        seq_v = read_fasta(s["viralrecon_fasta"])
        seq_a = read_fasta(s["artic_fasta"])

        snp_cv = snp_difference(seq_c, seq_v)
        snp_ca = snp_difference(seq_c, seq_a)

        ns_c = percent_ns(seq_c)
        ns_v = percent_ns(seq_v)
        ns_a = percent_ns(seq_a)

        concord_v = int(s["lineage_cecret"] == s["lineage_viralrecon"])
        concord_a = int(s["lineage_cecret"] == s["lineage_artic"])

        results.append({
            "sample_id": s["sample_id"],
            "sample_type": s["sample_type"],
            "snp_cv": snp_cv,
            "snp_ca": snp_ca,
            "ns_c": ns_c,
            "ns_v": ns_v,
            "ns_a": ns_a,
            "concord_v": concord_v,
            "concord_a": concord_a
        })

    return results


def summarize(values):
    return {
        "median": median(values),
        "min": min(values),
        "max": max(values)
    }


def overall_summary(results):
    snp_cv = [r["snp_cv"] for r in results]
    snp_ca = [r["snp_ca"] for r in results]

    ns_c = [r["ns_c"] for r in results]
    ns_v = [r["ns_v"] for r in results]
    ns_a = [r["ns_a"] for r in results]

    concord_v = sum(r["concord_v"] for r in results) / len(results) * 100
    concord_a = sum(r["concord_a"] for r in results) / len(results) * 100

    return {
        "n": len(results),
        "snp_cv": summarize(snp_cv),
        "snp_ca": summarize(snp_ca),
        "ns_c": summarize(ns_c),
        "ns_v": summarize(ns_v),
        "ns_a": summarize(ns_a),
        "concord_v": concord_v,
        "concord_a": concord_a
    }


def stratified_summary(results):
    grouped = defaultdict(list)
    for r in results:
        grouped[r["sample_type"]].append(r)

    summaries = {}

    for group, items in grouped.items():
        summaries[group] = overall_summary(items)

    return summaries


def print_overall(summary):
    print("\n=== Overall Summary (Table 1) ===")
    print(f"Samples: {summary['n']}")
    print(f"Median SNP diff (vs viralrecon): {summary['snp_cv']['median']} "
          f"(range {summary['snp_cv']['min']}-{summary['snp_cv']['max']})")
    print(f"Median SNP diff (vs ARTIC): {summary['snp_ca']['median']} "
          f"(range {summary['snp_ca']['min']}-{summary['snp_ca']['max']})")
    print(f"Lineage concordance vs viralrecon: {summary['concord_v']:.1f}%")
    print(f"Lineage concordance vs ARTIC: {summary['concord_a']:.1f}%")
    print(f"Median % Ns (CECRET): {summary['ns_c']['median']:.2f}")
    print(f"Median % Ns (viralrecon): {summary['ns_v']['median']:.2f}")
    print(f"Median % Ns (ARTIC): {summary['ns_a']['median']:.2f}")


def print_stratified(summaries):
    print("\n=== Stratified Summary (Table 2) ===")
    for group, s in summaries.items():
        print(f"\n[{group}]")
        print(f"  N = {s['n']}")
        print(f"  Median SNP diff (vs viralrecon): {s['snp_cv']['median']}")
        print(f"  Median SNP diff (vs ARTIC): {s['snp_ca']['median']}")
        print(f"  Median % Ns (CECRET): {s['ns_c']['median']:.2f}")
        print(f"  Median % Ns (viralrecon): {s['ns_v']['median']:.2f}")
        print(f"  Median % Ns (ARTIC): {s['ns_a']['median']:.2f}")
        print(f"  Concordance vs viralrecon: {s['concord_v']:.1f}%")
        print(f"  Concordance vs ARTIC: {s['concord_a']:.1f}%")


def main():
    parser = argparse.ArgumentParser(description="Benchmark CECRET vs viralrecon and ARTIC")
    parser.add_argument("--metadata", required=True, help="CSV file with sample metadata")
    args = parser.parse_args()

    samples = load_metadata(args.metadata)
    results = compute_metrics(samples)

    overall = overall_summary(results)
    stratified = stratified_summary(results)

    print_overall(overall)
    print_stratified(stratified)


if __name__ == "__main__":
    main()
