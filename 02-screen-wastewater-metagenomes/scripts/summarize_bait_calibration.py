#!/usr/bin/env python3
"""Derive the bait and calibration measurements used in the Discussion.

The default invocation discovers the accepted stage 01 bait evidence and stage
02 calibration evidence from this repository and writes a deterministic JSON
summary. It needs no BLAST database, sequencing archive, cluster, or network.

The report intentionally excludes the superseded threshold-20 production reads
and live accession lookups. Those have different provenance from the committed
calibration evidence summarized here.
"""

from __future__ import annotations

import collections
import csv
import gzip
import json
import math
import sys
from pathlib import Path
from typing import TextIO


KMER_LENGTH = 31
CALIBRATION_WINDOW = 150
VALID_CLASSES = ("target", "non_target", "top_tie", "no_hit")
NON_TARGET_CLASSES = frozenset(("non_target", "top_tie", "no_hit"))
TIE_EPS = 0.1
COMPLEMENT = str.maketrans("ACGT", "TGCA")


class AnalysisError(ValueError):
    """The committed evidence does not satisfy the analysis contract."""


def open_text(path: Path) -> TextIO:
    if not path.is_file():
        raise AnalysisError(f"missing input: {path}")
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="")
    return path.open(newline="")


def canonical(sequence: str) -> str:
    normalized = sequence.upper()
    if set(normalized) - set("ACGT"):
        raise AnalysisError(f"non-ACGT bait sequence: {sequence!r}")
    reverse = normalized.translate(COMPLEMENT)[::-1]
    return min(normalized, reverse)


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []
    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line.startswith(">"):
                if name is not None:
                    if name in records:
                        raise AnalysisError(f"duplicate FASTA record {name!r} in {path}")
                    records[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            elif line:
                if name is None:
                    raise AnalysisError(f"sequence before first FASTA header in {path}")
                chunks.append(line)
    if name is not None:
        if name in records:
            raise AnalysisError(f"duplicate FASTA record {name!r} in {path}")
        records[name] = "".join(chunks).upper()
    if not records:
        raise AnalysisError(f"no FASTA records in {path}")
    return records


def read_tsv(path: Path, required: set[str]) -> list[dict[str, str]]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise AnalysisError(f"no TSV header in {path}")
        missing = required - set(reader.fieldnames)
        if missing:
            raise AnalysisError(f"{path} is missing columns: {sorted(missing)}")
        rows = list(reader)
    if not rows:
        raise AnalysisError(f"no rows in {path}")
    return rows


def parse_nonnegative_int(value: str, field: str, path: Path) -> int:
    try:
        parsed = int(value)
    except ValueError as error:
        raise AnalysisError(f"invalid {field}={value!r} in {path}") from error
    if parsed < 0:
        raise AnalysisError(f"negative {field}={value!r} in {path}")
    return parsed


def bait_matches(sequence: str, bait_sequences: set[str]) -> set[str]:
    matches: set[str] = set()
    normalized = sequence.upper()
    for start in range(len(normalized) - KMER_LENGTH + 1):
        window = normalized[start : start + KMER_LENGTH]
        if set(window) <= set("ACGT"):
            candidate = min(window, window.translate(COMPLEMENT)[::-1])
            if candidate in bait_sequences:
                matches.add(candidate)
    return matches


def quantile(values: list[float] | list[int], probability: float) -> float | None:
    """Type-7/linear quantile: interpolate at p * (n - 1)."""
    if not values:
        return None
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    low = math.floor(position)
    high = math.ceil(position)
    if low == high:
        return round(float(ordered[low]), 6)
    value = ordered[low] * (high - position) + ordered[high] * (position - low)
    return round(float(value), 6)


def fraction(numerator: int, denominator: int) -> float | None:
    return round(numerator / denominator, 6) if denominator else None


def population_gini(values: list[int]) -> float:
    ordered = sorted(values)
    total = sum(ordered)
    if not ordered or total == 0:
        return 0.0
    weighted = sum(index * value for index, value in enumerate(ordered, start=1))
    result = 2 * weighted / (len(ordered) * total) - (len(ordered) + 1) / len(ordered)
    return round(result, 6)


def histogram_summary(histogram: collections.Counter[int]) -> dict[str, object]:
    expanded = [value for value, count in histogram.items() for _ in range(count)]
    total = sum(histogram.values())
    return {
        "n": total,
        "minimum": min(expanded) if expanded else None,
        "median": quantile(expanded, 0.5),
        "maximum": max(expanded) if expanded else None,
        "exactly_one": histogram.get(1, 0),
        "exactly_one_fraction": fraction(histogram.get(1, 0), total),
    }


def usage_summary(
    usage: collections.Counter[str],
    bait_sequence_to_id: dict[str, str],
) -> dict[str, object]:
    values = [usage.get(sequence, 0) for sequence in bait_sequence_to_id]
    ranked = sorted(
        (
            (count, bait_sequence_to_id[sequence])
            for sequence, count in usage.items()
            if count > 0
        ),
        key=lambda item: (-item[0], item[1]),
    )
    total = sum(values)
    return {
        "total_distinct_bait_read_incidences": total,
        "unused_baits": sum(value == 0 for value in values),
        "baits_used_once": sum(value == 1 for value in values),
        "median_observations_per_bait": quantile(values, 0.5),
        "population_gini_including_unused_baits": population_gini(values),
        "top_10_share": fraction(sum(count for count, _ in ranked[:10]), total),
    }


def select_threshold(observations: list[tuple[str, str, int, str]]) -> int | None:
    maximum = max((hits for _, _, hits, _ in observations), default=0)
    for threshold in range(1, maximum + 1):
        retained = collections.Counter(
            read_class
            for _, _, hits, read_class in observations
            if hits >= threshold
        )
        if retained["target"] > 0 and all(retained[name] == 0 for name in NON_TARGET_CLASSES):
            return threshold
    return None


def contiguous_runs(positions: list[tuple[int, str]]) -> list[list[tuple[int, str]]]:
    if not positions:
        return []
    ordered = sorted(positions)
    runs = [[ordered[0]]]
    for position, bait in ordered[1:]:
        if position != runs[-1][-1][0] + 1:
            runs.append([])
        runs[-1].append((position, bait))
    return runs


def summarize_topology(
    target_fasta: Path,
    bait_sequences: set[str],
    threshold: int,
) -> list[dict[str, object]]:
    summaries = []
    for target, sequence in read_fasta(target_fasta).items():
        positions = []
        for start in range(len(sequence) - KMER_LENGTH + 1):
            window = sequence[start : start + KMER_LENGTH]
            if set(window) <= set("ACGT"):
                candidate = min(window, window.translate(COMPLEMENT)[::-1])
                if candidate in bait_sequences:
                    positions.append((start, candidate))
        runs = contiguous_runs(positions)
        run_bait_counts = [len({bait for _, bait in run}) for run in runs]
        run_spans = [run[-1][0] - run[0][0] + KMER_LENGTH for run in runs]
        possible_windows = max(len(sequence) - CALIBRATION_WINDOW + 1, 0)
        qualifying_windows = 0
        if possible_windows:
            qualifying_windows = sum(
                len(
                    bait_matches(
                        sequence[start : start + CALIBRATION_WINDOW],
                        bait_sequences,
                    )
                )
                >= threshold
                for start in range(possible_windows)
            )
        summaries.append(
            {
                "target": target,
                "target_length_nt": len(sequence),
                "matching_start_positions": len(positions),
                "distinct_baits": len({bait for _, bait in positions}),
                "contiguous_runs": len(runs),
                "runs_with_at_least_threshold_baits": sum(
                    count >= threshold for count in run_bait_counts
                ),
                "longest_run_baits": max(run_bait_counts, default=0),
                "longest_run_span_nt": max(run_spans, default=0),
                "possible_150nt_windows": possible_windows,
                "qualifying_150nt_windows": qualifying_windows,
                "qualifying_150nt_window_fraction": fraction(
                    qualifying_windows,
                    possible_windows,
                ),
            }
        )
    return summaries


def update_best_hsp(
    records: dict[str, dict[str, dict[str, object]]],
    query: str,
    side: str,
    candidate: dict[str, object],
) -> None:
    current = records.setdefault(query, {}).get(side)
    candidate_rank = (
        candidate["bitscore"],
        candidate["query_coverage"],
        candidate["alignment_length"],
        candidate["identity"],
        candidate["accession"],
    )
    if current is None:
        records[query][side] = candidate
        return
    current_rank = (
        current["bitscore"],
        current["query_coverage"],
        current["alignment_length"],
        current["identity"],
        current["accession"],
    )
    if candidate_rank > current_rank:
        records[query][side] = candidate


def load_best_hsps(
    blast_path: Path,
    target_taxids: set[str],
) -> tuple[
    dict[str, dict[str, dict[str, object]]],
    collections.Counter[str],
    dict[str, set[str]],
]:
    best: dict[str, dict[str, dict[str, object]]] = {}
    row_counts: collections.Counter[str] = collections.Counter()
    accessions: dict[str, set[str]] = collections.defaultdict(set)
    with open_text(blast_path) as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row or row[0] == "qseqid":
                continue
            if len(row) < 16:
                raise AnalysisError(f"short BLAST row for query {row[0]!r} in {blast_path}")
            query = row[0]
            try:
                identity = float(row[4])
                alignment_length = int(row[5])
                bitscore = float(row[13])
                query_coverage = float(row[14])
            except ValueError as error:
                raise AnalysisError(f"invalid numeric BLAST field for query {query!r}") from error
            accession = row[2]
            taxids = {
                value
                for value in row[3].replace(",", ";").split(";")
                if value and value != "N/A"
            }
            candidate: dict[str, object] = {
                "accession": accession,
                "identity": identity,
                "alignment_length": alignment_length,
                "bitscore": bitscore,
                "query_coverage": query_coverage,
                "title": row[15],
            }
            row_counts[query] += 1
            accessions[query].add(accession)
            if taxids & target_taxids:
                update_best_hsp(best, query, "target", candidate)
            if taxids - target_taxids or not taxids:
                update_best_hsp(best, query, "other", candidate)
    return best, row_counts, accessions


def classify_hsps(record: dict[str, dict[str, object]] | None) -> str:
    if not record:
        return "no_hit"
    target = record.get("target")
    other = record.get("other")
    if target is None:
        return "non_target"
    if other is None:
        return "target"
    target_bits = float(target["bitscore"])
    other_bits = float(other["bitscore"])
    if abs(target_bits - other_bits) < TIE_EPS:
        return "top_tie"
    return "target" if target_bits > other_bits else "non_target"


def summarize_alignment_group(
    representatives: list[str],
    best_hsps: dict[str, dict[str, dict[str, object]]],
    representative_weights: collections.Counter[str],
    weighted: bool,
) -> dict[str, object]:
    identities: list[float] = []
    lengths: list[int] = []
    coverages: list[float] = []
    margins: list[float] = []
    no_competitor = 0
    for representative in representatives:
        weight = representative_weights[representative] if weighted else 1
        target = best_hsps[representative]["target"]
        identities.extend([float(target["identity"])] * weight)
        lengths.extend([int(target["alignment_length"])] * weight)
        coverages.extend([float(target["query_coverage"])] * weight)
        other = best_hsps[representative].get("other")
        if other is None:
            no_competitor += weight
        else:
            margin = float(target["bitscore"]) - float(other["bitscore"])
            margins.extend([margin] * weight)
    return {
        "n": len(identities),
        "minimum_identity": min(identities) if identities else None,
        "median_identity": quantile(identities, 0.5),
        "maximum_identity": max(identities) if identities else None,
        "identity_at_least_99_fraction": fraction(
            sum(value >= 99 for value in identities),
            len(identities),
        ),
        "identity_100_fraction": fraction(
            sum(value == 100 for value in identities),
            len(identities),
        ),
        "median_alignment_length": quantile(lengths, 0.5),
        "minimum_query_coverage": min(coverages) if coverages else None,
        "median_query_coverage": quantile(coverages, 0.5),
        "maximum_query_coverage": max(coverages) if coverages else None,
        "query_coverage_at_least_95_fraction": fraction(
            sum(value >= 95 for value in coverages),
            len(coverages),
        ),
        "minimum_reported_score_margin": min(margins) if margins else None,
        "median_reported_score_margin": quantile(margins, 0.5),
        "maximum_reported_score_margin": max(margins) if margins else None,
        "no_reported_non_target_competitor": no_competitor,
    }


def build_summary(study_root: Path) -> dict[str, object]:
    stage1 = study_root / "01-identify-cyclospora-specific-kmers"
    calibration = study_root / "02-screen-wastewater-metagenomes" / "results" / "calibration"
    paths = {
        "bait_fasta": stage1 / "baits/cyclospora_cayetanensis_rrna_core_nt_validated_baits.fasta",
        "target_fasta": stage1 / "curated/target_rrna.fasta",
        "bait_manifest": stage1 / "results/cyclospora_cayetanensis_rrna_kmer_manifest.tsv",
        "core_nt_validation": stage1 / "results/cyclospora_cayetanensis_core_nt_validation.tsv",
        "calibration_queries": calibration / "calibration_read_blast_queries.fasta.gz",
        "read_map": calibration / "read_blast_map.tsv",
        "read_classification": calibration / "read_blast_deacon.tsv",
        "target_taxids": calibration / "cyclospora_genus_taxids.txt",
        "read_blast": calibration / "read_blast_hits.tsv.gz",
    }
    for path in paths.values():
        if not path.is_file():
            raise AnalysisError(f"missing input: {path}")

    bait_records = read_fasta(paths["bait_fasta"])
    bait_id_to_sequence: dict[str, str] = {}
    bait_sequence_to_id: dict[str, str] = {}
    for bait_id, sequence in bait_records.items():
        if len(sequence) != KMER_LENGTH:
            raise AnalysisError(
                f"validated bait {bait_id!r} is {len(sequence)} nt, expected {KMER_LENGTH}"
            )
        normalized = canonical(sequence)
        if normalized in bait_sequence_to_id:
            raise AnalysisError(
                f"duplicate canonical bait sequence for {bait_id!r} and "
                f"{bait_sequence_to_id[normalized]!r}"
            )
        bait_id_to_sequence[bait_id] = normalized
        bait_sequence_to_id[normalized] = bait_id
    bait_sequences = set(bait_sequence_to_id)

    manifest_rows = read_tsv(
        paths["bait_manifest"],
        {"kmer", "rrna_classes", "target_copy_count", "core_nt_bait_id", "core_nt_status"},
    )
    manifest_by_sequence: dict[str, dict[str, str]] = {}
    for row in manifest_rows:
        if row["core_nt_status"] != "PASS_CORE_NT":
            continue
        sequence = canonical(row["kmer"])
        if sequence in manifest_by_sequence:
            raise AnalysisError(f"duplicate retained manifest k-mer {sequence}")
        manifest_by_sequence[sequence] = row
    if set(manifest_by_sequence) != bait_sequences:
        raise AnalysisError("retained manifest k-mers do not equal the validated bait FASTA")
    for sequence, row in manifest_by_sequence.items():
        if row["core_nt_bait_id"] != bait_sequence_to_id[sequence]:
            raise AnalysisError(f"manifest bait ID disagrees for {row['core_nt_bait_id']!r}")

    classifications: dict[str, tuple[str, int, str]] = {}
    for row in read_tsv(
        paths["read_classification"],
        {"read_id", "sample", "deacon_hits", "blast_class"},
    ):
        read_id = row["read_id"]
        if not read_id:
            raise AnalysisError(f"blank read_id in {paths['read_classification']}")
        if read_id in classifications:
            raise AnalysisError(f"duplicate read_id {read_id!r} in read classification")
        read_class = row["blast_class"]
        if read_class not in VALID_CLASSES:
            raise AnalysisError(f"unknown blast_class {read_class!r} for {read_id!r}")
        classifications[read_id] = (
            row["sample"],
            parse_nonnegative_int(row["deacon_hits"], "deacon_hits", paths["read_classification"]),
            read_class,
        )

    representative_sequences = read_fasta(paths["calibration_queries"])
    representative_matches = {
        representative: bait_matches(sequence, bait_sequences)
        for representative, sequence in representative_sequences.items()
    }
    observations: list[tuple[str, str, int, str]] = []
    representative_weights: collections.Counter[str] = collections.Counter()
    representative_classes: dict[str, str] = {}
    representative_hits: dict[str, int] = {}
    representative_samples: dict[str, collections.Counter[str]] = collections.defaultdict(collections.Counter)
    seen_read_ids: set[str] = set()
    for row in read_tsv(
        paths["read_map"],
        {"sample", "read_id", "representative_id", "deacon_hits", "read_length"},
    ):
        read_id = row["read_id"]
        if read_id in seen_read_ids:
            raise AnalysisError(f"duplicate read_id {read_id!r} in read map")
        seen_read_ids.add(read_id)
        if read_id not in classifications:
            raise AnalysisError(f"read {read_id!r} is missing from read classification")
        sample = row["sample"]
        hits = parse_nonnegative_int(row["deacon_hits"], "deacon_hits", paths["read_map"])
        class_sample, class_hits, read_class = classifications[read_id]
        if (sample, hits) != (class_sample, class_hits):
            raise AnalysisError(f"map/classification disagreement for read {read_id!r}")
        representative = row["representative_id"]
        if representative not in representative_matches:
            raise AnalysisError(f"representative {representative!r} is missing from query FASTA")
        recounted = len(representative_matches[representative])
        if recounted != hits:
            raise AnalysisError(
                f"bait recount mismatch for {representative!r}: evidence={hits}, recounted={recounted}"
            )
        if representative in representative_hits and representative_hits[representative] != hits:
            raise AnalysisError(f"inconsistent bait counts for representative {representative!r}")
        if representative in representative_classes and representative_classes[representative] != read_class:
            raise AnalysisError(f"inconsistent classes for representative {representative!r}")
        representative_hits[representative] = hits
        representative_classes[representative] = read_class
        representative_weights[representative] += 1
        representative_samples[representative][sample] += 1
        observations.append((representative, sample, hits, read_class))
    extra_classifications = set(classifications) - seen_read_ids
    if extra_classifications:
        raise AnalysisError(
            f"{len(extra_classifications)} classified reads are missing from the read map"
        )
    unused_representatives = set(representative_sequences) - set(representative_weights)
    if unused_representatives:
        raise AnalysisError(
            f"{len(unused_representatives)} query representatives are unused by the read map"
        )

    threshold = select_threshold(observations)
    if threshold is None:
        raise AnalysisError("calibration has no supported target-only threshold")
    boundary_samples = collections.Counter(
        sample
        for _, sample, hits, read_class in observations
        if hits == threshold - 1 and read_class in NON_TARGET_CLASSES
    )
    leave_one_boundary_sample = {
        sample: select_threshold(
            [observation for observation in observations if observation[1] != sample]
        )
        for sample in sorted(boundary_samples)
    }

    target_taxids = {
        line.split("#", 1)[0].strip()
        for line in paths["target_taxids"].read_text().splitlines()
        if line.split("#", 1)[0].strip()
    }
    if not target_taxids:
        raise AnalysisError(f"no target taxids in {paths['target_taxids']}")
    best_hsps, blast_row_counts, blast_accessions = load_best_hsps(
        paths["read_blast"],
        target_taxids,
    )
    for representative, read_class in representative_classes.items():
        recomputed = classify_hsps(best_hsps.get(representative))
        if recomputed != read_class:
            raise AnalysisError(
                f"BLAST recount class mismatch for {representative!r}: "
                f"evidence={read_class}, recomputed={recomputed}"
            )

    read_histogram: collections.Counter[int] = collections.Counter()
    class_histograms: dict[str, collections.Counter[int]] = {
        name: collections.Counter() for name in VALID_CLASSES
    }
    all_observation_usage: collections.Counter[str] = collections.Counter()
    all_unique_usage: collections.Counter[str] = collections.Counter()
    threshold_target_usage: collections.Counter[str] = collections.Counter()
    threshold_target_unique_usage: collections.Counter[str] = collections.Counter()
    for representative, _, hits, read_class in observations:
        read_histogram[hits] += 1
        class_histograms[read_class][hits] += 1
        matches = representative_matches[representative]
        all_observation_usage.update(matches)
        if read_class == "target" and hits >= threshold:
            threshold_target_usage.update(matches)
    for representative, matches in representative_matches.items():
        all_unique_usage.update(matches)
        if (
            representative_classes[representative] == "target"
            and representative_hits[representative] >= threshold
        ):
            threshold_target_unique_usage.update(matches)

    threshold_target_representatives = sorted(
        representative
        for representative, read_class in representative_classes.items()
        if read_class == "target" and representative_hits[representative] >= threshold
    )
    exact_threshold_representatives = [
        representative
        for representative in threshold_target_representatives
        if representative_hits[representative] == threshold
    ]
    threshold_target_observations = sum(
        representative_weights[representative]
        for representative in threshold_target_representatives
    )

    highest_usage_sequence = min(
        all_observation_usage,
        key=lambda sequence: (
            -all_observation_usage[sequence],
            bait_sequence_to_id[sequence],
        ),
    )
    highest_usage_by_class = collections.Counter()
    highest_usage_unique_by_class = collections.Counter()
    highest_usage_samples = collections.Counter()
    highest_usage_unique_samples: dict[str, set[str]] = collections.defaultdict(set)
    support_titles = collections.Counter()
    for representative, matches in representative_matches.items():
        if highest_usage_sequence not in matches:
            continue
        read_class = representative_classes[representative]
        weight = representative_weights[representative]
        highest_usage_by_class[read_class] += weight
        highest_usage_unique_by_class[read_class] += 1
        for sample, count in representative_samples[representative].items():
            highest_usage_samples[sample] += count
            highest_usage_unique_samples[sample].add(representative)
        other = best_hsps.get(representative, {}).get("other")
        if read_class == "non_target" and other is not None:
            label = f"{other['accession']} | {other['title']}"
            support_titles[label] += weight

    validation_rows = read_tsv(
        paths["core_nt_validation"],
        {"bait_id", "status", "exact_target_count", "exact_non_target_count"},
    )
    validation_by_id: dict[str, dict[str, str]] = {}
    for row in validation_rows:
        bait_id = row["bait_id"]
        if bait_id in validation_by_id:
            raise AnalysisError(f"duplicate validation bait_id {bait_id!r}")
        validation_by_id[bait_id] = row
    if not set(bait_id_to_sequence) <= set(validation_by_id):
        raise AnalysisError("validated bait FASTA contains IDs absent from validation table")
    absence_only_ids = {
        bait_id
        for bait_id, row in validation_by_id.items()
        if row["status"] == "PASS_CORE_NT"
        and parse_nonnegative_int(row["exact_target_count"], "exact_target_count", paths["core_nt_validation"]) == 0
        and parse_nonnegative_int(row["exact_non_target_count"], "exact_non_target_count", paths["core_nt_validation"]) == 0
    }
    if not absence_only_ids <= set(bait_id_to_sequence):
        raise AnalysisError("absence-only passers are missing from the validated bait FASTA")
    zero_target_reject_ids = {
        bait_id
        for bait_id, row in validation_by_id.items()
        if row["status"] != "PASS_CORE_NT"
        and parse_nonnegative_int(
            row["exact_target_count"],
            "exact_target_count",
            paths["core_nt_validation"],
        )
        == 0
        and parse_nonnegative_int(
            row["exact_non_target_count"],
            "exact_non_target_count",
            paths["core_nt_validation"],
        )
        > 0
    }
    absence_only_sequences = {bait_id_to_sequence[bait_id] for bait_id in absence_only_ids}
    reduced_observations = [
        (
            representative,
            sample,
            len(representative_matches[representative] - absence_only_sequences),
            read_class,
        )
        for representative, sample, _, read_class in observations
    ]
    reduced_threshold = select_threshold(reduced_observations)
    reduced_targets_at_threshold = (
        sum(
            read_class == "target" and hits >= reduced_threshold
            for _, _, hits, read_class in reduced_observations
        )
        if reduced_threshold is not None
        else 0
    )
    original_threshold_targets_still_passing = sum(
        read_class == "target"
        and original_hits >= threshold
        and reduced_hits >= threshold
        for (
            (_, _, original_hits, read_class),
            (_, _, reduced_hits, _),
        ) in zip(observations, reduced_observations, strict=True)
    )
    absence_threshold_incidences = sum(
        threshold_target_usage.get(sequence, 0)
        for sequence in absence_only_sequences
    )

    gene_counts = collections.Counter(
        row["rrna_classes"] for row in manifest_by_sequence.values()
    )
    target_copy_counts = collections.Counter(
        parse_nonnegative_int(
            row["target_copy_count"],
            "target_copy_count",
            paths["bait_manifest"],
        )
        for row in manifest_by_sequence.values()
    )

    summary: dict[str, object] = {
        "schema_version": 1,
        "parameters": {
            "kmer_length": KMER_LENGTH,
            "calibration_window_nt": CALIBRATION_WINDOW,
            "selected_threshold": threshold,
            "span_of_consecutive_threshold_kmers_nt": KMER_LENGTH + threshold - 1,
            "shared_intersection_of_consecutive_threshold_kmers_nt": (
                KMER_LENGTH - threshold + 1
            ),
            "bait_incidence_unit": "one distinct bait identity per read observation",
            "quantile_method": "linear interpolation at p * (n - 1)",
        },
        "bait_set": {
            "baits": len(bait_sequences),
            "by_rrna": {
                name: gene_counts[name]
                for name in ("18S", "28S", "5S", "5.8S")
            },
            "target_copy_count": {
                str(copy_count): count
                for copy_count, count in sorted(target_copy_counts.items())
            },
            "target_topology": summarize_topology(
                paths["target_fasta"],
                bait_sequences,
                threshold,
            ),
        },
        "calibration": {
            "read_observations": len(observations),
            "unique_sequences": len(representative_sequences),
            "reads_with_one_bait": read_histogram[1],
            "reads_with_one_bait_fraction": fraction(read_histogram[1], len(observations)),
            "read_bait_count": histogram_summary(read_histogram),
            "read_bait_count_by_class": {
                name: histogram_summary(class_histograms[name])
                for name in VALID_CLASSES
            },
            "threshold_target_observations": threshold_target_observations,
            "threshold_target_unique_sequences": len(threshold_target_representatives),
            "threshold_boundary": {
                "non_target_or_ambiguous_observations_at_threshold_minus_one": sum(
                    boundary_samples.values()
                ),
                "samples": dict(sorted(boundary_samples.items())),
                "selected_threshold_without_each_boundary_sample": leave_one_boundary_sample,
            },
            "exact_threshold_target_observations": sum(
                representative_weights[representative]
                for representative in exact_threshold_representatives
            ),
            "exact_threshold_target_unique_sequences": len(exact_threshold_representatives),
            "bait_usage": {
                "all_observations": usage_summary(
                    all_observation_usage,
                    bait_sequence_to_id,
                ),
                "all_unique_sequences": usage_summary(
                    all_unique_usage,
                    bait_sequence_to_id,
                ),
                "threshold_target_observations": usage_summary(
                    threshold_target_usage,
                    bait_sequence_to_id,
                ),
                "threshold_target_unique_sequences": usage_summary(
                    threshold_target_unique_usage,
                    bait_sequence_to_id,
                ),
                "most_used_bait": {
                    "bait_id": bait_sequence_to_id[highest_usage_sequence],
                    "observations_by_class": {
                        name: highest_usage_by_class[name]
                        for name in VALID_CLASSES
                    },
                    "unique_sequences_by_class": {
                        name: highest_usage_unique_by_class[name]
                        for name in VALID_CLASSES
                    },
                    "samples": len(highest_usage_samples),
                    "largest_sample_by_observations": highest_usage_samples.most_common(1)[0],
                    "largest_sample_by_unique_sequences": sorted(
                        (
                            (sample, len(representatives))
                            for sample, representatives in highest_usage_unique_samples.items()
                        ),
                        key=lambda item: (-item[1], item[0]),
                    )[0],
                    "top_non_target_alignment_support": support_titles.most_common(3),
                },
            },
            "target_alignment_support": {
                "threshold_observations": summarize_alignment_group(
                    threshold_target_representatives,
                    best_hsps,
                    representative_weights,
                    weighted=True,
                ),
                "threshold_unique_sequences": summarize_alignment_group(
                    threshold_target_representatives,
                    best_hsps,
                    representative_weights,
                    weighted=False,
                ),
                "exact_threshold_observations": summarize_alignment_group(
                    exact_threshold_representatives,
                    best_hsps,
                    representative_weights,
                    weighted=True,
                ),
                "exact_threshold_unique_sequences": summarize_alignment_group(
                    exact_threshold_representatives,
                    best_hsps,
                    representative_weights,
                    weighted=False,
                ),
                "threshold_queries_with_at_least_100_accessions": sum(
                    len(blast_accessions[representative]) >= 100
                    for representative in threshold_target_representatives
                ),
                "threshold_queries_with_at_least_100_rows": sum(
                    blast_row_counts[representative] >= 100
                    for representative in threshold_target_representatives
                ),
                "search_reporting_limit": "BLASTN used -max_target_seqs 100; competitor margins use the best reported HSPs",
            },
        },
        "database_support": {
            "absence_only_passers": len(absence_only_ids),
            "absence_only_fraction_of_baits": fraction(
                len(absence_only_ids),
                len(bait_sequences),
            ),
            "absence_only_by_rrna": dict(
                sorted(
                    collections.Counter(
                        manifest_by_sequence[bait_id_to_sequence[bait_id]]["rrna_classes"]
                        for bait_id in absence_only_ids
                    ).items()
                )
            ),
            "rejected_baits_without_exact_target_support": len(zero_target_reject_ids),
            "rejected_baits_without_exact_target_support_by_rrna": dict(
                sorted(
                    collections.Counter(
                        bait_id.rsplit("|", 1)[-1]
                        for bait_id in zero_target_reject_ids
                    ).items()
                )
            ),
            "absence_only_threshold_target_incidences": absence_threshold_incidences,
            "absence_only_fraction_of_threshold_target_incidences": fraction(
                absence_threshold_incidences,
                sum(threshold_target_usage.values()),
            ),
            "exclude_absence_only_sensitivity": {
                "selected_threshold": reduced_threshold,
                "target_observations_at_selected_threshold": reduced_targets_at_threshold,
                "original_threshold_targets_still_passing": original_threshold_targets_still_passing,
            },
        },
    }
    return summary


def main() -> int:
    study_root = Path(__file__).resolve().parents[2]
    output = (
        study_root
        / "02-screen-wastewater-metagenomes"
        / "results"
        / "calibration"
        / "bait_calibration_summary.json"
    )
    try:
        summary = build_summary(study_root)
    except (AnalysisError, OSError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp")
    temporary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    temporary.replace(output)
    print(f"selected threshold: {summary['parameters']['selected_threshold']}")
    print(f"wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
