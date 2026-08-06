from __future__ import annotations

import gzip
import hashlib
from collections.abc import Iterator
from pathlib import Path
from typing import TextIO

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def normalize_dna(sequence: str) -> str:
    return sequence.upper().replace("U", "T")


def reverse_complement(sequence: str) -> str:
    return normalize_dna(sequence).translate(_COMPLEMENT)[::-1]


def canonical_kmer(sequence: str) -> str:
    normalized = normalize_dna(sequence)
    if not set(normalized) <= {"A", "C", "G", "T"}:
        raise ValueError("canonical_kmer requires an ACGT-only sequence; found non-ACGT character")
    reverse = reverse_complement(normalized)
    return min(normalized, reverse)


def iter_canonical_kmers(
    sequence: str, k: int = 31
) -> Iterator[tuple[int, str]]:
    normalized = normalize_dna(sequence)
    for start in range(max(0, len(normalized) - k + 1)):
        word = normalized[start : start + k]
        if set(word) <= {"A", "C", "G", "T"}:
            yield start, canonical_kmer(word)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()
