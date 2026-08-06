import gzip
from pathlib import Path

import pytest

from rrna_bait.core import (
    canonical_kmer,
    iter_canonical_kmers,
    normalize_dna,
    open_text,
    reverse_complement,
    sha256_file,
)


def test_normalize_dna_uppercases_and_converts_uracil():
    assert normalize_dna("acguN") == "ACGTN"


def test_reverse_complement():
    assert reverse_complement("AACCGT") == "ACGGTT"


def test_canonical_kmer_uses_lexicographically_smaller_strand():
    forward = "AACCGT"
    assert canonical_kmer(forward) == min(forward, reverse_complement(forward))


def test_canonical_kmer_rejects_non_acgt_input_after_normalization():
    with pytest.raises(ValueError, match="non-ACGT"):
        canonical_kmer("acguN")


def test_iter_canonical_kmers_rejects_ambiguous_windows():
    sequence = "A" * 31 + "N" + "C" * 31
    observed = list(iter_canonical_kmers(sequence, 31))
    assert observed == [(0, "A" * 31), (32, "C" * 31)]


def test_sha256_file(tmp_path: Path):
    path = tmp_path / "value.txt"
    path.write_text("cyclospora\n")
    assert sha256_file(path) == (
        "3b090da6316bd6e5d82635b3c77a9ef4"
        "9a1419ff8e6705419666354ecbf93c96"
    )


def test_open_text_reads_plain_text(tmp_path: Path):
    path = tmp_path / "value.txt"
    path.write_text("cyclospora\\n")

    with open_text(path) as handle:
        assert handle.read() == "cyclospora\\n"


def test_open_text_reads_gzip_text(tmp_path: Path):
    path = tmp_path / "value.txt.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("cyclospora\\n")

    with open_text(path) as handle:
        assert handle.read() == "cyclospora\\n"
