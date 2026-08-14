import csv
import gzip
from pathlib import Path
import tomllib

from Bio import SeqIO
import pytest

from rrna_bait import sources
from rrna_bait.core import sha256_file
from rrna_bait.sources import (
    accession_root,
    atomic_download,
    filter_rfam,
    filter_silva,
    main,
    reconstruct_gzip_concat,
    verify_expected_identity,
    write_provenance_row,
)


def test_accession_root_removes_version_and_rfam_coordinates():
    assert accession_root("AF111183.1/10-100 description") == "AF111183"
    assert accession_root("AF111183.1.1795 lineage") == "AF111183"


def test_filter_silva_excludes_only_target_accessions(tmp_path: Path):
    source = tmp_path / "silva.fa.gz"
    output = tmp_path / "non_target.fa.gz"
    with gzip.open(source, "wt") as handle:
        handle.write(
            ">AF111183.1.1795 Eukaryota;Cyclospora cayetanensis\nACGT\n"
            ">AF111184.1.1798 Eukaryota;Cyclospora cercopitheci\nTGCA\n"
        )
    kept, removed = filter_silva(
        source, output, {"AF111183"}, "Cyclospora cayetanensis"
    )
    assert (kept, removed) == (1, 1)
    assert [r.id for r in SeqIO.parse(gzip.open(output, "rt"), "fasta")] == [
        "AF111184.1.1798"
    ]


def test_filter_rfam_excludes_target_accessions(tmp_path: Path):
    source = tmp_path / "rfam.fa.gz"
    output = tmp_path / "non_target.fa.gz"
    with gzip.open(source, "wt") as handle:
        handle.write(">XR_003297357.1/1-154\nACGU\n>OTHER.2/5-8\nUGCA\n")
    kept, removed = filter_rfam(source, output, {"XR_003297357"})
    assert (kept, removed) == (1, 1)


def test_filter_silva_only_excludes_exact_terminal_taxon_and_preserves_description(
    tmp_path: Path,
):
    source = tmp_path / "silva.fa.gz"
    output = tmp_path / "non_target.fa.gz"
    with gzip.open(source, "wt") as handle:
        handle.write(
            ">EXACT.1 lineage;Cyclospora cayetanensis\nacgu\n"
            ">TARGET.1 lineage;Microcyclospora cayetanensis\nacgu\n"
            ">OTHER.1 lineage;Cyclospora cayetanensis strain X\nacgu\n"
            ">RETAIN.1 original description;Microcyclospora cayetanensis\nugca\n"
        )

    assert filter_silva(source, output, set()) == (3, 1)
    records = list(SeqIO.parse(gzip.open(output, "rt"), "fasta"))
    assert [(record.description, str(record.seq)) for record in records] == [
        ("TARGET.1 lineage;Microcyclospora cayetanensis", "ACGT"),
        ("OTHER.1 lineage;Cyclospora cayetanensis strain X", "ACGT"),
        ("RETAIN.1 original description;Microcyclospora cayetanensis", "TGCA"),
    ]


def test_filter_outputs_deterministic_gzip_bytes(tmp_path: Path):
    source = tmp_path / "rfam.fa.gz"
    first_output = tmp_path / "first.fa.gz"
    second_output = tmp_path / "second.fa.gz"
    with gzip.open(source, "wt") as handle:
        handle.write(">RETAIN.1/1-4 description\nACGU\n")

    filter_rfam(source, first_output, set())
    filter_rfam(source, second_output, set())

    assert first_output.read_bytes() == second_output.read_bytes()


def test_atomic_download_replaces_destination_only_after_nonempty_download(
    tmp_path: Path,
):
    destination = tmp_path / "downloaded.fa.gz"
    destination.write_text("old data\n")
    source = tmp_path / "source.fa.gz"
    source.write_text("new data\n")

    atomic_download(source.as_uri(), destination)
    assert destination.read_text() == "new data\n"

    empty_source = tmp_path / "empty.fa.gz"
    empty_source.touch()
    with pytest.raises(ValueError, match="empty"):
        atomic_download(empty_source.as_uri(), destination)
    assert destination.read_text() == "new data\n"


def test_atomic_download_rejects_identity_mismatch_without_promotion(
    tmp_path: Path,
) -> None:
    destination = tmp_path / "cache.fa"
    destination.write_bytes(b"known-good")
    source = tmp_path / "remote.fa"
    source.write_bytes(b"changed")

    with pytest.raises(ValueError, match="expected identity mismatch"):
        atomic_download(
            source.as_uri(),
            destination,
            expected_bytes=10,
            expected_sha256=sha256_file(destination),
        )

    assert destination.read_bytes() == b"known-good"


def test_expected_identity_is_independent_of_observed_provenance(
    tmp_path: Path,
) -> None:
    cached = tmp_path / "cached.fa"
    cached.write_bytes(b"changed")

    with pytest.raises(ValueError, match="expected 10 bytes"):
        verify_expected_identity(
            cached, expected_bytes=10, expected_sha256="0" * 64
        )


def test_failed_gzip_concat_removes_temp_and_preserves_destination(
    tmp_path: Path,
) -> None:
    first = tmp_path / "first.fa.gz"
    broken = tmp_path / "broken.fa.gz"
    destination = tmp_path / "combined.fa"
    with gzip.open(first, "wt") as handle:
        handle.write(">first\nACGT\n")
    broken.write_bytes(b"not gzip")
    destination.write_bytes(b"known-good")

    with pytest.raises((OSError, ValueError)):
        reconstruct_gzip_concat(
            [first, broken],
            destination,
            expected_bytes=12,
            expected_sha256="0" * 64,
        )

    assert destination.read_bytes() == b"known-good"
    assert list(tmp_path.glob("combined.fa.tmp*")) == []


def test_production_build_uses_only_frozen_versioned_accessions() -> None:
    project = Path(__file__).parents[1]
    script = (project / "scripts/build_rrna_bait.sh").read_text()
    target_snapshot = project / "config/target_rrna_accessions.txt"
    other_snapshot = project / "config/other_cyclospora_accessions.txt"

    assert "esearch" not in script
    assert "target_rrna_accession_snapshot" in script
    assert "other_cyclospora_accession_snapshot" in script
    for snapshot, expected_count in (
        (target_snapshot, 329),
        (other_snapshot, 94),
    ):
        accessions = snapshot.read_text().splitlines()
        assert len(accessions) == expected_count
        assert len(accessions) == len(set(accessions))
        assert all(
            accession.rsplit(".", 1)[-1].isdigit() for accession in accessions
        )


def test_validate_task_rebuilds_and_searches_the_current_bait_set() -> None:
    project = Path(__file__).parents[1]
    workspace = tomllib.loads((project / "pixi.toml").read_text())
    task = workspace["tasks"]["validate"]
    script = (project / "scripts/validate_core_nt.sh").read_text()

    assert task["depends-on"] == ["build"]
    assert task["args"] == ["core_nt_database"]
    assert 'baits=kmers/cyclospora_cayetanensis_rrna_baits.fasta' in script
    assert 'bash scripts/build_core_nt_validated_index.sh "$blast_tsv"' in script
    assert "results/core_nt_bait_exact_match_blast.tsv" not in script

def test_write_provenance_replaces_source_row_and_sorts_by_name(tmp_path: Path):
    provenance = tmp_path / "provenance.tsv"
    silva = tmp_path / "silva.fa.gz"
    rfam = tmp_path / "rfam.fa.gz"
    silva.write_bytes(b"silva")
    rfam.write_bytes(b"rfam")

    write_provenance_row(
        provenance,
        "silva_ssu",
        "138.2",
        "https://silva",
        silva,
        project_root=tmp_path,
    )
    write_provenance_row(provenance, "rfam_5s", "15.1", "https://rfam", rfam)
    write_provenance_row(
        provenance,
        "silva_ssu",
        "138.3",
        "https://new-silva",
        silva,
        project_root=tmp_path,
    )

    with provenance.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["name"] for row in rows] == ["rfam_5s", "silva_ssu"]
    assert rows[1] == {
        "name": "silva_ssu",
        "release": "138.3",
        "url": "https://new-silva",
        "local_path": "silva.fa.gz",
        "retrieved_at_utc": rows[1]["retrieved_at_utc"],
        "bytes": str(silva.stat().st_size),
        "sha256": sha256_file(silva),
    }


def test_verify_provenance_accepts_matching_cached_file(tmp_path: Path):
    provenance = tmp_path / "provenance.tsv"
    cached = tmp_path / "cached.fasta"
    cached.write_text(">record\nACGT\n")
    write_provenance_row(
        provenance, "target_queries", "NCBI-live", "edirect:accessions:A1", cached
    )
    original_provenance = provenance.read_text()

    assert sources.verify_provenance_row(
        provenance,
        "target_queries",
        "NCBI-live",
        "edirect:accessions:A1",
        cached,
    )
    assert provenance.read_text() == original_provenance


def test_write_provenance_replaces_superseded_name_for_same_cache(
    tmp_path: Path,
) -> None:
    provenance = tmp_path / "provenance.tsv"
    cached = tmp_path / "cached.fasta"
    cached.write_text(">record\nACGT\n")
    write_provenance_row(
        provenance, "old_transport", "live", "edirect:accessions:A1", cached
    )
    write_provenance_row(
        provenance, "new_transport", "live", "https://eutils/A1", cached
    )

    with provenance.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["name"] for row in rows] == ["new_transport"]


def test_verify_provenance_rejects_tampered_cache_and_missing_row(tmp_path: Path):
    provenance = tmp_path / "provenance.tsv"
    cached = tmp_path / "cached.fasta"
    cached.write_text(">record\nACGT\n")
    write_provenance_row(
        provenance, "target_queries", "NCBI-live", "edirect:accessions:A1", cached
    )
    cached.write_text(">record\nTGCA\n")

    assert not sources.verify_provenance_row(
        provenance,
        "target_queries",
        "NCBI-live",
        "edirect:accessions:A1",
        cached,
    )
    assert not sources.verify_provenance_row(
        provenance,
        "other_cyclospora",
        "NCBI-live",
        "edirect:query:other",
        cached,
    )


def test_cli_subcommands_filter_sources_and_record_provenance(tmp_path: Path):
    silva_source = tmp_path / "silva.fa.gz"
    silva_output = tmp_path / "silva-output.fa.gz"
    rfam_source = tmp_path / "rfam.fa.gz"
    rfam_output = tmp_path / "rfam-output.fa.gz"
    provenance = tmp_path / "provenance.tsv"
    with gzip.open(silva_source, "wt") as handle:
        handle.write(">TARGET.1 lineage;Cyclospora cayetanensis\nACGT\n")
    with gzip.open(rfam_source, "wt") as handle:
        handle.write(">TARGET.1/1-4\nACGU\n")

    main(["filter-silva", str(silva_source), str(silva_output)])
    main(
        [
            "filter-rfam",
            str(rfam_source),
            str(rfam_output),
            "--excluded-accession",
            "TARGET",
        ]
    )
    main(
        [
            "record-provenance",
            str(provenance),
            "--name",
            "rfam_5s",
            "--release",
            "15.1",
            "--url",
            "https://rfam",
            "--local-path",
            str(rfam_source),
        ]
    )
    main(
        [
            "verify-provenance",
            str(provenance),
            "--name",
            "rfam_5s",
            "--release",
            "15.1",
            "--url",
            "https://rfam",
            "--local-path",
            str(rfam_source),
        ]
    )

    assert list(SeqIO.parse(gzip.open(silva_output, "rt"), "fasta")) == []
    assert list(SeqIO.parse(gzip.open(rfam_output, "rt"), "fasta")) == []
    assert provenance.read_text().splitlines()[1].startswith("rfam_5s\t15.1\t")

    rfam_source.write_bytes(b"tampered")
    with pytest.raises(SystemExit) as error:
        main(
            [
                "verify-provenance",
                str(provenance),
                "--name",
                "rfam_5s",
                "--release",
                "15.1",
                "--url",
                "https://rfam",
                "--local-path",
                str(rfam_source),
            ]
        )
    assert error.value.code == 1
