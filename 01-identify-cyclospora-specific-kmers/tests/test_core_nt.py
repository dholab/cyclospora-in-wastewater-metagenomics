from pathlib import Path

from rrna_bait.core_nt import BlastHit, classify_hits, write_decisions


def test_exact_non_target_hit_rejects_bait():
    hit = BlastHit.from_tsv(
        "bait_1\tMK946294.1\t342097\tuncultured soil eukaryote"
        "\t100.000\t31\t0\t0\t1\t31\t54\t24\tsoil SSU"
    )
    assert classify_hits({"bait_1"}, [hit])["bait_1"].status == (
        "REJECT_CORE_NT_EXACT_NON_TARGET"
    )


def test_target_exact_is_allowed_and_30_of_31_is_informational():
    target = BlastHit.from_tsv(
        "bait_1\tX.1\t88456\tCyclospora cayetanensis\t100.000\t31\t0\t0"
        "\t1\t31\t1\t31\ttarget"
    )
    near = BlastHit.from_tsv(
        "bait_1\tY.1\t1\tother\t96.774\t31\t1\t0\t1\t31\t1\t31\tnear"
    )
    decision = classify_hits({"bait_1"}, [target, near])["bait_1"]
    assert (decision.status, decision.near_match_count) == ("PASS_CORE_NT", 1)


def test_exact_hit_without_a_taxid_rejects_bait():
    hit = BlastHit.from_tsv(
        "bait_1\tX.1\tN/A\tunknown\t100.000\t31\t0\t0\t1\t31\t1\t31\tunknown"
    )

    assert classify_hits({"bait_1"}, [hit])["bait_1"].status == (
        "REJECT_CORE_NT_EXACT_NON_TARGET"
    )


def test_write_decisions_writes_sorted_tsv(tmp_path: Path):
    decisions = classify_hits({"bait_b", "bait_a"}, [])
    path = tmp_path / "decisions.tsv"

    write_decisions(decisions, path)

    assert path.read_text().splitlines() == [
        "bait_id\tstatus\texact_target_count\texact_non_target_count\tnear_match_count",
        "bait_a\tPASS_CORE_NT\t0\t0\t0",
        "bait_b\tPASS_CORE_NT\t0\t0\t0",
    ]
