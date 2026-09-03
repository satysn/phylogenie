from pipeline.core import find_orfs


def test_find_orfs_detects_start_and_stop():
    # ATG (M) + 9x GCT (A) + TAA (stop) -> 10 aa protein (M + 9 Ala)
    seq = "ATG" + "GCT" * 9 + "TAA"
    orfs = find_orfs(seq, min_len=9)
    assert orfs, "expected at least one ORF"
    best = orfs[0]
    assert best["protein"].startswith("M")
    assert best["length_aa"] == 10
    assert best["has_stop"] is True


def test_find_orfs_respects_min_length():
    seq = "ATG" + "GCT" * 2 + "TAA"  # 3 aa ORF, below typical min_len
    orfs = find_orfs(seq, min_len=100)
    assert orfs == []


def test_find_orfs_handles_empty_sequence():
    assert find_orfs("", min_len=100) == []
