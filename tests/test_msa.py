import pytest

from pipeline.core import run_msa


@pytest.mark.asyncio
async def test_msa_aligns_identical_sequences_at_100_percent():
    seqs = [
        {"accession": "SEQ1", "sequence": "ACGTACGTACGT"},
        {"accession": "SEQ2", "sequence": "ACGTACGTACGT"},
    ]
    result = await run_msa(seqs)
    assert "error" not in result
    assert result["num_sequences"] == 2
    identities = [a["identity_to_ref"] for a in result["aligned"]]
    assert identities[0] == 100.0
    assert identities[1] == 100.0


@pytest.mark.asyncio
async def test_msa_skipped_for_single_sequence():
    result = await run_msa([{"accession": "SEQ1", "sequence": "ACGT"}])
    assert result["skipped"] is True
