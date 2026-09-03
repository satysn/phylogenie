import pytest

from pipeline.core import run_protein_analysis

# A real short protein sequence (human insulin B chain, no network needed)
INSULIN_B = "FVNQHLCGSHLVEALYLVCGERGFFYTPKT"


@pytest.mark.asyncio
async def test_protein_analysis_returns_physicochemical_properties():
    result = await run_protein_analysis(INSULIN_B, gene_id="insulin_b")
    assert "error" not in result
    assert result["gene_id"] == "insulin_b"
    assert result["protein_length"] == len(INSULIN_B)
    assert result["molecular_weight"] > 0
    assert 0 < result["isoelectric_point"] < 14
    ss = result["secondary_structure"]
    total = ss["helix"] + ss["turn"] + ss["sheet"] + ss["coil"]
    assert 95 <= total <= 105  # allow small rounding drift


@pytest.mark.asyncio
async def test_protein_analysis_rejects_too_short_sequence():
    result = await run_protein_analysis("MK", gene_id="tiny")
    assert "error" in result
