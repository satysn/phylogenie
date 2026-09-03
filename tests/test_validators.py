import pytest

from pipeline.validators import (
    OPTIONAL_TOOLS,
    ValidationError,
    validate_accession,
    validate_accessions,
    validate_organism_model,
    validate_tools,
)


def test_valid_accession_formats():
    assert validate_accession("nm_001301717") == "NM_001301717"
    assert validate_accession(" NM_000546.2 ") == "NM_000546.2"
    assert validate_accession("AB012345") == "AB012345"


def test_invalid_accession_rejected():
    with pytest.raises(ValidationError):
        validate_accession("")
    with pytest.raises(ValidationError):
        validate_accession("'; DROP TABLE jobs; --")
    with pytest.raises(ValidationError):
        validate_accession("not an accession!!")


def test_validate_accessions_dedupes_and_caps():
    result = validate_accessions(["NM_001301717", "nm_001301717", "NM_000546"], max_sequences=10)
    assert result == ["NM_001301717", "NM_000546"]

    with pytest.raises(ValidationError):
        validate_accessions(["NM_001301717", "NM_000546", "NM_004333"], max_sequences=2)

    with pytest.raises(ValidationError):
        validate_accessions([], max_sequences=10)


def test_validate_organism_model_falls_back_to_human():
    assert validate_organism_model("Mouse") == "mouse"
    assert validate_organism_model("not-a-real-organism") == "human"
    assert validate_organism_model("") == "human"


def test_validate_tools_defaults_all_enabled():
    result = validate_tools(None)
    assert set(result) == set(OPTIONAL_TOOLS)
    assert all(result.values())


def test_validate_tools_respects_disabled_selection():
    result = validate_tools({"blastn": False, "blastp": False})
    assert result["blastn"] is False
    assert result["blastp"] is False
    assert result["augustus"] is True  # untouched keys default to enabled


def test_validate_tools_drops_unknown_keys():
    result = validate_tools({"blastn": False, "not_a_real_tool": True})
    assert "not_a_real_tool" not in result
    assert set(result) == set(OPTIONAL_TOOLS)
