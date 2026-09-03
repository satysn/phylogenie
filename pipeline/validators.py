"""
Input validation helpers.
"""

import re

# NCBI nucleotide accession formats, e.g. NM_001301717, NM_001301717.2, AB012345, U12345
_ACCESSION_RE = re.compile(r"^[A-Z]{1,4}_?\d{5,9}(\.\d+)?$")

VALID_ORGANISM_MODELS = {
    "human", "mouse", "arabidopsis", "drosophila", "zebrafish", "rice",
}


class ValidationError(ValueError):
    pass


def validate_accession(accession: str) -> str:
    acc = accession.strip().upper()
    if not acc:
        raise ValidationError("Accession cannot be empty")
    if len(acc) > 30:
        raise ValidationError(f"Accession '{acc[:30]}...' is too long")
    if not _ACCESSION_RE.match(acc):
        raise ValidationError(
            f"'{acc}' does not look like a valid NCBI nucleotide accession "
            "(expected e.g. NM_001301717 or AB012345)"
        )
    return acc


def validate_accessions(accessions: list[str], max_sequences: int) -> list[str]:
    cleaned = []
    seen = set()
    for a in accessions:
        acc = validate_accession(a)
        if acc not in seen:
            seen.add(acc)
            cleaned.append(acc)
    if not cleaned:
        raise ValidationError("At least one valid accession is required")
    if len(cleaned) > max_sequences:
        raise ValidationError(f"Maximum {max_sequences} sequences per job (got {len(cleaned)})")
    return cleaned


def validate_organism_model(name: str) -> str:
    name = (name or "human").strip().lower()
    return name if name in VALID_ORGANISM_MODELS else "human"
