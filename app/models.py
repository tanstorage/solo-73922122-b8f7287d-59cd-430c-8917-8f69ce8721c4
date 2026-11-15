"""Pydantic request/response models for the pedigree haplotype audit API."""

from __future__ import annotations

import re
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

_PHASED_RE = re.compile(r"^[01]\|[01]$")
_UNPHASED_RE = re.compile(r"^[01]/[01]$")


class SiteInput(BaseModel):
    """One biallelic site, ordered by increasing position in the request."""

    model_config = ConfigDict(
        json_schema_extra={
            "examples": [
                {
                    "position": 1000,
                    "father": "0|1",
                    "mother": "1|0",
                    "child": "0/1",
                }
            ]
        }
    )

    position: int = Field(
        ge=0,
        description="Chromosomal position; must be strictly increasing across sites.",
    )
    father: str = Field(description="Phased paternal genotype, e.g. '0|1'.")
    mother: str = Field(description="Phased maternal genotype, e.g. '1|0'.")
    child: str = Field(description="Unphased child genotype, e.g. '0/1'.")

    @field_validator("father", "mother")
    @classmethod
    def _validate_phased(cls, value: str) -> str:
        value = value.strip()
        if not _PHASED_RE.fullmatch(value):
            raise ValueError(
                "phased genotype must be exactly two alleles separated by '|', "
                "e.g. '0|1' or '1|0'"
            )
        return value

    @field_validator("child")
    @classmethod
    def _validate_unphased(cls, value: str) -> str:
        value = value.strip()
        if not _UNPHASED_RE.fullmatch(value):
            raise ValueError(
                "child genotype must be exactly two alleles separated by '/', "
                "e.g. '0/1' or '1/1'"
            )
        return value

    @property
    def father_alleles(self) -> tuple[int, int]:
        return (int(self.father[0]), int(self.father[2]))

    @property
    def mother_alleles(self) -> tuple[int, int]:
        return (int(self.mother[0]), int(self.mother[2]))

    @property
    def child_alleles(self) -> tuple[int, int]:
        return tuple(sorted((int(self.child[0]), int(self.child[2]))))  # type: ignore[return-value]


class AnalysisRequest(BaseModel):
    """Request body for POST /api/v1/analyze."""

    model_config = ConfigDict(
        json_schema_extra={
            "examples": [
                {
                    "chromosome": "chr1",
                    "max_switches": 1,
                    "sites": [
                        {
                            "position": 1000,
                            "father": "0|1",
                            "mother": "1|0",
                            "child": "0/1",
                        },
                        {
                            "position": 2000,
                            "father": "0|1",
                            "mother": "0|1",
                            "child": "1/1",
                        },
                    ],
                }
            ]
        }
    )

    chromosome: str | None = Field(
        default=None,
        max_length=128,
        description="Optional chromosome label shared by all sites.",
    )
    max_switches: int = Field(
        ge=0,
        le=78,
        description=(
            "Maximum allowed total of paternal and maternal haplotype-index "
            "changes (0-78). Exceeding it marks the whole segment for review."
        ),
    )
    sites: list[SiteInput] = Field(
        min_length=2,
        max_length=40,
        description="2 to 40 biallelic sites, sorted by strictly increasing position.",
    )

    @model_validator(mode="after")
    def _validate_order(self) -> "AnalysisRequest":
        positions = [site.position for site in self.sites]
        for previous, current in zip(positions, positions[1:]):
            if current <= previous:
                raise ValueError(
                    "sites must be ordered by strictly increasing position; "
                    f"found {previous} before {current}"
                )
        return self


class SourceAllele(BaseModel):
    source: int = Field(description="Haplotype index transmitted by this parent.")
    allele: int = Field(description="Allele carried by that haplotype at this site.")


class ParentInference(BaseModel):
    status: Literal["determined", "ambiguous", "conflict", "unavailable"]
    source: int | None = Field(
        default=None,
        description="Haplotype index when every optimal plan agrees, else null.",
    )
    allele: int | None = Field(
        default=None,
        description="Allele carried by the determined haplotype, else null.",
    )
    options: list[SourceAllele] = Field(
        description="All haplotype indices used by at least one optimal plan."
    )


class SiteResult(BaseModel):
    index: int
    position: int
    status: Literal["determined", "ambiguous", "conflict", "unavailable"]
    father: ParentInference
    mother: ParentInference


class ExampleSite(BaseModel):
    index: int
    position: int
    father_source: int
    mother_source: int
    father_allele: int
    mother_allele: int


class ExampleAssignment(BaseModel):
    father_switches: int
    mother_switches: int
    total_switches: int
    sites: list[ExampleSite]


class AnalysisResponse(BaseModel):
    chromosome: str | None
    site_count: int
    max_switches: int
    feasible: bool
    status: Literal["ok", "needs_review", "conflict"]
    min_switches: int | None = Field(
        default=None,
        description="Minimum total haplotype-index changes, null when infeasible.",
    )
    over_limit: bool
    limit_excess: int | None = Field(
        default=None,
        description="min_switches - max_switches when positive, else null.",
    )
    conflict_positions: list[int] = Field(
        description="Positions whose child genotype cannot be explained by the parents."
    )
    ambiguous_positions: list[int] = Field(
        description="Positions whose parental source is not unique among optimal plans."
    )
    determined_positions: list[int]
    sites: list[SiteResult]
    example_assignment: ExampleAssignment | None = Field(
        default=None,
        description=(
            "One minimum-switch plan, tie-broken deterministically. Sites marked "
            "ambiguous may also be explained by other equally optimal plans."
        ),
    )
    message: str
