"""FastAPI service for auditing parent-of-origin transmission before phasing."""

from __future__ import annotations

from fastapi import FastAPI
from fastapi.responses import RedirectResponse

from .models import (
    AnalysisRequest,
    AnalysisResponse,
    ExampleAssignment,
    ExampleSite,
    ParentInference,
    SiteResult,
    SourceAllele,
)
from .solver import Analysis, Site, analyze

app = FastAPI(
    title="Pedigree Haplotype Transmission Audit",
    version="1.0.0",
    description=(
        "Checks, before SNP phasing, whether a child's unphased genotypes can be "
        "explained by one transmitted haplotype from each phased parent, and "
        "computes the minimum number of haplotype-index changes across a "
        "chromosome segment.\n\n"
        "* Every site must transmit exactly one allele from each parent.\n"
        "* Cost = paternal source changes + maternal source changes.\n"
        "* A site is `determined` only when all minimum-cost plans agree on the "
        "source index for both parents; otherwise it is `ambiguous`.\n"
        "* Sites that no transmission can explain are reported as `conflict`.\n"
        "* If the minimum exceeds `max_switches`, the whole segment is "
        "`needs_review`.\n\n"
        "Open **POST /api/v1/analyze** below and press *Try it out* to submit a "
        "request and inspect the JSON result."
    ),
    docs_url="/docs",
    redoc_url="/redoc",
)


def _to_site_result(result) -> SiteResult:
    def parent(model) -> ParentInference:
        return ParentInference(
            status=model.status,
            source=model.source,
            allele=model.allele,
            options=[SourceAllele(source=s, allele=a) for s, a in model.options],
        )

    return SiteResult(
        index=result.index,
        position=result.site.position,
        status=result.status,
        father=parent(result.father),
        mother=parent(result.mother),
    )


def _message(request: AnalysisRequest, analysis: Analysis) -> str:
    if not analysis.feasible:
        positions = ", ".join(str(p) for p in analysis.conflict_sites)
        return (
            "No valid transmission exists: the child genotype cannot be formed "
            f"from one allele per parent at position(s) {positions}."
        )
    if analysis.over_limit:
        return (
            f"Minimum switches {analysis.min_switches} exceed the limit of "
            f"{request.max_switches}; the whole segment is flagged for review."
        )
    if analysis.ambiguous_sites:
        positions = ", ".join(str(p) for p in analysis.ambiguous_sites)
        return (
            f"Minimum switches {analysis.min_switches} within the limit, but the "
            f"parental source is not unique at position(s) {positions}."
        )
    return (
        f"Minimum switches {analysis.min_switches} within the limit; every site "
        "has a unique parental source."
    )


@app.get("/", include_in_schema=False)
def root() -> RedirectResponse:
    return RedirectResponse(url="/docs")


@app.get("/healthz", tags=["ops"], summary="Liveness probe")
def healthz() -> dict[str, str]:
    return {"status": "ok"}


@app.post(
    "/api/v1/analyze",
    response_model=AnalysisResponse,
    summary="Audit parent-of-origin transmission across a segment",
    responses={
        422: {"description": "Malformed genotype, bad ordering, or limit out of range."}
    },
)
def analyze_endpoint(request: AnalysisRequest) -> AnalysisResponse:
    sites = [
        Site(
            position=site.position,
            father_alleles=site.father_alleles,
            mother_alleles=site.mother_alleles,
            child_alleles=site.child_alleles,
        )
        for site in request.sites
    ]
    analysis = analyze(sites, request.max_switches)

    determined = [
        result.site.position
        for result in analysis.sites
        if result.status == "determined"
    ]
    limit_excess = (
        analysis.min_switches - request.max_switches
        if analysis.min_switches is not None and analysis.over_limit
        else None
    )

    example = None
    if analysis.example_assignment is not None:
        example = ExampleAssignment(
            father_switches=analysis.example_assignment.father_switches,
            mother_switches=analysis.example_assignment.mother_switches,
            total_switches=analysis.example_assignment.total_switches,
            sites=[
                ExampleSite(
                    index=item.index,
                    position=item.position,
                    father_source=item.father_source,
                    mother_source=item.mother_source,
                    father_allele=item.father_allele,
                    mother_allele=item.mother_allele,
                )
                for item in analysis.example_assignment.sites
            ],
        )

    return AnalysisResponse(
        chromosome=request.chromosome,
        site_count=len(request.sites),
        max_switches=request.max_switches,
        feasible=analysis.feasible,
        status=analysis.status,
        min_switches=analysis.min_switches,
        over_limit=analysis.over_limit,
        limit_excess=limit_excess,
        conflict_positions=list(analysis.conflict_sites),
        ambiguous_positions=list(analysis.ambiguous_sites),
        determined_positions=determined,
        sites=[_to_site_result(result) for result in analysis.sites],
        example_assignment=example,
        message=_message(request, analysis),
    )
