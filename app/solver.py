"""Core algorithm for pedigree haplotype-transmission analysis.

Each site has a phased father genotype, a phased mother genotype and an
unphased child genotype.  At every site one haplotype is transmitted by each
parent.  The haplotype index (0 or 1) is called the *source*.  The cost of a
transmission plan is the number of source changes between consecutive sites,
summed over both parents.

The solver uses dynamic programming over the four possible (father, mother)
source pairs and a forward/backward pass to decide which source values occur in
at least one minimum-cost plan.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

# State encoding: state = 2 * father_source + mother_source.
# Allowed states are (father, mother) in the order (0,0), (0,1), (1,0), (1,1).
STATE_COUNT = 4
INF = 10**9

# TRANSITION_COST[a][b] is the number of parents whose source changes when the
# transmission pair changes from state ``a`` to state ``b``.
TRANSITION_COST: tuple[tuple[int, ...], ...] = tuple(
    tuple(
        int((a >> 1) != (b >> 1)) + int((a & 1) != (b & 1))
        for b in range(STATE_COUNT)
    )
    for a in range(STATE_COUNT)
)


@dataclass(frozen=True)
class Site:
    """A single biallelic site with phased parental genotypes."""

    position: int
    father_alleles: tuple[int, int]
    mother_alleles: tuple[int, int]
    child_alleles: tuple[int, int]  # sorted, unphased


@dataclass(frozen=True)
class ParentResult:
    """Transmitted-haplotype inference for one parent at one site."""

    status: str  # determined | ambiguous | conflict | unavailable
    source: int | None
    allele: int | None
    options: tuple[tuple[int, int], ...]  # (source, allele), sorted by source


@dataclass(frozen=True)
class SiteResult:
    index: int
    site: Site
    status: str  # determined | ambiguous | conflict | unavailable
    father: ParentResult
    mother: ParentResult


@dataclass(frozen=True)
class ExampleSite:
    index: int
    position: int
    father_source: int
    mother_source: int
    father_allele: int
    mother_allele: int


@dataclass(frozen=True)
class ExampleAssignment:
    """One minimum-cost transmission plan (ties broken deterministically)."""

    father_switches: int
    mother_switches: int
    total_switches: int
    sites: tuple[ExampleSite, ...]


@dataclass(frozen=True)
class Analysis:
    feasible: bool
    status: str  # ok | needs_review | conflict
    min_switches: int | None
    over_limit: bool
    conflict_sites: tuple[int, ...]  # positions
    ambiguous_sites: tuple[int, ...]  # positions
    sites: tuple[SiteResult, ...]
    example_assignment: ExampleAssignment | None


def legal_states(
    father_alleles: Sequence[int],
    mother_alleles: Sequence[int],
    child_alleles: Sequence[int],
) -> list[int]:
    """Return all (father_source, mother_source) pairs explaining the child.

    A parent contributes the allele carried by the haplotype it transmits.  The
    child genotype is unphased, so its two alleles are compared as a multiset.
    """

    target = tuple(sorted(child_alleles))
    states: list[int] = []
    for father_source in (0, 1):
        for mother_source in (0, 1):
            transmitted = (
                father_alleles[father_source],
                mother_alleles[mother_source],
            )
            if tuple(sorted(transmitted)) == target:
                states.append(2 * father_source + mother_source)
    return states


def _parent_result(source_alleles: Sequence[tuple[int, int]]) -> ParentResult:
    """Build a parent result from the optimal (source, allele) options.

    A source is reported as determined only when every optimal plan uses the
    same haplotype index.  Homozygous parents carry two identical alleles but
    still have distinct haplotype indices, so they remain ambiguous unless the
    optimal plans agree on the index.
    """

    options = tuple(sorted(set(source_alleles)))
    sources = {source for source, _ in options}
    alleles = {allele for _, allele in options}
    allele = next(iter(alleles)) if len(alleles) == 1 else None
    if len(sources) == 1:
        source = next(iter(sources))
        return ParentResult(
            status="determined", source=source, allele=allele, options=options
        )
    return ParentResult(
        status="ambiguous", source=None, allele=allele, options=options
    )


def _unavailable_parent() -> ParentResult:
    return ParentResult(status="unavailable", source=None, allele=None, options=())


def _conflict_parent() -> ParentResult:
    return ParentResult(status="conflict", source=None, allele=None, options=())


def analyze(sites: Sequence[Site], max_switches: int) -> Analysis:
    """Solve the minimum-switch transmission problem.

    ``max_switches`` is only used to classify the segment; the minimum cost and
    per-site ambiguity are always computed from the transmission model.
    """

    if not sites:
        raise ValueError("at least one site is required")

    n = len(sites)
    legal = [
        legal_states(site.father_alleles, site.mother_alleles, site.child_alleles)
        for site in sites
    ]
    conflict_indices = [i for i, states in enumerate(legal) if not states]

    if conflict_indices:
        conflict_set = set(conflict_indices)
        site_results = []
        for i, site in enumerate(sites):
            if i in conflict_set:
                site_results.append(
                    SiteResult(
                        index=i,
                        site=site,
                        status="conflict",
                        father=_conflict_parent(),
                        mother=_conflict_parent(),
                    )
                )
            else:
                site_results.append(
                    SiteResult(
                        index=i,
                        site=site,
                        status="unavailable",
                        father=_unavailable_parent(),
                        mother=_unavailable_parent(),
                    )
                )
        return Analysis(
            feasible=False,
            status="conflict",
            min_switches=None,
            over_limit=False,
            conflict_sites=tuple(sites[i].position for i in conflict_indices),
            ambiguous_sites=(),
            sites=tuple(site_results),
            example_assignment=None,
        )

    # forward[i][state] = minimum cost for sites 0..i ending in ``state``.
    forward = [[INF] * STATE_COUNT for _ in range(n)]
    for state in legal[0]:
        forward[0][state] = 0

    for i in range(1, n):
        previous = legal[i - 1]
        for state in legal[i]:
            best = INF
            for prev_state in previous:
                candidate = (
                    forward[i - 1][prev_state] + TRANSITION_COST[prev_state][state]
                )
                if candidate < best:
                    best = candidate
            forward[i][state] = best

    # backward[i][state] = minimum cost for transitions after site i, starting
    # from ``state`` at site i.  The transition into site i is already counted
    # by the forward pass, so there is no double counting.
    backward = [[INF] * STATE_COUNT for _ in range(n)]
    for state in legal[n - 1]:
        backward[n - 1][state] = 0

    for i in range(n - 2, -1, -1):
        following = legal[i + 1]
        for state in legal[i]:
            best = INF
            for next_state in following:
                candidate = (
                    TRANSITION_COST[state][next_state] + backward[i + 1][next_state]
                )
                if candidate < best:
                    best = candidate
            backward[i][state] = best

    min_switches = min(forward[n - 1])
    over_limit = min_switches > max_switches

    # States at site i that belong to at least one globally optimal plan.
    optimal_states: list[list[int]] = []
    for i in range(n):
        optimal_states.append(
            [
                state
                for state in legal[i]
                if forward[i][state] + backward[i][state] == min_switches
            ]
        )

    site_results: list[SiteResult] = []
    ambiguous_positions: list[int] = []

    for i, site in enumerate(sites):
        states = optimal_states[i]
        father_options = [(state >> 1, site.father_alleles[state >> 1]) for state in states]
        mother_options = [(state & 1, site.mother_alleles[state & 1]) for state in states]

        father = _parent_result(father_options)
        mother = _parent_result(mother_options)

        site_status = (
            "determined"
            if father.status == "determined" and mother.status == "determined"
            else "ambiguous"
        )
        if site_status == "ambiguous":
            ambiguous_positions.append(site.position)

        site_results.append(
            SiteResult(
                index=i,
                site=site,
                status=site_status,
                father=father,
                mother=mother,
            )
        )

    # Reconstruct one optimal plan deterministically (smallest state index on
    # ties).  This is illustrative only: ambiguous sites may have other
    # equally optimal sources.
    state = min(legal[n - 1], key=lambda s: (forward[n - 1][s], s))
    path = [state]
    for i in range(n - 1, 0, -1):
        target = forward[i][state]
        state = min(
            previous
            for previous in legal[i - 1]
            if forward[i - 1][previous] + TRANSITION_COST[previous][state] == target
        )
        path.append(state)
    path.reverse()

    example_sites = []
    father_switches = 0
    mother_switches = 0
    for i, state in enumerate(path):
        if i > 0:
            previous = path[i - 1]
            father_switches += int((state >> 1) != (previous >> 1))
            mother_switches += int((state & 1) != (previous & 1))
        father_source = state >> 1
        mother_source = state & 1
        example_sites.append(
            ExampleSite(
                index=i,
                position=sites[i].position,
                father_source=father_source,
                mother_source=mother_source,
                father_allele=sites[i].father_alleles[father_source],
                mother_allele=sites[i].mother_alleles[mother_source],
            )
        )

    example = ExampleAssignment(
        father_switches=father_switches,
        mother_switches=mother_switches,
        total_switches=father_switches + mother_switches,
        sites=tuple(example_sites),
    )

    return Analysis(
        feasible=True,
        status="needs_review" if over_limit else "ok",
        min_switches=min_switches,
        over_limit=over_limit,
        conflict_sites=(),
        ambiguous_sites=tuple(ambiguous_positions),
        sites=tuple(site_results),
        example_assignment=example,
    )
