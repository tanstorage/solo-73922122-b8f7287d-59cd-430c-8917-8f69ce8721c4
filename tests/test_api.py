"""Tests for the transmission solver and the JSON API."""

from __future__ import annotations

import itertools

import pytest
from fastapi.testclient import TestClient

from app.main import app
from app.solver import Site, analyze

client = TestClient(app)


def make_site(position: int, father: str, mother: str, child: str) -> Site:
    return Site(
        position=position,
        father_alleles=(int(father[0]), int(father[2])),
        mother_alleles=(int(mother[0]), int(mother[2])),
        child_alleles=tuple(sorted((int(child[0]), int(child[2])))),
    )


def brute_force(sites: list[Site]):
    """Enumerate every transmission plan; return (min cost, optimal plans).

    An optimal plan is a tuple of (father_source, mother_source) per site.
    Returns (None, []) when any site has no legal transmission.
    """

    choices = []
    for site in sites:
        legal = []
        for f in (0, 1):
            for m in (0, 1):
                if sorted((site.father_alleles[f], site.mother_alleles[m])) == list(
                    site.child_alleles
                ):
                    legal.append((f, m))
        if not legal:
            return None, []
        choices.append(legal)

    best = None
    plans = []
    for plan in itertools.product(*choices):
        cost = 0
        for (f0, m0), (f1, m1) in zip(plan, plan[1:]):
            cost += int(f0 != f1) + int(m0 != m1)
        if best is None or cost < best:
            best = cost
            plans = [plan]
        elif cost == best:
            plans.append(plan)
    return best, plans


def test_simple_switch_free_segment():
    sites = [
        make_site(100, "0|1", "1|0", "0/0"),
        make_site(200, "0|1", "1|0", "0/0"),
    ]
    result = analyze(sites, 0)
    assert result.feasible and result.min_switches == 0
    assert result.status == "ok"
    # Only (father=0, mother=1) is legal at each site.
    assert [s.father.source for s in result.sites] == [0, 0]
    assert [s.mother.source for s in result.sites] == [1, 1]
    assert result.ambiguous_sites == ()


def test_both_parents_switch_required():
    # Site 1 forces (father=0, mother=0); site 2 forces (father=1, mother=1).
    # Both parental sources change, so the cost is 2.
    sites = [
        make_site(100, "0|1", "0|1", "0/0"),
        make_site(200, "0|1", "0|1", "1/1"),
    ]
    result = analyze(sites, 2)
    assert result.min_switches == 2
    assert result.status == "ok"


def test_limit_exceeded_flags_whole_segment():
    sites = [
        make_site(100, "0|1", "0|1", "0/0"),
        make_site(200, "0|1", "0|1", "1/1"),
    ]
    result = analyze(sites, 1)
    assert result.min_switches == 2
    assert result.over_limit is True
    assert result.status == "needs_review"


def test_conflict_site_reported():
    sites = [
        make_site(100, "0|1", "1|0", "0/0"),
        make_site(200, "0|0", "0|0", "0/1"),
    ]
    result = analyze(sites, 10)
    assert result.feasible is False
    assert result.conflict_sites == (200,)
    assert result.min_switches is None
    assert result.status == "conflict"
    assert result.sites[0].status == "unavailable"
    assert result.sites[1].status == "conflict"


def test_ambiguity_when_optimal_plans_disagree():
    # Site 100 forces state (father=1, mother=1).  Site 200 allows (0,1) and
    # (1,1); site 300 forces (0,1).  Switching at site 200 or at site 300 both
    # cost one paternal change, so the father's source at site 200 is ambiguous
    # (he is homozygous there, so the allele is still known) while the mother's
    # source stays determined.
    sites = [
        make_site(100, "0|1", "0|1", "1/1"),
        make_site(200, "0|0", "0|1", "0/1"),
        make_site(300, "0|1", "1|0", "0/0"),
    ]
    result = analyze(sites, 5)
    assert result.min_switches == 1
    assert result.sites[0].status == "determined"
    assert result.sites[1].status == "ambiguous"
    assert result.sites[1].father.source is None
    assert result.sites[1].father.allele == 0
    assert [source for source, _ in result.sites[1].father.options] == [0, 1]
    assert result.sites[1].mother.status == "determined"
    assert result.sites[1].mother.source == 1
    assert result.sites[2].status == "determined"
    assert result.ambiguous_sites == (200,)


def test_ambiguity_is_global_not_per_site():
    # Site 1 allows states 0 and 2; site 2 allows states 0 and 1.  Per site the
    # source is ambiguous, but only state 0 at both sites costs 0 switches, so
    # the global optimum determines every source.
    sites = [
        make_site(100, "0|0", "0|1", "0/0"),  # states (0,0) and (1,0)
        make_site(200, "0|1", "0|0", "0/0"),  # states (0,0) and (0,1)
    ]
    result = analyze(sites, 5)
    assert result.min_switches == 0
    assert [s.status for s in result.sites] == ["determined", "determined"]
    assert [s.father.source for s in result.sites] == [0, 0]
    assert [s.mother.source for s in result.sites] == [0, 0]


def test_example_assignment_is_optimal_and_consistent():
    sites = [
        make_site(100, "0|1", "0|1", "0/0"),
        make_site(200, "0|1", "0|1", "0/1"),
        make_site(300, "1|0", "0|1", "0/0"),
    ]
    result = analyze(sites, 5)
    example = result.example_assignment
    assert example is not None
    assert example.total_switches == result.min_switches
    assert example.total_switches == example.father_switches + example.mother_switches
    for item, site_result in zip(example.sites, result.sites):
        assert item.father_allele == site_result.site.father_alleles[item.father_source]
        assert item.mother_allele == site_result.site.mother_alleles[item.mother_source]
        assert sorted((item.father_allele, item.mother_allele)) == list(
            site_result.site.child_alleles
        )
        if site_result.father.status == "determined":
            assert item.father_source == site_result.father.source
        if site_result.mother.status == "determined":
            assert item.mother_source == site_result.mother.source


def test_dp_matches_brute_force_on_random_segments():
    import random

    rng = random.Random(20261115)
    for _ in range(400):
        n = rng.randint(2, 7)
        sites = []
        for i in range(n):
            father = f"{rng.randint(0, 1)}|{rng.randint(0, 1)}"
            mother = f"{rng.randint(0, 1)}|{rng.randint(0, 1)}"
            child = f"{rng.randint(0, 1)}/{rng.randint(0, 1)}"
            sites.append(make_site((i + 1) * 100, father, mother, child))

        best, plans = brute_force(sites)
        result = analyze(sites, 40)

        if best is None:
            assert result.feasible is False
            assert result.conflict_sites
            continue

        assert result.feasible is True
        assert result.min_switches == best
        for i, site_result in enumerate(result.sites):
            expected_father = {plan[i][0] for plan in plans}
            expected_mother = {plan[i][1] for plan in plans}
            actual_father = {source for source, _ in site_result.father.options}
            actual_mother = {source for source, _ in site_result.mother.options}
            assert actual_father == expected_father
            assert actual_mother == expected_mother
            assert (site_result.father.source is not None) == (
                len(expected_father) == 1
            )
            assert (site_result.mother.source is not None) == (
                len(expected_mother) == 1
            )


def test_site_count_bounds():
    one = {"chromosome": "chr1", "max_switches": 0, "sites": [
        {"position": 1, "father": "0|1", "mother": "0|1", "child": "0/1"}
    ]}
    assert client.post("/api/v1/analyze", json=one).status_code == 422

    forty_one = {
        "chromosome": "chr1",
        "max_switches": 40,
        "sites": [
            {"position": i, "father": "0|1", "mother": "0|1", "child": "0/1"}
            for i in range(41)
        ],
    }
    assert client.post("/api/v1/analyze", json=forty_one).status_code == 422


def test_api_happy_path_and_docs():
    assert client.get("/docs").status_code == 200
    assert client.get("/openapi.json").status_code == 200
    assert client.get("/healthz").json() == {"status": "ok"}

    payload = {
        "chromosome": "chr7",
        "max_switches": 1,
        "sites": [
            {"position": 1000, "father": "0|1", "mother": "1|0", "child": "0/0"},
            {"position": 2000, "father": "0|1", "mother": "1|0", "child": "0/0"},
        ],
    }
    response = client.post("/api/v1/analyze", json=payload)
    assert response.status_code == 200
    body = response.json()
    assert body["feasible"] is True
    assert body["min_switches"] == 0
    assert body["status"] == "ok"
    assert body["determined_positions"] == [1000, 2000]
    assert body["example_assignment"]["total_switches"] == 0
    assert body["site_count"] == 2


def test_api_rejects_bad_genotypes_and_ordering():
    base = {
        "max_switches": 1,
        "sites": [
            {"position": 2000, "father": "0|1", "mother": "1|0", "child": "0/1"},
            {"position": 1000, "father": "0|1", "mother": "1|0", "child": "0/1"},
        ],
    }
    assert client.post("/api/v1/analyze", json=base).status_code == 422

    bad = {
        "max_switches": 1,
        "sites": [
            {"position": 1, "father": "0/1", "mother": "1|0", "child": "0/1"},
            {"position": 2, "father": "0|1", "mother": "1|0", "child": "0/1"},
        ],
    }
    assert client.post("/api/v1/analyze", json=bad).status_code == 422

    out_of_range = {
        "max_switches": 79,
        "sites": [
            {"position": 1, "father": "0|1", "mother": "1|0", "child": "0/1"},
            {"position": 2, "father": "0|1", "mother": "1|0", "child": "0/1"},
        ],
    }
    assert client.post("/api/v1/analyze", json=out_of_range).status_code == 422


def test_api_conflict_and_over_limit_responses():
    conflict = {
        "max_switches": 5,
        "sites": [
            {"position": 10, "father": "0|0", "mother": "0|0", "child": "0/1"},
            {"position": 20, "father": "0|1", "mother": "1|0", "child": "0/1"},
        ],
    }
    body = client.post("/api/v1/analyze", json=conflict).json()
    assert body["feasible"] is False
    assert body["status"] == "conflict"
    assert body["conflict_positions"] == [10]
    assert body["min_switches"] is None
    assert body["example_assignment"] is None

    over = {
        "max_switches": 1,
        "sites": [
            {"position": 10, "father": "0|1", "mother": "0|1", "child": "0/0"},
            {"position": 20, "father": "0|1", "mother": "0|1", "child": "1/1"},
        ],
    }
    body = client.post("/api/v1/analyze", json=over).json()
    assert body["feasible"] is True
    assert body["status"] == "needs_review"
    assert body["over_limit"] is True
    assert body["min_switches"] == 2
    assert body["limit_excess"] == 1


def test_api_marks_homozygous_parent_ambiguity_by_source_index():
    payload = {
        "max_switches": 3,
        "sites": [
            {"position": 10, "father": "0|0", "mother": "0|0", "child": "0/0"},
            {"position": 20, "father": "0|0", "mother": "0|0", "child": "0/0"},
        ],
    }
    body = client.post("/api/v1/analyze", json=payload).json()
    assert body["status"] == "ok"
    # Both haplotype indices carry allele 0, so every plan is optimal and the
    # index cannot be determined even though the allele is unambiguous.
    assert body["ambiguous_positions"] == [10, 20]
    assert body["sites"][0]["father"]["source"] is None
    assert body["sites"][0]["father"]["options"] == [
        {"source": 0, "allele": 0},
        {"source": 1, "allele": 0},
    ]
