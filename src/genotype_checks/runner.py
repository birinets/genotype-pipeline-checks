# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Ruslan Plakhuta
"""Load the cases and run a matcher over them."""

from __future__ import annotations

import copy
import json
from dataclasses import dataclass
from importlib import resources
from typing import Any, Callable

Matcher = Callable[[dict, dict], bool]


def load_cases() -> list[dict[str, Any]]:
    """Return the bundled cases."""
    path = resources.files("genotype_checks").joinpath("cases.json")
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)["cases"]


@dataclass(frozen=True)
class Result:
    case_id: str
    expect_match: bool
    got_match: bool
    why: str

    @property
    def passed(self) -> bool:
        return self.expect_match == self.got_match


def run_cases(matcher: Matcher, cases: list[dict[str, Any]] | None = None) -> list[Result]:
    """Call ``matcher(record, genotype)`` once per case. Inputs are copied so a matcher cannot alter them."""
    results = []
    for case in cases if cases is not None else load_cases():
        got = bool(matcher(copy.deepcopy(case["record"]), copy.deepcopy(case["genotype"])))
        results.append(Result(case["id"], case["expect"] == "match", got, case["why"]))
    return results


def failed_ids(results: list[Result]) -> set[str]:
    return {r.case_id for r in results if not r.passed}


def format_report(results: list[Result]) -> str:
    width = max(len(r.case_id) for r in results)
    lines = []
    for r in results:
        mark = "PASS" if r.passed else "FAIL"
        want = "match" if r.expect_match else "no match"
        got = "match" if r.got_match else "no match"
        detail = "" if r.passed else f"  (expected {want}, got {got})"
        lines.append(f"  [{mark}] {r.case_id.ljust(width)}{detail}".rstrip())
    passed = sum(r.passed for r in results)
    lines.append(f"\n{passed}/{len(results)} cases passed.")
    return "\n".join(lines)
