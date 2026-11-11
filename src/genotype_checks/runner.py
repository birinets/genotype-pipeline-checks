# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Ruslan Plakhuta
"""Load the cases and run a matcher over them."""

from __future__ import annotations

import copy
import json
from dataclasses import dataclass
from importlib import resources
from typing import Any, Callable, Union

from .matchers import MATCH, OUTCOMES

Matcher = Callable[[dict, dict], Any]
Answer = Union[bool, str]


class UnknownOutcome(ValueError):
    """A matcher answered with a name that is not one of the outcomes."""


def load_cases() -> list[dict[str, Any]]:
    """Return the bundled cases."""
    path = resources.files("genotype_checks").joinpath("cases.json")
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)["cases"]


@dataclass(frozen=True)
class Result:
    """One case's result. ``got`` is True/False for a yes/no matcher, or an outcome name."""

    case_id: str
    expect: str
    got: Answer
    why: str

    @property
    def expect_match(self) -> bool:
        return self.expect == MATCH

    @property
    def passed(self) -> bool:
        if isinstance(self.got, bool):
            return self.got == self.expect_match
        return self.got == self.expect


def _as_answer(value: Any) -> Answer:
    """Read a matcher's answer: an outcome name, or anything else as yes/no.

    A tuple, such as an ``(outcome, reason)`` pair, is read by its first element. A name that is not
    one of ``OUTCOMES`` raises ``UnknownOutcome`` rather than failing every case quietly.
    """
    if isinstance(value, tuple) and value:
        value = value[0]
    if isinstance(value, str):
        if value not in OUTCOMES:
            raise UnknownOutcome(f"unknown outcome {value!r}; expected one of {', '.join(OUTCOMES)}")
        return value
    return bool(value)


def run_cases(matcher: Matcher, cases: list[dict[str, Any]] | None = None) -> list[Result]:
    """Call ``matcher(record, genotype)`` once per case. Inputs are copied so a matcher cannot alter them.

    A yes/no answer passes when it is yes exactly for the cases whose outcome is ``match``. An
    outcome name passes only when it is the case's outcome. Any other name raises ``UnknownOutcome``.
    """
    results = []
    for case in cases if cases is not None else load_cases():
        got = _as_answer(matcher(copy.deepcopy(case["record"]), copy.deepcopy(case["genotype"])))
        results.append(Result(case["id"], case["expect"], got, case["why"]))
    return results


def failed_ids(results: list[Result]) -> set[str]:
    return {r.case_id for r in results if not r.passed}


def _describe(got: Answer) -> str:
    if isinstance(got, bool):
        return "match" if got else "no match"
    return got


def format_report(results: list[Result]) -> str:
    width = max(len(r.case_id) for r in results)
    lines = []
    for r in results:
        mark = "PASS" if r.passed else "FAIL"
        detail = "" if r.passed else f"  (expected {r.expect}, got {_describe(r.got)})"
        lines.append(f"  [{mark}] {r.case_id.ljust(width)}{detail}".rstrip())
    if all(isinstance(r.got, bool) for r in results):
        lines.append("  Yes/no answers: a case expects yes only when its outcome is match.")
    passed = sum(r.passed for r in results)
    lines.append(f"\n{passed}/{len(results)} cases passed.")
    return "\n".join(lines)
