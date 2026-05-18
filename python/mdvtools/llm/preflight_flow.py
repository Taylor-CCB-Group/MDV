from __future__ import annotations

from collections.abc import Callable
from typing import Any

from .code_preflight import format_preflight_issues, validate_generated_code_preflight


def preflight_with_single_retry(
    *,
    initial_code: str,
    regenerate_once: Callable[[str], str],
    log: Callable[[str], Any],
    datasource_fields: dict[str, set[str]] | None = None,
    allowed_wrapper_subgroup_keys: set[str] | None = None,
) -> tuple[str, dict[str, object]]:
    """
    Validate generated code before execution and allow one guided regeneration.
    """
    first = validate_generated_code_preflight(
        initial_code,
        datasource_fields=datasource_fields,
        allowed_wrapper_subgroup_keys=allowed_wrapper_subgroup_keys,
    )
    metadata: dict[str, object] = {
        "preflight_ok": first.ok,
        "preflight_retried": False,
        "preflight_errors": [i.message for i in first.issues],
    }
    if first.ok:
        return initial_code, metadata

    issue_text = format_preflight_issues(first.issues)
    log(f"Preflight failed before execution:\n{issue_text}")
    retry_code = regenerate_once(issue_text)
    second = validate_generated_code_preflight(
        retry_code,
        datasource_fields=datasource_fields,
        allowed_wrapper_subgroup_keys=allowed_wrapper_subgroup_keys,
    )
    metadata["preflight_retried"] = True
    metadata["preflight_ok"] = second.ok
    metadata["preflight_errors"] = [i.message for i in second.issues]
    if not second.ok:
        raise Exception(
            "Code preflight failed after one retry:\n"
            f"{format_preflight_issues(second.issues)}"
        )
    return retry_code, metadata

