"""Declarative workflows (Task 16.1.2)."""

from __future__ import annotations

import json
from typing import Iterable

from ._xtalforge import workflow_run as _run
from ._xtalforge import workflow_schema as _schema
from ._xtalforge import workflow_validate as _validate


class Workflow:
    """A versioned xtalforge workflow document."""

    def __init__(self, text: str) -> None:
        self.text = text

    @staticmethod
    def schema() -> str:
        """Returns the published JSON Schema."""
        return _schema()

    def validate(self) -> dict[str, object]:
        """Parses and validates the workflow."""
        return json.loads(_validate(self.text))

    def run(
        self,
        inputs: Iterable[str],
        *,
        jobs: int = 1,
        checkpoint: str | None = None,
    ) -> dict[str, object]:
        """Runs the workflow over files and returns the JSON report."""
        return json.loads(_run(self.text, list(inputs), jobs, checkpoint))


__all__ = ["Workflow"]
