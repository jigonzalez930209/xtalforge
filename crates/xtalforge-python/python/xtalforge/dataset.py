"""JSON Lines dataset view (Task 16.3.1)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable, Iterator

Record = dict[str, object]


class Dataset:
    """A JSON Lines dataset of structures and properties."""

    def __init__(self, records: Iterable[Record]) -> None:
        self.records = [dict(record) for record in records]

    @classmethod
    def from_jsonl(cls, text: str) -> "Dataset":
        records = [json.loads(line) for line in text.splitlines() if line.strip()]
        return cls(records)

    @classmethod
    def load(cls, path: str | Path) -> "Dataset":
        return cls.from_jsonl(Path(path).read_text(encoding="utf-8"))

    def to_jsonl(self) -> str:
        return "".join(json.dumps(record, sort_keys=True) + "\n" for record in self.records)

    def save(self, path: str | Path) -> None:
        Path(path).write_text(self.to_jsonl(), encoding="utf-8")

    def __len__(self) -> int:
        return len(self.records)

    def __iter__(self) -> Iterator[Record]:
        return iter(self.records)

    def to_pandas(self) -> object:
        """Returns a pandas DataFrame; pandas is optional."""
        import pandas

        return pandas.DataFrame(self.records)

    def to_polars(self) -> object:
        """Returns a polars DataFrame; polars is optional."""
        import polars

        return polars.DataFrame(self.records)


__all__ = ["Dataset", "Record"]
