from typing import Any

class StructureWidget:
    scene: dict[str, object]
    frame: int
    frames: int

    def __init__(self, structure: str, trajectory: list[str] | None = None) -> None: ...
    def show_frame(self, index: int) -> None: ...
