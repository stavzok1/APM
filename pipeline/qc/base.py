from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any, Dict, Optional


class Level(str, Enum):
    PASS = "pass"
    WARN = "warn"
    FAIL = "fail"


@dataclass(frozen=True)
class Finding:
    check: str
    level: Level
    message: str
    context: Dict[str, Any]
    hint: Optional[str] = None

