"""DEPRECATED shim — pathway interpretation moved into evidence.py.

Kept only so old imports keep working. `rm src/pathways.py` locally to remove it
(the sandbox cannot delete files). New code should import from evidence.
"""
from __future__ import annotations

from .evidence import DEFAULT_LIBRARIES, enrich  # noqa: F401
