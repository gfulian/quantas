# -*- coding: utf-8 -*-

"""Central CRYSTAL output markers used by Quantas interface parsers."""

from __future__ import annotations

import re


FLOAT = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"
FLOAT_RE = re.compile(FLOAT)
ATOM_COUNT_RE = re.compile(r"ATOMS IN THE UNIT CELL:\s*(?P<count>\d+)")
SPACE_GROUP_RE = re.compile(r"SPACE\s+GROUP\s+N\.\s*:\s*(?P<number>\d+)")

NORMAL_TERMINATION_RE = re.compile(r"^\s*E{10}\s+TERMINATION\b", re.IGNORECASE)

SCF_CYCLE_RE = re.compile(
    rf"^\s*CYC\s+(?P<cycle>\d+)\s+ETOT\(AU\)\s+"
    rf"(?P<energy>{FLOAT})\s+DETOT\s*(?P<delta>{FLOAT})\b",
    re.IGNORECASE,
)
SCF_END_RE = re.compile(
    r"^\s*==\s*SCF\s+ENDED\s*-\s*(?P<reason>.+)$", re.IGNORECASE
)
SCF_END_ENERGY_RE = re.compile(
    rf"\bE\(AU\)\s+(?P<energy>{FLOAT})\s+CYCLES\s+(?P<cycles>\d+)\b",
    re.IGNORECASE,
)

TOTAL_DFT_ENERGY_RE = re.compile(
    rf"TOTAL\s+ENERGY\(DFT\)\(AU\)\(\s*(?P<cycles>\d+)\)\s*"
    rf"(?P<energy>{FLOAT})",
    re.IGNORECASE,
)
TOTAL_DFT_DELTA_RE = re.compile(
    rf"\bDE(?:\s*\(AU\))?\s*=?\s*(?P<delta>{FLOAT})",
    re.IGNORECASE,
)
CORRECTED_TOTAL_ENERGY_RE = re.compile(
    rf"^\s*TOTAL\s+ENERGY\s*\+.*?(?P<energy>{FLOAT})\s*$",
    re.IGNORECASE,
)
CENTRAL_POINT_RE = re.compile(
    rf"^\s*CENTRAL\s+POINT\s+(?P<energy>{FLOAT})\b",
    re.IGNORECASE,
)

OPT_POINT_RE = re.compile(
    r"^\s*(?P<label>[A-Z][A-Z ]*?)\s+OPTIMIZATION\s*-\s*POINT\s+"
    r"(?P<point>\d+)\b",
    re.IGNORECASE,
)
OPT_END_RE = re.compile(
    rf"^\s*\*\s*OPT\s+END\s*-\s*(?P<reason>[^*]+?)\s*\*"
    rf".*?E\(AU\):\s*(?P<energy>{FLOAT}).*?POINTS\s+(?P<points>\d+)\s*\*",
    re.IGNORECASE,
)
MAX_GRADIENT_RE = re.compile(
    rf"^\s*MAX\s+GRADIENT\s+(?P<value>{FLOAT})\b", re.IGNORECASE
)
RMS_GRADIENT_RE = re.compile(
    rf"^\s*RMS\s+GRADIENT\s+(?P<value>{FLOAT})\b", re.IGNORECASE
)
MAX_DISPLACEMENT_RE = re.compile(
    rf"^\s*MAX\s+DISPLAC\.\s+(?P<value>{FLOAT})\b", re.IGNORECASE
)
RMS_DISPLACEMENT_RE = re.compile(
    rf"^\s*RMS\s+DISPLAC\.\s+(?P<value>{FLOAT})\b", re.IGNORECASE
)
