# -*- coding: utf-8 -*-

"""Canonical bibliography and citation sets for Quantas."""

from .models import Citation, CitationKind
from .registry import CITATIONS, get_citation
from .render import (
    render_citation,
    render_citation_inline,
    render_citation_list,
    render_citation_notice,
    render_rst_bibliography,
    render_rst_footnote,
)
from .sets import (
    METHOD_CITATION_KEYS,
    MODULE_CITATION_KEYS,
    method_citation_keys,
    module_citation_keys,
)

__all__ = [
    "CITATIONS",
    "METHOD_CITATION_KEYS",
    "MODULE_CITATION_KEYS",
    "Citation",
    "CitationKind",
    "get_citation",
    "method_citation_keys",
    "module_citation_keys",
    "render_citation",
    "render_citation_inline",
    "render_citation_list",
    "render_citation_notice",
    "render_rst_bibliography",
    "render_rst_footnote",
]
