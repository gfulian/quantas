# -*- coding: utf-8 -*-

"""Deterministic renderers for canonical Quantas citations."""

from __future__ import annotations

from collections.abc import Iterable

from .models import Citation, CitationKind
from .registry import get_citation


def render_citation(citation: Citation) -> str:
    """Render one citation as deterministic plain text."""
    authors = ", ".join(citation.authors)
    source_parts: list[str] = []
    if citation.journal:
        source_parts.append(citation.journal)
    if citation.volume:
        source_parts.append(citation.volume)
    if citation.pages:
        source_parts.append(citation.pages)
    if citation.publisher:
        source_parts.append(citation.publisher)
    source = ", ".join(source_parts)
    lines = [f"{authors} ({citation.year}).", f"{citation.title}."]
    if source:
        lines.append(source + ".")
    if citation.doi:
        lines.append(f"https://doi.org/{citation.doi}")
    elif citation.url:
        lines.append(citation.url)
    return "\n".join(lines)


def render_citation_inline(citation: Citation | str) -> str:
    """Render one canonical citation on a single deterministic line.

    Parameters
    ----------
    citation : Citation or str
        Citation record or registered citation key.

    Returns
    -------
    str
        Single-line bibliographic representation.
    """
    record = get_citation(citation) if isinstance(citation, str) else citation
    return " ".join(render_citation(record).splitlines())


def render_citation_list(keys: Iterable[str]) -> str:
    """Render an ordered list of registered citations."""
    records = [render_citation(get_citation(key)) for key in keys]
    return "\n\n".join(records)


def _render_rst_authors(authors: tuple[str, ...]) -> str:
    """Render an ordered author list for an RST bibliography entry."""
    if not authors:
        return ""
    if len(authors) == 1:
        return authors[0]
    if len(authors) == 2:
        return f"{authors[0]} and {authors[1]}"
    return ", ".join(authors[:-1]) + f", and {authors[-1]}"


def render_rst_footnote(citation: Citation | str) -> str:
    """Render one canonical citation as a labelled auto-numbered RST footnote.

    Parameters
    ----------
    citation : Citation or str
        Citation record or registered citation key.

    Returns
    -------
    str
        A reStructuredText footnote definition.  The stable citation key is
        used as the label, while Sphinx displays a page-local number.
    """
    record = get_citation(citation) if isinstance(citation, str) else citation
    authors = _render_rst_authors(record.authors)
    if record.kind is CitationKind.BOOK:
        body = f"{authors}. *{record.title}*."
        if record.publisher:
            body += f" {record.publisher}"
        body += f" ({record.year})."
    elif record.kind is CitationKind.SOFTWARE:
        body = f"{authors}. *{record.title}*."
        if record.publisher:
            body += f" {record.publisher}"
        body += f" ({record.year})."
    else:
        body = f'{authors}. “{record.title}”.'
        source: list[str] = []
        if record.journal:
            source.append(f"*{record.journal}*")
        if record.volume:
            source.append(f"**{record.volume}**")
        if record.pages:
            source.append(record.pages)
        if source:
            body += " " + ", ".join(source)
        body += f" ({record.year})."
    if record.doi:
        body += (
            f" `DOI: {record.doi} "
            f"<https://doi.org/{record.doi}>`_."
        )
    elif record.url:
        body += f" `External link <{record.url}>`_."
    return f".. [#{record.key}] {body}"


def render_rst_bibliography(
    keys: Iterable[str],
    *,
    heading: str = "Bibliographic references",
) -> str:
    """Render a page-local numbered RST bibliography from canonical keys.

    Parameters
    ----------
    keys : iterable of str
        Canonical citation keys in first-appearance order.
    heading : str, optional
        Section heading placed above the footnote definitions.

    Returns
    -------
    str
        Complete reStructuredText bibliography section.
    """
    records = [render_rst_footnote(key) for key in keys]
    underline = "-" * len(heading)
    return f"{heading}\n{underline}\n\n" + "\n\n".join(records) + "\n"


def render_citation_notice(keys: Iterable[str]) -> str:
    """Render the standard Quantas citation notice and reference list."""
    separator = "_" * 80
    body = render_citation_list(keys)
    return (
        f"{separator}\n"
        "The methods used in this calculation are described by the following "
        "references.\n"
        "Please cite Quantas and the relevant scientific methods when publishing "
        "results\n"
        "derived from this calculation.\n\n"
        "Quantas is academic, open-source software. Scientific recognition within "
        "the\n"
        "community is essential to its continued development. Your support is "
        "greatly\n"
        "appreciated. Thank you!\n\n"
        "References\n"
        "----------\n"
        f"{body}\n"
        f"{separator}"
    )


__all__ = [
    "render_citation",
    "render_citation_inline",
    "render_citation_list",
    "render_citation_notice",
    "render_rst_bibliography",
    "render_rst_footnote",
]
