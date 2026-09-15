# -*- coding: utf-8 -*-

"""Deterministic renderers for canonical Quantas citations."""

from __future__ import annotations

from collections.abc import Iterable
import unicodedata

from .models import Citation, CitationKind
from .registry import get_citation


def _ascii_text(value: str) -> str:
    """Return a portable ASCII representation for plain-text surfaces."""
    punctuation = str.maketrans(
        {
            "–": "-",
            "—": "-",
            "−": "-",
            "‘": "'",
            "’": "'",
            "“": '"',
            "”": '"',
        }
    )
    normalized = unicodedata.normalize("NFKD", value.translate(punctuation))
    return normalized.encode("ascii", errors="ignore").decode("ascii")


def _plain_authors(authors: tuple[str, ...]) -> str:
    """Render canonical author strings for ASCII report output."""
    return ", ".join(_ascii_text(author) for author in authors)


def _plain_source(citation: Citation) -> str:
    """Render source metadata for deterministic plain-text output."""
    parts: list[str] = []
    if citation.kind is CitationKind.CHAPTER:
        if citation.container_title:
            parts.append(f"In {_ascii_text(citation.container_title)}")
        if citation.editors:
            parts.append("edited by " + _plain_authors(citation.editors))
    elif citation.journal:
        parts.append(_ascii_text(citation.journal))
    if citation.volume:
        parts.append(citation.volume)
    if citation.pages:
        parts.append(citation.pages)
    if citation.report_number:
        parts.append(citation.report_number)
    if citation.publisher:
        parts.append(_ascii_text(citation.publisher))
    return ", ".join(parts)


def render_citation(citation: Citation) -> str:
    """Render one citation as deterministic ASCII plain text.

    Parameters
    ----------
    citation : Citation
        Canonical bibliographic record.

    Returns
    -------
    str
        Portable plain-text representation suitable for terminal reports and
        HDF5-embedded report text.
    """
    authors = _plain_authors(citation.authors)
    source = _plain_source(citation)
    lines = [f"{authors} ({citation.year}).", f"{_ascii_text(citation.title)}."]
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


def _render_rst_source(record: Citation) -> str:
    """Render bibliographic source metadata for an RST footnote."""
    if record.kind is CitationKind.CHAPTER:
        parts: list[str] = []
        if record.container_title:
            parts.append(f"In *{record.container_title}*")
        if record.editors:
            parts.append(f"edited by {_render_rst_authors(record.editors)}")
        if record.volume:
            parts.append(f"**{record.volume}**")
        if record.pages:
            parts.append(record.pages)
        if record.publisher:
            parts.append(record.publisher)
        return ", ".join(parts)

    if record.kind is CitationKind.REPORT:
        parts = []
        if record.report_number:
            parts.append(record.report_number)
        if record.publisher:
            parts.append(record.publisher)
        return ", ".join(parts)

    parts = []
    if record.journal:
        parts.append(f"*{record.journal}*")
    if record.volume:
        parts.append(f"**{record.volume}**")
    if record.pages:
        parts.append(record.pages)
    if record.publisher:
        parts.append(record.publisher)
    return ", ".join(parts)


def render_rst_footnote(citation: Citation | str) -> str:
    """Render one canonical citation as a labelled auto-numbered RST footnote.

    Parameters
    ----------
    citation : Citation or str
        Citation record or registered citation key.

    Returns
    -------
    str
        A reStructuredText footnote definition. The stable citation key is
        used as the label, while Sphinx displays a page-local number.
    """
    record = get_citation(citation) if isinstance(citation, str) else citation
    authors = _render_rst_authors(record.authors)
    source = _render_rst_source(record)

    if record.kind in {CitationKind.BOOK, CitationKind.REPORT, CitationKind.SOFTWARE}:
        body = f"{authors}. *{record.title}*."
    else:
        body = f'{authors}. “{record.title}”.'

    if source:
        body += f" {source}"
    body += f" ({record.year})."

    if record.doi:
        body += f" `DOI: {record.doi} <https://doi.org/{record.doi}>`_."
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
