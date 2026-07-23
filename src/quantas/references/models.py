# -*- coding: utf-8 -*-

"""Passive bibliographic data contracts used by Quantas."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class CitationKind(str, Enum):
    """Supported bibliographic record categories."""

    ARTICLE = "article"
    BOOK = "book"
    PREPRINT = "preprint"
    SOFTWARE = "software"


@dataclass(frozen=True, slots=True)
class Citation:
    """Canonical bibliographic record.

    Parameters
    ----------
    key : str
        Stable machine-readable identifier.
    authors : tuple of str
        Ordered author names.
    title : str
        Publication title.
    year : int
        Publication year.
    kind : CitationKind, optional
        Bibliographic record category.
    journal : str or None, optional
        Journal or series title.
    volume : str or None, optional
        Journal volume.
    pages : str or None, optional
        Page range or article number.
    publisher : str or None, optional
        Publisher for books or software.
    doi : str or None, optional
        Digital object identifier without a URL prefix.
    url : str or None, optional
        Stable external URL when no DOI is available.
    """

    key: str
    authors: tuple[str, ...]
    title: str
    year: int
    kind: CitationKind = CitationKind.ARTICLE
    journal: str | None = None
    volume: str | None = None
    pages: str | None = None
    publisher: str | None = None
    doi: str | None = None
    url: str | None = None


__all__ = ["Citation", "CitationKind"]
