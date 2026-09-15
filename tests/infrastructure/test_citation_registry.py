"""Characterization tests for the canonical Quantas citation registry."""

from __future__ import annotations

from dataclasses import replace

import pytest

from quantas.references import (
    CITATIONS,
    CitationKind,
    get_citation,
    render_citation,
    render_rst_footnote,
    validate_citation_registry,
)
from quantas.references.sets import METHOD_CITATION_KEYS


def test_registry_is_valid_and_complete_for_b12() -> None:
    """The b12 registry passes its structural validation contract."""
    validate_citation_registry()
    assert len(CITATIONS) == 44


def test_author_style_uses_initials_and_unicode() -> None:
    """Canonical metadata uses initials plus scientific-name typography."""
    quantas = get_citation("quantas_2022")
    eosfit = get_citation("eosfit7_angel_gonzalez_platas_alvaro_2014")

    assert quantas.authors == ("G. Ulian", "G. Valdrè")
    assert eosfit.authors == ("R. J. Angel", "J. Gonzalez-Platas", "M. Alvaro")
    assert eosfit.journal == "Zeitschrift für Kristallographie"


def test_plain_text_renderer_is_ascii_but_rst_preserves_typography() -> None:
    """Report text stays portable while scientific RST keeps Unicode metadata."""
    plain = render_citation(get_citation("quantas_2022"))
    rst = render_rst_footnote("eosfit7_angel_gonzalez_platas_alvaro_2014")

    assert plain.isascii()
    assert "G. Valdre" in plain
    assert "Valdrè" not in plain
    assert "Zeitschrift für Kristallographie" in rst
    stixrude = render_citation(get_citation("stixrude_lithgow_bertelloni_2005"))
    assert "minerals-I. Physical properties" in stixrude


def test_chapter_and_report_kinds_are_explicit() -> None:
    """Book chapters and technical reports are not forced into article records."""
    fei = get_citation("fei_1995")
    odrpack = get_citation("boggs_byrd_rogers_schnabel_1992")

    assert fei.kind is CitationKind.CHAPTER
    assert fei.container_title == (
        "Mineral Physics & Crystallography: A Handbook of Physical Constants"
    )
    assert fei.editors == ("T. J. Ahrens",)
    assert fei.doi == "10.1029/RF002p0029"

    assert odrpack.kind is CitationKind.REPORT
    assert odrpack.report_number == "NISTIR 4834"
    assert odrpack.doi == "10.6028/NIST.IR.4834"


def test_b12_inventory_contains_manual_documentation_references() -> None:
    """References queued for the documentation migration already resolve."""
    expected = {
        "orear_1982",
        "boggs_byrd_rogers_schnabel_1992",
        "zwolak_boggs_watson_2007",
        "berman_1988",
        "fei_1995",
        "pawley_redfern_holland_1996",
        "salje_wruck_thomas_1991",
        "holland_powell_2011",
        "helffrich_connolly_2009",
        "kroll_kirfel_heinemann_barbier_2012",
        "erba_mahmoud_belmonte_dovesi_2014",
        "karki_stixrude_clark_warren_ackland_crain_1997",
        "sinogeikin_bass_1999",
        "jiang_speziale_duffy_2006",
    }
    assert expected <= CITATIONS.keys()


def test_new_method_sets_reference_only_canonical_records() -> None:
    """New solver and CRYSTAL method sets are registry-backed."""
    assert METHOD_CITATION_KEYS["effective_variance_weighting"] == ("orear_1982",)
    assert METHOD_CITATION_KEYS["orthogonal_distance_regression"] == (
        "boggs_byrd_rogers_schnabel_1992",
        "zwolak_boggs_watson_2007",
    )
    assert METHOD_CITATION_KEYS["crystal_finite_pressure_elasticity"] == (
        "erba_mahmoud_belmonte_dovesi_2014",
    )


def test_registry_validation_rejects_duplicate_doi() -> None:
    """Two canonical keys cannot silently identify the same publication."""
    original = get_citation("quantas_2022")
    duplicate = replace(original, key="duplicate_2022")

    with pytest.raises(ValueError, match="duplicates DOI"):
        validate_citation_registry(
            {
                original.key: original,
                duplicate.key: duplicate,
            }
        )
