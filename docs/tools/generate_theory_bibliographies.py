#!/usr/bin/env python3
"""Generate page-local scientific bibliographies from ``quantas.references``."""

from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))


THEORY_REFERENCE_KEYS: dict[str, tuple[str, ...]] = {
    "ha": (
        "mcquarrie_simon_1997",
    ),
    "qha": (
        "anderson_1995",
        "anderson_masuda_isaak_1995",
        "erba_2014",
        "erba_shahrokhi_moradian_dovesi_2015",
    ),
    "elasticity": (
        "nye_1985",
        "mouhat_coudert_2014",
        "hill_1952",
        "elate_gaillac_pullumbi_coudert_2016",
        "stixrude_lithgow_bertelloni_2005",
    ),
    "seismic": (
        "nye_1985",
        "jaeken_cottenier_2016",
        "seismic_ulian_valdre_2024",
    ),
    "eos": (
        "eosfit7_angel_gonzalez_platas_alvaro_2014",
        "anderson_1995",
        "stixrude_lithgow_bertelloni_2005",
    ),
    "thermoelasticity": (
        "stixrude_lithgow_bertelloni_2005",
        "wallace_1972",
        "barron_klein_1965",
        "destefanis_ravoux_cossard_erba_2019",
        "davies_1974",
        "waters_bielawski_2016",
    ),
    "earth_profiles": (
        "prem_dziewonski_anderson_1981",
        "hasterok_chapman_2011",
        "parsons_sclater_1977",
        "katsura_2022",
        "katsura_software_2022",
    ),
}


def generate(output_root: Path | None = None) -> tuple[Path, ...]:
    """Generate all scientific-background bibliography fragments."""
    from quantas.references.render import render_rst_bibliography

    root = output_root or ROOT / "docs" / "source" / "_generated" / "references"
    root.mkdir(parents=True, exist_ok=True)
    generated: list[Path] = []
    for page, keys in THEORY_REFERENCE_KEYS.items():
        path = root / f"{page}.inc"
        path.write_text(render_rst_bibliography(keys), encoding="utf-8")
        generated.append(path)
    return tuple(generated)


if __name__ == "__main__":
    for generated_path in generate():
        print(generated_path.relative_to(ROOT))
