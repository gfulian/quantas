#!/usr/bin/env python3
"""Generate page-local scientific bibliographies from ``quantas.references``."""

from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))


PAGE_REFERENCE_KEYS: dict[str, tuple[str, ...]] = {
    "theory/ha.rst": (
        "mcquarrie_simon_1997",
    ),
    "theory/qha.rst": (
        "anderson_1995",
        "anderson_masuda_isaak_1995",
        "erba_2014",
        "erba_shahrokhi_moradian_dovesi_2015",
    ),
    "theory/elasticity.rst": (
        "nye_1985",
        "mouhat_coudert_2014",
        "hill_1952",
        "elate_gaillac_pullumbi_coudert_2016",
        "stixrude_lithgow_bertelloni_2005",
    ),
    "theory/seismic.rst": (
        "nye_1985",
        "jaeken_cottenier_2016",
        "seismic_ulian_valdre_2024",
    ),
    "theory/eos.rst": (
        "eosfit7_angel_gonzalez_platas_alvaro_2014",
        "sjeos_alchagirov_perdew_boettger_albers_fiolhais_2001",
        "staroverov_scuseria_tao_perdew_2004",
        "anderson_1995",
        "stixrude_lithgow_bertelloni_2005",
    ),
    "theory/thermoelasticity.rst": (
        "stixrude_lithgow_bertelloni_2005",
        "wallace_1972",
        "barron_klein_1965",
        "destefanis_ravoux_cossard_erba_2019",
        "davies_1974",
        "waters_bielawski_2016",
    ),
    "theory/earth_profiles.rst": (
        "prem_dziewonski_anderson_1981",
        "hasterok_chapman_2011",
        "parsons_sclater_1977",
        "katsura_2022",
        "katsura_software_2022",
    ),
    "introduction/citing_quantas.rst": (
        "quantas_2022",
        "seismic_ulian_valdre_2024",
        "mcquarrie_simon_1997",
        "anderson_1995",
        "anderson_masuda_isaak_1995",
        "erba_2014",
        "erba_shahrokhi_moradian_dovesi_2015",
        "qha_ulian_valdre_2018",
        "eosfit7_angel_gonzalez_platas_alvaro_2014",
        "sjeos_alchagirov_perdew_boettger_albers_fiolhais_2001",
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
        "prem_dziewonski_anderson_1981",
        "hasterok_chapman_2011",
        "parsons_sclater_1977",
        "katsura_2022",
        "katsura_software_2022",
        "destefanis_ravoux_cossard_erba_2019",
        "erba_mahmoud_belmonte_dovesi_2014",
        "stixrude_lithgow_bertelloni_2005",
        "davies_1974",
        "waters_bielawski_2016",
        "wallace_1972",
        "karki_stixrude_clark_warren_ackland_crain_1997",
        "sinogeikin_bass_1999",
        "jiang_speziale_duffy_2006",
        "nye_1985",
        "hill_1952",
        "elate_gaillac_pullumbi_coudert_2016",
        "jaeken_cottenier_2016",
    ),
    "validation/thermoelasticity.rst": (
        "karki_stixrude_clark_warren_ackland_crain_1997",
        "sinogeikin_bass_1999",
        "jiang_speziale_duffy_2006",
    ),
    "developer/interfaces.rst": (
        "erba_mahmoud_belmonte_dovesi_2014",
    ),
}


def fragment_name(page: str) -> str:
    """Return the generated fragment name for one documentation page."""
    path = Path(page)
    if path.parts[0] == "theory":
        return f"{path.stem}.inc"
    return f"{'_'.join(path.with_suffix('').parts)}.inc"


def generate(output_root: Path | None = None) -> tuple[Path, ...]:
    """Generate all page-local bibliography fragments."""
    from quantas.references.render import render_rst_bibliography

    root = output_root or ROOT / "docs" / "source" / "_generated" / "references"
    root.mkdir(parents=True, exist_ok=True)
    generated: list[Path] = []
    for page, keys in PAGE_REFERENCE_KEYS.items():
        path = root / fragment_name(page)
        path.write_text(render_rst_bibliography(keys), encoding="utf-8")
        generated.append(path)
    return tuple(generated)


if __name__ == "__main__":
    for generated_path in generate():
        print(generated_path.relative_to(ROOT))
