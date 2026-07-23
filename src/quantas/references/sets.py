# -*- coding: utf-8 -*-

"""Named citation sets for Quantas modules and scientific methods."""

from __future__ import annotations

MODULE_CITATION_KEYS: dict[str, tuple[str, ...]] = {
    "ha": (
        "quantas_2022",
        "mcquarrie_simon_1997",
    ),
    "qha": (
        "quantas_2022",
        "qha_ulian_valdre_2018",
        "anderson_1995",
        "anderson_masuda_isaak_1995",
        "erba_2014",
        "erba_shahrokhi_moradian_dovesi_2015",
    ),
    "eos": (
        "quantas_2022",
        "eos_ulian_tosoni_valdre_2014",
        "eosfit7_angel_gonzalez_platas_alvaro_2014",
    ),
    "elasticity": (
        "quantas_2022",
        "elasticity_ulian_moro_valdre_2018",
        "nye_1985",
        "hill_1952",
        "elate_gaillac_pullumbi_coudert_2016",
        "mouhat_coudert_2014",
    ),
    "seismic": (
        "quantas_2022",
        "seismic_ulian_valdre_2024",
        "jaeken_cottenier_2016",
        "nye_1985",
    ),
    "thermoelasticity": (
        "quantas_2022",
        "stixrude_lithgow_bertelloni_2005",
        "destefanis_ravoux_cossard_erba_2019",
        "barron_klein_1965",
        "wallace_1972",
        "davies_1974",
        "waters_bielawski_2016",
        "mouhat_coudert_2014",
    ),
    "earth_profiles": (
        "prem_dziewonski_anderson_1981",
        "hasterok_chapman_2011",
        "parsons_sclater_1977",
        "katsura_2022",
        "katsura_software_2022",
    ),
}

METHOD_CITATION_KEYS: dict[str, tuple[str, ...]] = {
    "linear_elasticity": ("nye_1985",),
    "voigt_reuss_hill": ("hill_1952",),
    "directional_elasticity": (
        "nye_1985",
        "elate_gaillac_pullumbi_coudert_2016",
    ),
    "christoffel_acoustics": (
        "jaeken_cottenier_2016",
        "seismic_ulian_valdre_2024",
    ),
    "harmonic_statistical_thermodynamics": ("mcquarrie_simon_1997",),
    "quasi_harmonic_approximation": (
        "anderson_1995",
        "anderson_masuda_isaak_1995",
        "erba_2014",
        "erba_shahrokhi_moradian_dovesi_2015",
    ),
    "isothermal_equations_of_state": (
        "eosfit7_angel_gonzalez_platas_alvaro_2014",
        "anderson_1995",
    ),
    "thermal_expansion_equations": (
        "eosfit7_angel_gonzalez_platas_alvaro_2014",
    ),
    "pvt_equations_of_state": (
        "eosfit7_angel_gonzalez_platas_alvaro_2014",
        "anderson_1995",
    ),
    "mie_gruneisen_debye": (
        "anderson_1995",
        "stixrude_lithgow_bertelloni_2005",
    ),
    "wallace_stress_strain": ("barron_klein_1965", "wallace_1972"),
    "cold_finite_strain": ("stixrude_lithgow_bertelloni_2005",),
    "quasi_static_thermoelasticity": (
        "destefanis_ravoux_cossard_erba_2019",
        "stixrude_lithgow_bertelloni_2005",
    ),
    "adiabatic_elasticity": (
        "wallace_1972",
        "davies_1974",
        "waters_bielawski_2016",
    ),
    "mechanical_stability": (
        "mouhat_coudert_2014",
        "barron_klein_1965",
        "wallace_1972",
    ),
}


def module_citation_keys(module: str) -> tuple[str, ...]:
    """Return the canonical citation keys for one Quantas module."""
    return MODULE_CITATION_KEYS.get(module, ("quantas_2022",))


def method_citation_keys(method: str) -> tuple[str, ...]:
    """Return canonical citation keys for one scientific method."""
    return METHOD_CITATION_KEYS.get(method, ())


__all__ = [
    "METHOD_CITATION_KEYS",
    "MODULE_CITATION_KEYS",
    "method_citation_keys",
    "module_citation_keys",
]
