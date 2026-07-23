# -*- coding: utf-8 -*-

"""Canonical Quantas bibliography registry."""

from __future__ import annotations

from .models import Citation, CitationKind

QUANTAS_2022 = Citation(
    key="quantas_2022",
    authors=("Gianfranco Ulian", "Giovanni Valdre'"),
    title=(
        "QUANTAS, a Python software for the analysis of solids from ab initio "
        "quantum mechanical simulations and experimental data"
    ),
    year=2022,
    journal="Journal of Applied Crystallography",
    volume="55",
    pages="386-396",
    doi="10.1107/S1600576722000085",
)

QHA_ULIAN_VALDRE_2018 = Citation(
    key="qha_ulian_valdre_2018",
    authors=("Gianfranco Ulian", "Giovanni Valdre'"),
    title=(
        "Equation of state of hexagonal hydroxylapatite (P63) as obtained from "
        "density functional theory simulations"
    ),
    year=2018,
    journal="International Journal of Quantum Chemistry",
    volume="118",
    pages="e25553",
    doi="10.1002/qua.25553",
)

MCQUARRIE_SIMON_1997 = Citation(
    key="mcquarrie_simon_1997",
    authors=("Donald A. McQuarrie", "John D. Simon"),
    title="Physical Chemistry: A Molecular Approach",
    year=1997,
    kind=CitationKind.BOOK,
    publisher="University Science Books, Sausalito, California",
)

ANDERSON_1995 = Citation(
    key="anderson_1995",
    authors=("Orson L. Anderson",),
    title="Equations of State of Solids for Geophysics and Ceramic Science",
    year=1995,
    kind=CitationKind.BOOK,
    volume="31",
    journal="Oxford Monographs on Geology and Geophysics",
    publisher="Oxford University Press, New York",
)

ANDERSON_MASUDA_ISAAK_1995 = Citation(
    key="anderson_masuda_isaak_1995",
    authors=("Orson L. Anderson", "Koji Masuda", "Donald G. Isaak"),
    title="A new thermodynamic approach for high-pressure physics",
    year=1995,
    journal="Physics of the Earth and Planetary Interiors",
    volume="91",
    pages="3-16",
    doi="10.1016/0031-9201(95)03044-W",
)

ERBA_2014 = Citation(
    key="erba_2014",
    authors=("Alessandro Erba",),
    title=(
        "On combining temperature and pressure effects on structural properties "
        "of crystals with standard ab initio techniques"
    ),
    year=2014,
    journal="Journal of Chemical Physics",
    volume="141",
    pages="124115",
    doi="10.1063/1.4896228",
)

ERBA_SHAHROKHI_MORADIAN_DOVESI_2015 = Citation(
    key="erba_shahrokhi_moradian_dovesi_2015",
    authors=(
        "Alessandro Erba",
        "M. Shahrokhi",
        "R. Moradian",
        "Roberto Dovesi",
    ),
    title=(
        "On how differently the quasi-harmonic approximation works for two "
        "isostructural crystals: Thermal properties of periclase and lime"
    ),
    year=2015,
    journal="Journal of Chemical Physics",
    volume="142",
    pages="044114",
    doi="10.1063/1.4906422",
)

EOS_ULIAN_ET_AL_2014 = Citation(
    key="eos_ulian_tosoni_valdre_2014",
    authors=("Gianfranco Ulian", "Sergio Tosoni", "Giovanni Valdre'"),
    title=(
        "The compressional behaviour and the mechanical properties of talc "
        "[Mg3Si4O10(OH)2]: a density functional theory investigation"
    ),
    year=2014,
    journal="Physics and Chemistry of Minerals",
    volume="41",
    pages="639-650",
    doi="10.1007/s00269-014-0677-x",
)

EOSFIT7_ANGEL_ET_AL_2014 = Citation(
    key="eosfit7_angel_gonzalez_platas_alvaro_2014",
    authors=("Ross J. Angel", "Javier Gonzalez-Platas", "Matteo Alvaro"),
    title="EosFit7c and a Fortran module (library) for equation of state calculations",
    year=2014,
    journal="Zeitschrift fuer Kristallographie",
    volume="229",
    pages="405-419",
    doi="10.1515/zkri-2013-1711",
)


NYE_1985 = Citation(
    key="nye_1985",
    authors=("John F. Nye",),
    title="Physical Properties of Crystals: Their Representation by Tensors and Matrices",
    year=1985,
    kind=CitationKind.BOOK,
    publisher="Oxford University Press, Oxford, 2nd edition",
)

HILL_1952 = Citation(
    key="hill_1952",
    authors=("Rodney Hill",),
    title="The elastic behaviour of a crystalline aggregate",
    year=1952,
    journal="Proceedings of the Physical Society. Section A",
    volume="65",
    pages="349-354",
    doi="10.1088/0370-1298/65/5/307",
)

ELASTICITY_ULIAN_ET_AL_2018 = Citation(
    key="elasticity_ulian_moro_valdre_2018",
    authors=("Gianfranco Ulian", "Daniele Moro", "Giovanni Valdre'"),
    title=(
        "First principle investigation of the mechanical properties of natural "
        "mineral layered nanocomposite: clinochlore as a model system"
    ),
    year=2018,
    journal="Composite Structures",
    volume="202",
    pages="551-558",
    doi="10.1016/j.compstruct.2018.02.089",
)

ELATE_GAILLAC_ET_AL_2016 = Citation(
    key="elate_gaillac_pullumbi_coudert_2016",
    authors=("Romain Gaillac", "Pluton Pullumbi", "Francois-Xavier Coudert"),
    title=(
        "ELATE: an open-source online application for analysis and visualization "
        "of elastic tensors"
    ),
    year=2016,
    journal="Journal of Physics: Condensed Matter",
    volume="28",
    pages="275201",
    doi="10.1088/0953-8984/28/27/275201",
)

SEISMIC_ULIAN_VALDRE_2024 = Citation(
    key="seismic_ulian_valdre_2024",
    authors=("Gianfranco Ulian", "Giovanni Valdre'"),
    title=(
        "SEISMIC, a Python-based code of the Quantas package to calculate the "
        "phase and group acoustic velocities in crystals"
    ),
    year=2024,
    journal="Computers & Geosciences",
    volume="188",
    pages="105615",
    doi="10.1016/j.cageo.2024.105615",
)

JAEKEN_COTTENIER_2016 = Citation(
    key="jaeken_cottenier_2016",
    authors=("Jan W. Jaeken", "Stefaan Cottenier"),
    title="Solving the Christoffel equation: Phase and group velocities",
    year=2016,
    journal="Computer Physics Communications",
    volume="207",
    pages="445-451",
    doi="10.1016/j.cpc.2016.06.014",
)

BARRON_KLEIN_1965 = Citation(
    key="barron_klein_1965",
    authors=("T. H. K. Barron", "M. L. Klein"),
    title=("Second-order elastic constants of a solid under stress"),
    year=1965,
    journal="Proceedings of the Physical Society",
    volume="85",
    pages="523-532",
    doi="10.1088/0370-1328/85/3/313",
)

WALLACE_1972 = Citation(
    key="wallace_1972",
    authors=("D. C. Wallace",),
    title="Thermodynamics of Crystals",
    year=1972,
    kind=CitationKind.BOOK,
    publisher="John Wiley & Sons, New York",
)

STIXRUDE_LITHGOW_BERTELLONI_2005 = Citation(
    key="stixrude_lithgow_bertelloni_2005",
    authors=("Lars Stixrude", "Carolina Lithgow-Bertelloni"),
    title="Thermodynamics of mantle minerals - I. Physical properties",
    year=2005,
    journal="Geophysical Journal International",
    volume="162",
    pages="610-632",
    doi="10.1111/j.1365-246X.2005.02642.x",
)

WATERS_BIELAWSKI_2016 = Citation(
    key="waters_bielawski_2016",
    authors=("T. J. Waters", "M. Bielawski"),
    title="Isothermal and adiabatic elastic tensors",
    year=2016,
    kind=CitationKind.PREPRINT,
    journal="arXiv",
    pages="1605.06548",
    doi="10.48550/arXiv.1605.06548",
)


MOUHAT_COUDERT_2014 = Citation(
    key="mouhat_coudert_2014",
    authors=("Felix Mouhat", "Francois-Xavier Coudert"),
    title=(
        "Necessary and sufficient elastic stability conditions in various "
        "crystal systems"
    ),
    year=2014,
    journal="Physical Review B",
    volume="90",
    pages="224104",
    doi="10.1103/PhysRevB.90.224104",
)

DESTEFANIS_RAVOUX_COSSARD_ERBA_2019 = Citation(
    key="destefanis_ravoux_cossard_erba_2019",
    authors=(
        "Maurizio Destefanis",
        "Corentin Ravoux",
        "Alessandro Cossard",
        "Alessandro Erba",
    ),
    title="Thermo-Elasticity of Materials from Quasi-Harmonic Calculations",
    year=2019,
    journal="Minerals",
    volume="9",
    pages="16",
    doi="10.3390/min9010016",
)


DAVIES_1974 = Citation(
    key="davies_1974",
    authors=("G. F. Davies",),
    title="Effective elastic moduli under hydrostatic stress-I. Quasi-harmonic theory",
    year=1974,
    journal="Journal of Physics and Chemistry of Solids",
    volume="35",
    pages="1513-1520",
    doi="10.1016/S0022-3697(74)80279-9",
)

PREM_DZIEWONSKI_ANDERSON_1981 = Citation(
    key="prem_dziewonski_anderson_1981",
    authors=("Adam M. Dziewonski", "Don L. Anderson"),
    title="Preliminary reference Earth model",
    year=1981,
    journal="Physics of the Earth and Planetary Interiors",
    volume="25",
    pages="297-356",
    doi="10.1016/0031-9201(81)90046-7",
)

HASTEROK_CHAPMAN_2011 = Citation(
    key="hasterok_chapman_2011",
    authors=("Derrick Hasterok", "David S. Chapman"),
    title="Heat production and geotherms for the continental lithosphere",
    year=2011,
    journal="Earth and Planetary Science Letters",
    volume="307",
    pages="59-70",
    doi="10.1016/j.epsl.2011.04.034",
)

PARSONS_SCLATER_1977 = Citation(
    key="parsons_sclater_1977",
    authors=("Barry Parsons", "John G. Sclater"),
    title=(
        "An analysis of the variation of ocean floor bathymetry and heat flow with age"
    ),
    year=1977,
    journal="Journal of Geophysical Research",
    volume="82",
    pages="803-827",
    doi="10.1029/JB082i005p00803",
)

KATSURA_2022 = Citation(
    key="katsura_2022",
    authors=("Tomoo Katsura",),
    title="A revised adiabatic temperature profile for the mantle",
    year=2022,
    journal="Journal of Geophysical Research: Solid Earth",
    volume="127",
    pages="e2021JB023562",
    doi="10.1029/2021JB023562",
)

KATSURA_SOFTWARE_2022 = Citation(
    key="katsura_software_2022",
    authors=("Tomoo Katsura",),
    title=(
        "Matlab scripts of 'A revised adiabatic temperature profile for the "
        "mantle' (Version 1.1.0)"
    ),
    year=2022,
    kind=CitationKind.SOFTWARE,
    publisher="Zenodo",
    doi="10.5281/zenodo.5903286",
)

CITATIONS = {
    citation.key: citation
    for citation in (
        QUANTAS_2022,
        QHA_ULIAN_VALDRE_2018,
        MCQUARRIE_SIMON_1997,
        ANDERSON_1995,
        ANDERSON_MASUDA_ISAAK_1995,
        ERBA_2014,
        ERBA_SHAHROKHI_MORADIAN_DOVESI_2015,
        EOS_ULIAN_ET_AL_2014,
        EOSFIT7_ANGEL_ET_AL_2014,
        NYE_1985,
        HILL_1952,
        ELASTICITY_ULIAN_ET_AL_2018,
        ELATE_GAILLAC_ET_AL_2016,
        SEISMIC_ULIAN_VALDRE_2024,
        JAEKEN_COTTENIER_2016,
        BARRON_KLEIN_1965,
        WALLACE_1972,
        STIXRUDE_LITHGOW_BERTELLONI_2005,
        WATERS_BIELAWSKI_2016,
        MOUHAT_COUDERT_2014,
        DESTEFANIS_RAVOUX_COSSARD_ERBA_2019,
        DAVIES_1974,
        PREM_DZIEWONSKI_ANDERSON_1981,
        HASTEROK_CHAPMAN_2011,
        PARSONS_SCLATER_1977,
        KATSURA_2022,
        KATSURA_SOFTWARE_2022,
    )
}


def get_citation(key: str) -> Citation:
    """Return one canonical bibliographic record.

    Raises
    ------
    KeyError
        If ``key`` is not registered.
    """
    return CITATIONS[key]


__all__ = [
    "ANDERSON_1995",
    "ANDERSON_MASUDA_ISAAK_1995",
    "BARRON_KLEIN_1965",
    "CITATIONS",
    "ELASTICITY_ULIAN_ET_AL_2018",
    "ELATE_GAILLAC_ET_AL_2016",
    "EOSFIT7_ANGEL_ET_AL_2014",
    "EOS_ULIAN_ET_AL_2014",
    "ERBA_2014",
    "ERBA_SHAHROKHI_MORADIAN_DOVESI_2015",
    "HILL_1952",
    "JAEKEN_COTTENIER_2016",
    "MCQUARRIE_SIMON_1997",
    "NYE_1985",
    "QHA_ULIAN_VALDRE_2018",
    "QUANTAS_2022",
    "SEISMIC_ULIAN_VALDRE_2024",
    "STIXRUDE_LITHGOW_BERTELLONI_2005",
    "WALLACE_1972",
    "WATERS_BIELAWSKI_2016",
    "MOUHAT_COUDERT_2014",
    "DESTEFANIS_RAVOUX_COSSARD_ERBA_2019",
    "DAVIES_1974",
    "PREM_DZIEWONSKI_ANDERSON_1981",
    "HASTEROK_CHAPMAN_2011",
    "PARSONS_SCLATER_1977",
    "KATSURA_2022",
    "KATSURA_SOFTWARE_2022",
    "get_citation",
]
