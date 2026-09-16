# -*- coding: utf-8 -*-

"""Canonical Quantas bibliography registry."""

from __future__ import annotations

from collections.abc import Mapping
import re

from .models import Citation, CitationKind

QUANTAS_2022 = Citation(
    key="quantas_2022",
    authors=("G. Ulian", "G. Valdrè"),
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
    authors=("G. Ulian", "G. Valdrè"),
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
    authors=("D. A. McQuarrie", "J. D. Simon"),
    title="Physical Chemistry: A Molecular Approach",
    year=1997,
    kind=CitationKind.BOOK,
    publisher="University Science Books, Sausalito, California",
)

KIEFFER_1979 = Citation(
    key="kieffer_1979",
    authors=("S. W. Kieffer",),
    title=(
        "Thermodynamics and lattice vibrations of minerals: 3. Lattice "
        "dynamics and an approximation for minerals with application to "
        "simple substances and framework silicates"
    ),
    year=1979,
    journal="Reviews of Geophysics and Space Physics",
    volume="17",
    pages="35-59",
    doi="10.1029/RG017i001p00035",
)

ANDERSON_1995 = Citation(
    key="anderson_1995",
    authors=("O. L. Anderson",),
    title="Equations of State of Solids for Geophysics and Ceramic Science",
    year=1995,
    kind=CitationKind.BOOK,
    volume="31",
    journal="Oxford Monographs on Geology and Geophysics",
    publisher="Oxford University Press, New York",
)

ANDERSON_MASUDA_ISAAK_1995 = Citation(
    key="anderson_masuda_isaak_1995",
    authors=("O. L. Anderson", "K. Masuda", "D. G. Isaak"),
    title="A new thermodynamic approach for high-pressure physics",
    year=1995,
    journal="Physics of the Earth and Planetary Interiors",
    volume="91",
    pages="3-16",
    doi="10.1016/0031-9201(95)03044-W",
)

ERBA_2014 = Citation(
    key="erba_2014",
    authors=("A. Erba",),
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
        "A. Erba",
        "M. Shahrokhi",
        "R. Moradian",
        "R. Dovesi",
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
    authors=("G. Ulian", "S. Tosoni", "G. Valdrè"),
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
    authors=("R. J. Angel", "J. Gonzalez-Platas", "M. Alvaro"),
    title="EosFit7c and a Fortran module (library) for equation of state calculations",
    year=2014,
    journal="Zeitschrift für Kristallographie",
    volume="229",
    pages="405-419",
    doi="10.1515/zkri-2013-1711",
)

SJEOS_ALCHAGIROV_ET_AL_2001 = Citation(
    key="sjeos_alchagirov_perdew_boettger_albers_fiolhais_2001",
    authors=(
        "A. B. Alchagirov",
        "J. P. Perdew",
        "J. C. Boettger",
        "R. C. Albers",
        "C. Fiolhais",
    ),
    title=(
        "Energy and pressure versus volume: Equations of state motivated by "
        "the stabilized jellium model"
    ),
    year=2001,
    journal="Physical Review B",
    volume="63",
    pages="224115",
    doi="10.1103/PhysRevB.63.224115",
)


STAROVEROV_SCUSERIA_TAO_PERDEW_2004 = Citation(
    key="staroverov_scuseria_tao_perdew_2004",
    authors=(
        "V. N. Staroverov",
        "G. E. Scuseria",
        "J. Tao",
        "J. P. Perdew",
    ),
    title="Tests of a ladder of density functionals for bulk solids and surfaces",
    year=2004,
    journal="Physical Review B",
    volume="69",
    pages="075102",
    doi="10.1103/PhysRevB.69.075102",
)


NYE_1985 = Citation(
    key="nye_1985",
    authors=("J. F. Nye",),
    title=(
        "Physical Properties of Crystals: Their Representation by Tensors and "
        "Matrices"
    ),
    year=1985,
    kind=CitationKind.BOOK,
    publisher="Oxford University Press, Oxford, 2nd edition",
)

HILL_1952 = Citation(
    key="hill_1952",
    authors=("R. Hill",),
    title="The elastic behaviour of a crystalline aggregate",
    year=1952,
    journal="Proceedings of the Physical Society. Section A",
    volume="65",
    pages="349-354",
    doi="10.1088/0370-1298/65/5/307",
)

ELASTICITY_ULIAN_ET_AL_2018 = Citation(
    key="elasticity_ulian_moro_valdre_2018",
    authors=("G. Ulian", "D. Moro", "G. Valdrè"),
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
    authors=("R. Gaillac", "P. Pullumbi", "F.-X. Coudert"),
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
    authors=("G. Ulian", "G. Valdrè"),
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
    authors=("J. W. Jaeken", "S. Cottenier"),
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
    authors=("L. Stixrude", "C. Lithgow-Bertelloni"),
    title="Thermodynamics of mantle minerals—I. Physical properties",
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
    authors=("F. Mouhat", "F.-X. Coudert"),
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
        "M. Destefanis",
        "C. Ravoux",
        "A. Cossard",
        "A. Erba",
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
    authors=("A. M. Dziewonski", "D. L. Anderson"),
    title="Preliminary reference Earth model",
    year=1981,
    journal="Physics of the Earth and Planetary Interiors",
    volume="25",
    pages="297-356",
    doi="10.1016/0031-9201(81)90046-7",
)

HASTEROK_CHAPMAN_2011 = Citation(
    key="hasterok_chapman_2011",
    authors=("D. Hasterok", "D. S. Chapman"),
    title="Heat production and geotherms for the continental lithosphere",
    year=2011,
    journal="Earth and Planetary Science Letters",
    volume="307",
    pages="59-70",
    doi="10.1016/j.epsl.2011.04.034",
)

PARSONS_SCLATER_1977 = Citation(
    key="parsons_sclater_1977",
    authors=("B. Parsons", "J. G. Sclater"),
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
    authors=("T. Katsura",),
    title="A revised adiabatic temperature profile for the mantle",
    year=2022,
    journal="Journal of Geophysical Research: Solid Earth",
    volume="127",
    pages="e2021JB023562",
    doi="10.1029/2021JB023562",
)

KATSURA_SOFTWARE_2022 = Citation(
    key="katsura_software_2022",
    authors=("T. Katsura",),
    title=(
        "Matlab scripts of 'A revised adiabatic temperature profile for the "
        "mantle' (Version 1.1.0)"
    ),
    year=2022,
    kind=CitationKind.SOFTWARE,
    publisher="Zenodo",
    doi="10.5281/zenodo.5903286",
)

OREAR_1982 = Citation(
    key="orear_1982",
    authors=("J. Orear",),
    title="Least squares when both variables have uncertainties",
    year=1982,
    journal="American Journal of Physics",
    volume="50",
    pages="912-916",
    doi="10.1119/1.12972",
)

BOGGS_BYRD_ROGERS_SCHNABEL_1992 = Citation(
    key="boggs_byrd_rogers_schnabel_1992",
    authors=("P. T. Boggs", "R. H. Byrd", "J. E. Rogers", "R. B. Schnabel"),
    title=(
        "User's Reference Guide for ODRPACK Version 2.01: Software for "
        "Weighted Orthogonal Distance Regression"
    ),
    year=1992,
    kind=CitationKind.REPORT,
    publisher="National Institute of Standards and Technology, Gaithersburg, MD",
    report_number="NISTIR 4834",
    doi="10.6028/NIST.IR.4834",
)

ZWOLAK_BOGGS_WATSON_2007 = Citation(
    key="zwolak_boggs_watson_2007",
    authors=("J. W. Zwolak", "P. T. Boggs", "L. T. Watson"),
    title=(
        "Algorithm 869: ODRPACK95: A weighted orthogonal distance regression "
        "code with bound constraints"
    ),
    year=2007,
    journal="ACM Transactions on Mathematical Software",
    volume="33",
    pages="Article 27",
    doi="10.1145/1268776.1268782",
)

BERMAN_1988 = Citation(
    key="berman_1988",
    authors=("R. G. Berman",),
    title=(
        "Internally-Consistent Thermodynamic Data for Minerals in the System "
        "Na2O-K2O-CaO-MgO-FeO-Fe2O3-Al2O3-SiO2-TiO2-H2O-CO2"
    ),
    year=1988,
    journal="Journal of Petrology",
    volume="29",
    pages="445-522",
    doi="10.1093/petrology/29.2.445",
)

FEI_1995 = Citation(
    key="fei_1995",
    authors=("Y. Fei",),
    title="Thermal Expansion",
    year=1995,
    kind=CitationKind.CHAPTER,
    container_title=(
        "Mineral Physics & Crystallography: A Handbook of Physical Constants"
    ),
    editors=("T. J. Ahrens",),
    volume="2",
    pages="29-44",
    publisher="American Geophysical Union, Washington, DC",
    doi="10.1029/RF002p0029",
)

PAWLEY_REDFERN_HOLLAND_1996 = Citation(
    key="pawley_redfern_holland_1996",
    authors=("A. R. Pawley", "S. A. T. Redfern", "T. J. B. Holland"),
    title=(
        "Volume behavior of hydrous minerals at high pressure and temperature: "
        "I. Thermal expansion of lawsonite, zoisite, clinozoisite, and diaspore"
    ),
    year=1996,
    journal="American Mineralogist",
    volume="81",
    pages="335-340",
    doi="10.2138/am-1996-3-407",
)

SALJE_WRUCK_THOMAS_1991 = Citation(
    key="salje_wruck_thomas_1991",
    authors=("E. K. H. Salje", "B. Wruck", "H. Thomas"),
    title="Order-parameter saturation and low-temperature extension of Landau theory",
    year=1991,
    journal="Zeitschrift für Physik B Condensed Matter",
    volume="82",
    pages="399-404",
    doi="10.1007/BF01357186",
)

HOLLAND_POWELL_2011 = Citation(
    key="holland_powell_2011",
    authors=("T. J. B. Holland", "R. Powell"),
    title=(
        "An improved and extended internally consistent thermodynamic dataset "
        "for phases of petrological interest, involving a new equation of state "
        "for solids"
    ),
    year=2011,
    journal="Journal of Metamorphic Geology",
    volume="29",
    pages="333-383",
    doi="10.1111/j.1525-1314.2010.00923.x",
)

HELFFRICH_CONNOLLY_2009 = Citation(
    key="helffrich_connolly_2009",
    authors=("G. R. Helffrich", "J. A. D. Connolly"),
    title=(
        "Physical contradictions and remedies using simple polythermal equations "
        "of state"
    ),
    year=2009,
    journal="American Mineralogist",
    volume="94",
    pages="1616-1619",
    doi="10.2138/am.2009.3262",
)

KROLL_KIRFEL_HEINEMANN_BARBIER_2012 = Citation(
    key="kroll_kirfel_heinemann_barbier_2012",
    authors=("H. Kroll", "A. Kirfel", "R. Heinemann", "B. Barbier"),
    title=(
        "Volume thermal expansion and related thermophysical parameters in the "
        "Mg,Fe olivine solid-solution series"
    ),
    year=2012,
    journal="European Journal of Mineralogy",
    volume="24",
    pages="935-956",
    doi="10.1127/0935-1221/2012/0024-2235",
)

ERBA_MAHMOUD_BELMONTE_DOVESI_2014 = Citation(
    key="erba_mahmoud_belmonte_dovesi_2014",
    authors=("A. Erba", "A. Mahmoud", "D. Belmonte", "R. Dovesi"),
    title=(
        "High pressure elastic properties of minerals from ab initio simulations: "
        "The case of pyrope, grossular and andradite silicate garnets"
    ),
    year=2014,
    journal="Journal of Chemical Physics",
    volume="140",
    pages="124703",
    doi="10.1063/1.4869144",
)

KARKI_STIXRUDE_CLARK_WARREN_ACKLAND_CRAIN_1997 = Citation(
    key="karki_stixrude_clark_warren_ackland_crain_1997",
    authors=(
        "B. B. Karki",
        "L. Stixrude",
        "S. J. Clark",
        "M. C. Warren",
        "G. J. Ackland",
        "J. Crain",
    ),
    title="Structure and elasticity of MgO at high pressure",
    year=1997,
    journal="American Mineralogist",
    volume="82",
    pages="51-60",
    doi="10.2138/am-1997-1-207",
)

SINOGEIKIN_BASS_1999 = Citation(
    key="sinogeikin_bass_1999",
    authors=("S. V. Sinogeikin", "J. D. Bass"),
    title="Single-crystal elasticity of MgO at high pressure",
    year=1999,
    journal="Physical Review B",
    volume="59",
    pages="R14141-R14144",
    doi="10.1103/PhysRevB.59.R14141",
)

JIANG_SPEZIALE_DUFFY_2006 = Citation(
    key="jiang_speziale_duffy_2006",
    authors=("F. Jiang", "S. Speziale", "T. S. Duffy"),
    title=(
        "Elasticity of magnesite and dolomite from a genetic algorithm for "
        "inverting Brillouin spectroscopy measurements"
    ),
    year=2006,
    journal="Physics of the Earth and Planetary Interiors",
    volume="155",
    pages="1-20",
    doi="10.1016/j.pepi.2005.08.004",
)

CITATIONS = {
    citation.key: citation
    for citation in (
        QUANTAS_2022,
        QHA_ULIAN_VALDRE_2018,
        MCQUARRIE_SIMON_1997,
        KIEFFER_1979,
        ANDERSON_1995,
        ANDERSON_MASUDA_ISAAK_1995,
        ERBA_2014,
        ERBA_SHAHROKHI_MORADIAN_DOVESI_2015,
        EOS_ULIAN_ET_AL_2014,
        EOSFIT7_ANGEL_ET_AL_2014,
        SJEOS_ALCHAGIROV_ET_AL_2001,
        STAROVEROV_SCUSERIA_TAO_PERDEW_2004,
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
        OREAR_1982,
        BOGGS_BYRD_ROGERS_SCHNABEL_1992,
        ZWOLAK_BOGGS_WATSON_2007,
        BERMAN_1988,
        FEI_1995,
        PAWLEY_REDFERN_HOLLAND_1996,
        SALJE_WRUCK_THOMAS_1991,
        HOLLAND_POWELL_2011,
        HELFFRICH_CONNOLLY_2009,
        KROLL_KIRFEL_HEINEMANN_BARBIER_2012,
        ERBA_MAHMOUD_BELMONTE_DOVESI_2014,
        KARKI_STIXRUDE_CLARK_WARREN_ACKLAND_CRAIN_1997,
        SINOGEIKIN_BASS_1999,
        JIANG_SPEZIALE_DUFFY_2006,
    )
}


_KEY_PATTERN = re.compile(r"^[a-z0-9]+(?:_[a-z0-9]+)*$")
_AUTHOR_PATTERN = re.compile(
    r"^(?:[^\W\d_]\.(?:-[^\W\d_]\.)?\s+)+(?:\S.*)$"
)


def validate_citation_registry(
    citations: Mapping[str, Citation] = CITATIONS,
) -> None:
    """Validate canonical citation identifiers and bibliographic metadata.

    Parameters
    ----------
    citations : mapping of str to Citation, optional
        Citation mapping to validate. Defaults to the canonical Quantas
        registry.

    Raises
    ------
    ValueError
        If a citation key, author list, DOI, record-specific field, or DOI
        uniqueness constraint is invalid.
    """
    seen_doi: dict[str, str] = {}
    for key, citation in citations.items():
        if key != citation.key:
            raise ValueError(
                f"Citation mapping key {key!r} does not match record key "
                f"{citation.key!r}"
            )
        if _KEY_PATTERN.fullmatch(key) is None:
            raise ValueError(f"Invalid citation key: {key!r}")
        if not citation.authors:
            raise ValueError(f"Citation {key!r} must define at least one author")
        for author in citation.authors:
            if _AUTHOR_PATTERN.fullmatch(author) is None:
                raise ValueError(
                    f"Citation {key!r} author {author!r} must use initials + surname"
                )
        for editor in citation.editors:
            if _AUTHOR_PATTERN.fullmatch(editor) is None:
                raise ValueError(
                    f"Citation {key!r} editor {editor!r} must use initials + surname"
                )
        if not citation.title.strip():
            raise ValueError(f"Citation {key!r} must define a title")
        if citation.year <= 0:
            raise ValueError(f"Citation {key!r} has an invalid year")
        if citation.doi:
            if citation.doi.startswith(("http://", "https://", "doi:")):
                raise ValueError(
                    f"Citation {key!r} DOI must not include a URL or DOI prefix"
                )
            if citation.doi in seen_doi:
                raise ValueError(
                    f"Citation {key!r} duplicates DOI from {seen_doi[citation.doi]!r}"
                )
            seen_doi[citation.doi] = key
            if citation.url is not None:
                raise ValueError(
                    f"Citation {key!r} must not define url when a DOI is available"
                )
        if citation.kind in {CitationKind.ARTICLE, CitationKind.PREPRINT}:
            if not citation.journal:
                raise ValueError(f"Citation {key!r} must define a journal")
        elif citation.kind is CitationKind.CHAPTER:
            if not citation.container_title or not citation.publisher:
                raise ValueError(
                    f"Chapter citation {key!r} must define container_title "
                    "and publisher"
                )
        elif citation.kind is CitationKind.REPORT:
            if not citation.report_number or not citation.publisher:
                raise ValueError(
                    f"Report citation {key!r} must define report_number and publisher"
                )
        elif citation.kind in {CitationKind.BOOK, CitationKind.SOFTWARE}:
            if not citation.publisher:
                raise ValueError(f"Citation {key!r} must define a publisher")


validate_citation_registry()


def get_citation(key: str) -> Citation:
    """Return one canonical bibliographic record.

    Raises
    ------
    KeyError
        If ``key`` is not registered.
    """
    return CITATIONS[key]


__all__ = [
    "BERMAN_1988",
    "BOGGS_BYRD_ROGERS_SCHNABEL_1992",
    "ANDERSON_1995",
    "ANDERSON_MASUDA_ISAAK_1995",
    "BARRON_KLEIN_1965",
    "CITATIONS",
    "ELASTICITY_ULIAN_ET_AL_2018",
    "ELATE_GAILLAC_ET_AL_2016",
    "EOSFIT7_ANGEL_ET_AL_2014",
    "SJEOS_ALCHAGIROV_ET_AL_2001",
    "STAROVEROV_SCUSERIA_TAO_PERDEW_2004",
    "EOS_ULIAN_ET_AL_2014",
    "ERBA_2014",
    "ERBA_SHAHROKHI_MORADIAN_DOVESI_2015",
    "HILL_1952",
    "JAEKEN_COTTENIER_2016",
    "KIEFFER_1979",
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
    "ERBA_MAHMOUD_BELMONTE_DOVESI_2014",
    "FEI_1995",
    "HELFFRICH_CONNOLLY_2009",
    "HOLLAND_POWELL_2011",
    "JIANG_SPEZIALE_DUFFY_2006",
    "KARKI_STIXRUDE_CLARK_WARREN_ACKLAND_CRAIN_1997",
    "KROLL_KIRFEL_HEINEMANN_BARBIER_2012",
    "OREAR_1982",
    "PAWLEY_REDFERN_HOLLAND_1996",
    "SALJE_WRUCK_THOMAS_1991",
    "SINOGEIKIN_BASS_1999",
    "ZWOLAK_BOGGS_WATSON_2007",
    "get_citation",
    "validate_citation_registry",
]
