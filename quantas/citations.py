"""
This module contains some references that have be cited when using specific
features of Quantas.
"""


def biblio_header():
    msg = "{:_^80}\n".format("")
    msg += """
The methods employed for this Quantas run are described in the following list.
Please, make sure to provide adequate citation of these references if you are
going to publish the results obtained by Quantas.
Quantas is an academic, open-source software, and, to aid its development,
scientific credit in the community is of utmost importance.
Your help would be really appreciated, and we thank you for it!
"""
    return msg


def biblio_footer():
    msg = "{:_^80}\n".format("")
    return msg


def quantas_citation():
    msg = """
For the base use of Quantas, please cite:

  Gianfranco Ulian and Giovanni Valdre'
  'QUANTAS, a Python software for the analysis of solids from ab initio
  quantum mechanical simulations and experimental data'
  Journal of Applied Crystallography 55, 386-396 (2022)
  http://dx.doi.org/10.1107/S1600576722000085
"""
    return msg


def qha_citation():
    msg = """
For quasi-harmonic approximation calculations used in this run, please cite:

  Gianfranco Ulian and Giovanni Valdre'
  'Equation of state of hexagonal hydroxylapatite (P63) as obtained from
  density functional theory simulations.'
  International Journal of Quantum Chemistry, 15, e25553 (2018)
  https://doi.org/10.1002/qua.25553
"""
    return msg


def eos_citation():
    msg = """
For equation of state fitting used in this run, please cite:

  Gianfranco Ulian, Sergio Tosoni and Giovanni Valdrè
  'The compressional behaviour and the mechanical properties of talc
  [Mg3Si4O10(OH)2]: a density functional theory investigation.'
  Physics and Chemistry of Minerals 41, 639-650 (2014)
  https://doi.org/10.1007/s00269-014-0677-x
"""
    return msg


def soec_citation():
    msg = """
For the analysis of the second order elastic constants used in this run,
please cite:

  Gianfranco Ulian, Daniele Moro and Giovanni Valdre'
  'First principle investigation of the mechanical properties of natural
  mineral layered nanocomposite: clinochlore as a model system.'
  Composite Structures, 202, 551-558 (2018).
  https://doi.org/10.1016/j.compstruct.2018.02.089

  Romain Gaillac, Pluton Pullumbi and François-Xavier Coudert
  'ELATE: an open-source online application for analysis and visualization of
  elastic tensors.'
  Journal of Physics: Condensed Matter 28, 275201 (2016)
  https://doi.org/10.1088/0953-8984/28/27/275201
"""
    return msg


def christoffel_citation():
    msg = """
For the analysis of the wave velocities via solving the Christoffel equation,
please cite:

  Gianfranco Ulian, Daniele Moro and Giovanni Valdre'
  'First principle investigation of the mechanical properties of natural
  mineral layered nanocomposite: clinochlore as a model system.'
  Composite Structures, 202, 551-558 (2018).
  https://doi.org/10.1016/j.compstruct.2018.02.089

  Jan W. Jaeken, Stefaan Cottenier
  'Solving the Christoffel equation: Phase and group velocities.'
  Computer Physics Communications, 207, 445-51 (2016)
  https://doi.org/10.1016/j.cpc.2016.06.014

  Romain Gaillac, Pluton Pullumbi and François-Xavier Coudert
  'ELATE: an open-source online application for analysis and visualization of
  elastic tensors.'
  Journal of Physics: Condensed Matter 28, 275201 (2016)
  https://doi.org/10.1088/0953-8984/28/27/275201
"""
    return msg
