"""Characterization tests for generic CRYSTAL output parsing."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from quantas.interfaces.crystal.geometry import CrystalGeometryParser
from quantas.interfaces.crystal.output import CrystalOutputParser
from quantas.models.computation import (
    EnergyKind,
    OptimizationStatus,
    RunTerminationStatus,
    SCFStatus,
)


ROOT = Path(__file__).resolve().parents[2]
PHONON_OUTPUT = (
    ROOT / "examples/qha/crystal-phonons/dol_pbe0_crystal_pdisp_p01.out"
)
QHA_OUTPUT = ROOT / "examples/qha/crystal-qha/mgo-b3lyp-crystal-qha.out"
ELASTIC_OUTPUT = (
    ROOT / "examples/thermoelasticity/crystal_outputs/dol_pbe0_soec_p+00.out"
)


@pytest.mark.interfaces
def test_crystal_parser_characterizes_real_phonon_output() -> None:
    """Generic parsing should preserve the known phonon-run facts."""
    parser = CrystalOutputParser(PHONON_OUTPUT)

    termination = parser.termination()
    assert termination.status is RunTerminationStatus.NORMAL

    scf = parser.scf_results()
    assert len(scf) == 1
    assert scf[0].status is SCFStatus.CONVERGED
    assert scf[0].cycles == 22
    assert scf[0].cycle_numbers[0] == 0
    assert scf[0].cycle_numbers[-1] == 21
    assert scf[0].final_energy is not None
    assert scf[0].final_energy.value == pytest.approx(-37937.465424317)

    dft = parser.dft_energies()
    reference = parser.reference_energies()
    assert len(dft) == 1
    assert len(reference) == 1
    assert dft[0].kind is EnergyKind.DFT
    assert reference[0].kind is EnergyKind.REFERENCE
    assert reference[0].value == pytest.approx(-37937.46542432)


@pytest.mark.interfaces
def test_crystal_parser_characterizes_native_qha_multiple_runs() -> None:
    """One monolithic QHA output may contain many independent SCF/OPT runs."""
    parser = CrystalOutputParser(QHA_OUTPUT)

    assert parser.termination().status is RunTerminationStatus.NORMAL

    scf = parser.scf_results()
    assert len(scf) == 34
    assert all(result.status is SCFStatus.CONVERGED for result in scf)

    optimizations = parser.optimizations()
    assert len(optimizations) == 11
    assert all(
        result.status is OptimizationStatus.CONVERGED for result in optimizations
    )
    assert all(result.cycles == 1 for result in optimizations)

    assert len(parser.dft_energies()) == 34
    assert len(parser.reference_energies()) == 11


@pytest.mark.interfaces
def test_crystal_parser_characterizes_real_elastic_optimization() -> None:
    """Optimization history should retain convergence metrics and final energy."""
    parser = CrystalOutputParser(ELASTIC_OUTPUT)

    assert parser.termination().status is RunTerminationStatus.NORMAL
    assert len(parser.scf_results()) == 2

    optimizations = parser.optimizations()
    assert len(optimizations) == 1
    result = optimizations[0]
    assert result.status is OptimizationStatus.CONVERGED
    assert result.cycles == 9
    assert result.final_energy is not None
    assert result.final_energy.kind is EnergyKind.TOTAL
    assert result.final_energy.value == pytest.approx(-1405.091273573)

    last = result.steps[-1]
    assert last.index == 9
    assert last.energy is not None
    assert last.energy.value == pytest.approx(-1405.0912735728)
    assert last.delta_energy == pytest.approx(-1.979e-6)
    assert last.max_gradient == pytest.approx(9.0e-6)
    assert last.rms_gradient == pytest.approx(6.0e-6)
    assert last.max_displacement == pytest.approx(3.4e-5)
    assert last.rms_displacement == pytest.approx(2.4e-5)


@pytest.mark.interfaces
def test_crystal_parser_keeps_incomplete_scf_history() -> None:
    """A truncated SCF block should retain its observed energy history."""
    parser = CrystalOutputParser(
        [
            " CYC   0 ETOT(AU) -1.000000000000E+01 DETOT -1.00E+01 tst 0 PX 1",
            " CYC   1 ETOT(AU) -1.100000000000E+01 DETOT -1.00E+00 tst 0 PX 1",
        ]
    )

    assert parser.termination().status is RunTerminationStatus.INCOMPLETE
    scf = parser.scf_results()
    assert len(scf) == 1
    assert scf[0].status is SCFStatus.INCOMPLETE
    np.testing.assert_allclose(scf[0].energies, [-10.0, -11.0])


@pytest.mark.interfaces
def test_crystal_parser_distinguishes_failed_optimization() -> None:
    """Optimization failure is a parsed fact independent of run termination."""
    parser = CrystalOutputParser(
        [
            " COORDINATE AND CELL OPTIMIZATION - POINT    1",
            " TOTAL ENERGY(DFT)(AU)(  5) -1.234500000000E+02 DE-1.0E-06 tester 1E-8",
            " MAX GRADIENT      2.000000E-04",
            " RMS GRADIENT      1.000000E-04",
            " * OPT END - FAILED * E(AU):  -1.234500000000E+02  POINTS    1 *",
            " EEEEEEEEEE TERMINATION  DATE 01 01 2026 TIME 00:00:00.0",
        ]
    )

    assert parser.termination().status is RunTerminationStatus.NORMAL
    result = parser.optimizations()[0]
    assert result.status is OptimizationStatus.FAILED
    assert result.final_energy is not None
    assert result.final_energy.value == pytest.approx(-123.45)


@pytest.mark.interfaces
def test_crystal_geometry_normalizes_conventional_atomic_numbers() -> None:
    """CRYSTAL conventional Z values must not leak into neutral structures."""
    parser = CrystalGeometryParser(
        [
            " PRIMITIVE CELL - CENTRING CODE 1/0 VOLUME= 125.0",
            " A B C ALPHA BETA GAMMA",
            " 5.0 5.0 5.0 90.0 90.0 90.0",
            " ATOMS IN THE ASYMMETRIC UNIT 1 - ATOMS IN THE UNIT CELL: 1",
            " ATOM X/A Y/B Z/C",
            " 1 T 208 O 0.0 0.0 0.0",
        ]
    )

    structure = parser.initial_primitive_cell()

    np.testing.assert_array_equal(structure.atomic_numbers, [8])
    np.testing.assert_array_equal(
        structure.metadata["crystal_conventional_atomic_numbers"],
        [208],
    )
