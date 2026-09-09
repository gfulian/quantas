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
    assert len(parser.total_energies()) == 1
    assert dft[0].kind is EnergyKind.DFT
    assert reference[0].kind is EnergyKind.REFERENCE
    assert reference[0].value == pytest.approx(-37937.46542432)
    assert parser.total_energies()[0].value == pytest.approx(dft[0].value)


@pytest.mark.interfaces
@pytest.mark.parametrize(
    ("correction_lines", "expected_total", "expected_marker", "corrections"),
    [
        ((), -275.28465266097, "TOTAL ENERGY(DFT)(AU)", ()),
        (
            (
                " GRIMME DISPERSION ENERGY (AU) -1.6299623064495E-02",
                " TOTAL ENERGY + DISP (AU) -2.7530095228403E+02",
            ),
            -275.30095228403,
            "TOTAL ENERGY + DISP (AU)",
            ("DISP",),
        ),
        (
            (
                " D3 DISPERSION ENERGY (AU)       -2.1575387300529E-02",
                " TOTAL ENERGY + DISP (AU)        -2.7530622804827E+02",
            ),
            -275.30622804827,
            "TOTAL ENERGY + DISP (AU)",
            ("DISP",),
        ),
        (
            (
                " GCP ENERGY (AU)                  1.1981901029943E-02",
                " TOTAL ENERGY + GCP (AU)         -2.7527267075994E+02",
            ),
            -275.27267075994,
            "TOTAL ENERGY + GCP (AU)",
            ("GCP",),
        ),
        (
            (
                " D3 DISPERSION ENERGY (AU)      -9.0912767684237E-03",
                " GCP ENERGY (AU)                 7.6500852104158E-03",
                " TOTAL ENERGY + DISP + GCP (AU) -2.7513674677358E+02",
            ),
            -275.13674677358,
            "TOTAL ENERGY + DISP + GCP (AU)",
            ("DISP", "GCP"),
        ),
    ],
)
def test_crystal_parser_resolves_scf_and_total_energy(
    correction_lines: tuple[str, ...],
    expected_total: float,
    expected_marker: str,
    corrections: tuple[str, ...],
) -> None:
    """CRYSTAL total energy should include every printed a-posteriori correction."""
    scf_energy = -275.13530558202 if len(corrections) == 2 else -275.28465266097
    cycles = 6 if len(corrections) == 2 else 7
    parser = CrystalOutputParser(
        [
            (
                " == SCF ENDED - CONVERGENCE ON ENERGY      "
                f"E(AU) {scf_energy:.14E} CYCLES   {cycles}"
            ),
            (
                f" TOTAL ENERGY(DFT)(AU)({cycles:3d}) "
                f"{scf_energy:.14E} DE-3.8E-10 tester 3.5E-12"
            ),
            *correction_lines,
        ]
    )

    scf = parser.scf_energies()
    total = parser.total_energies()

    assert len(scf) == 1
    assert len(total) == 1
    assert scf[0].value == pytest.approx(scf_energy)
    assert total[0].kind is EnergyKind.TOTAL
    assert total[0].value == pytest.approx(expected_total)
    assert total[0].metadata["source_marker"] == expected_marker
    assert total[0].metadata["corrections"] == corrections
    assert total[0].metadata["scf_energy"] == pytest.approx(scf_energy)
    assert total[0].metadata["total_correction_energy"] == pytest.approx(
        expected_total - scf_energy
    )


@pytest.mark.interfaces
def test_crystal_parser_does_not_cross_energy_state_boundaries() -> None:
    """A corrected total must remain associated with its own SCF state."""
    parser = CrystalOutputParser(
        [
            " == SCF ENDED - CONVERGENCE ON ENERGY E(AU) -1.000E+01 CYCLES 3",
            " TOTAL ENERGY(DFT)(AU)(  3) -1.000E+01 DE-1E-8",
            " TOTAL ENERGY + DISP (AU) -1.010E+01",
            " == SCF ENDED - CONVERGENCE ON ENERGY E(AU) -2.000E+01 CYCLES 4",
            " TOTAL ENERGY(DFT)(AU)(  4) -2.000E+01 DE-1E-8",
        ]
    )

    totals = parser.total_energies()

    assert [record.value for record in totals] == pytest.approx([-10.1, -20.0])
    assert totals[0].metadata["corrections"] == ("DISP",)
    assert totals[1].metadata["corrections"] == ()


@pytest.mark.interfaces
def test_crystal_parser_prefers_most_complete_corrected_total() -> None:
    """A combined corrected total should supersede a partial one in one state."""
    parser = CrystalOutputParser(
        [
            " == SCF ENDED - CONVERGENCE ON ENERGY E(AU) -1.000E+01 CYCLES 3",
            " TOTAL ENERGY(DFT)(AU)(  3) -1.000E+01 DE-1E-8",
            " TOTAL ENERGY + DISP (AU) -1.010E+01",
            " TOTAL ENERGY + DISP + GCP (AU) -1.009E+01",
        ]
    )

    total = parser.total_energies()[0]

    assert total.value == pytest.approx(-10.09)
    assert total.metadata["source_marker"] == "TOTAL ENERGY + DISP + GCP (AU)"
    assert total.metadata["corrections"] == ("DISP", "GCP")


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
    assert len(parser.total_energies()) == 34
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
