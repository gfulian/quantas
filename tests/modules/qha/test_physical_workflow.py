from __future__ import annotations

import numpy as np

from quantas.modules.qha.calculator import QHACalculator
from quantas.cli.qha_observer import QHATextObserver
from quantas.cli.qha_render import render_phonon_frequency_fit_report
from quantas.modules.qha.models import QHAInput, QHAOptions


def make_physical_input() -> QHAInput:
    volume = np.linspace(9.5, 10.8, 7)
    energy = 2.0e-3 * (volume - 10.0) ** 2
    base = 450.0 * (10.0 / volume) ** 2
    frequencies = np.stack([base, 1.15 * base, 1.4 * base], axis=0)[None, :, :]
    return QHAInput(
        jobname="physical-workflow",
        natoms=1,
        qpoints=1,
        volume=volume,
        energy=energy,
        frequencies=frequencies,
        weights=np.array([1.0]),
    )


def make_options(**kwargs: object) -> QHAOptions:
    data = dict(
        temperature_min=10.0,
        temperature_max=1000.0,
        temperature_step=990.0,
        pressure_min=0.0,
        pressure_max=0.0,
        pressure_step=1.0,
        scheme="td",
        minimization="poly",
        free_energy_degree=3,
        debug=True,
    )
    data.update(kwargs)
    return QHAOptions(**data)


def test_qha_minimization_uses_temperature_dependent_free_energy() -> None:
    result = (
        QHACalculator(make_physical_input(), make_options()).execute().results["qha"]
    )

    assert result.free_energy is not None
    assert result.equilibrium_volume is not None
    assert result.equilibrium_volume.shape == (2, 1)
    assert result.equilibrium_volume[1, 0] > result.equilibrium_volume[0, 0]


def test_phonon_fit_report_does_not_require_q_coordinates() -> None:
    text = render_phonon_frequency_fit_report(
        make_physical_input(),
        make_options(scheme="freq"),
        include_debug=False,
    )

    assert "Fitted records" in text
    assert "Skipped zero records" in text
    assert " 0" not in text.split("Fitted records", 1)[1].split("\n", 1)[0]


def test_qha_text_observer_does_not_trim_debug_records_by_default() -> None:
    observer = QHATextObserver(silent=True, verbosity="debug")

    assert observer.max_debug_records is None
